/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest.h>
#include <t8_data/t8_element_array.h>
#include <t8_data/t8_locidx_list.h>
#include <t8_element_cxx.hxx>

/* We want to export the whole implementation to be callable from "C" */
T8_EXTERN_C_BEGIN ();

/* Check the lastly inserted elements of an array for recursive coarsening.
 * The last inserted element must be the last element of a family.
 * \param [in] forest  The new forest currently in construction.
 * \param [in] ltreeid The current local tree.
 * \param [in] lelement_id The id of the currently coarsened element in the tree of the original forest.
 * \param [in] ts      The scheme for this local tree.
 * \param [in] telements The array of newly created (adapted) elements.
 *                      The last inserted element must be the last child in its family.
 * \param [in] el_coarsen the index of the first element in \a telement
 *                        array which could be coarsened recursively.
 * \param [in,out] el_inserted On input the number of elements in \a telement, on output
 *                        the new number of elements (so it will be smaller or equal to its input).
 * \param [in] el_buffer Buffer space to store a family of elements.
 */
static void
t8_forest_adapt_coarsen_recursive (t8_forest_t forest, t8_locidx_t ltreeid,
                                   t8_locidx_t lelement_id,
                                   t8_eclass_scheme_c * ts,
                                   t8_element_array_t * telements,
                                   t8_locidx_t el_coarsen,
                                   t8_locidx_t * el_inserted,
                                   t8_element_t ** el_buffer)
{
  t8_element_t       *element;
  t8_locidx_t         pos;
  size_t              elements_in_array;
  int                 num_children, i, isfamily;
  /* el_inserted is the index of the last element in telements plus one.
   * el_coarsen is the index of the first element which could possibly
   * be coarsened. */

  elements_in_array = t8_element_array_get_count (telements);
  T8_ASSERT (*el_inserted == (t8_locidx_t) elements_in_array);
  T8_ASSERT (el_coarsen >= 0);
  element = t8_element_array_index_locidx (telements, *el_inserted - 1);
  /* TODO: This assumes that the number of children is the same for each
   *       element in that class. This may not be the case. */
  num_children = ts->t8_element_num_children (element);
  T8_ASSERT (ts->t8_element_child_id (element) == num_children - 1);

  t8_element_t      **fam = el_buffer;
  const t8_element_t **constfam = (const t8_element_t **) fam;
  pos = *el_inserted - num_children;
  isfamily = 1;
  while (isfamily && pos >= el_coarsen && ts->t8_element_child_id (element)
         == num_children - 1) {
    isfamily = 1;
    /* Get all elements at indices pos, pos + 1, ... ,pos + num_children - 1 */
    for (i = 0; i < num_children; i++) {
      fam[i] = t8_element_array_index_locidx (telements, pos + i);
      if (ts->t8_element_child_id (fam[i]) != i) {
        /* These elements cannot form a family. Stop coarsening. */
        isfamily = 0;
        break;
      }
    }
    T8_ASSERT (!isfamily || ts->t8_element_is_family (constfam));
    if (isfamily
        && forest->set_adapt_fn (forest, forest->set_from, ltreeid,
                                 lelement_id, ts, num_children,
                                 constfam) < 0) {
      /* Coarsen the element */
      *el_inserted -= num_children - 1;
      /* remove num_children - 1 elements from the array */
      T8_ASSERT (elements_in_array == t8_element_array_get_count (telements));
      ts->t8_element_parent (constfam[0], fam[0]);
      elements_in_array -= num_children - 1;
      t8_element_array_resize (telements, elements_in_array);
      /* Set element to the new constructed parent. Since resizing the array
       * may change the position in memory, we have to do it after resizing. */
      element = t8_element_array_index_locidx (telements, pos);
    }
    else {
      /* If the elements are no family or
       * the family is not to be coarsened we abort the coarsening process */
      isfamily = 0;
    }
    pos -= num_children - 1;
  }
}

/* Check the lastly inserted element of an array for recursive refining.
 * \param [in] forest  The new forest currently in construction.
 * \param [in] ltreeid The current local tree.
 * \param [in] lelement_id The id of the currently coarsened element in the tree of the original forest.
 * \param [in] ts      The scheme for this local tree.
 * \param [in] elem_list Helper list to temporarily insert the newly refined elements.
 *                       These will eventually get copied to \a telements.
 * \param [in] telements The array of newly created (adapted) elements.
 *                      The last inserted element must be the last child in its family.
 * \param [in,out] num_inserted On input the number of elements in \a telement, on output
 *                        the new number of elements (so it will be smaller or equal to its input).
 * \param [in] el_buffer Enough buffer space to store all children of the lastly created element.
 */
static void
t8_forest_adapt_refine_recursive (t8_forest_t forest, t8_locidx_t ltreeid,
                                  t8_locidx_t lelement_id,
                                  t8_eclass_scheme_c * ts,
                                  sc_list_t * elem_list,
                                  t8_element_array_t * telements,
                                  t8_locidx_t * num_inserted,
                                  t8_element_t ** el_buffer)
{
  t8_element_t       *insert_el;
  int                 num_children;
  int                 ci;

  if (elem_list->elem_count <= 0) {
    return;
  }
  while (elem_list->elem_count > 0) {
    /* Until the list is empty we
     * - remove the first element from the list.
     * - Check whether it should get refined
     * - If yes, we add all its children to the list
     * - If no, we add the element to the array of new elements
     */
    el_buffer[0] = (t8_element_t *) sc_list_pop (elem_list);
    num_children = ts->t8_element_num_children (el_buffer[0]);
    if (forest->set_adapt_fn (forest, forest->set_from, ltreeid, lelement_id,
                              ts, 1, (const t8_element_t **) el_buffer) > 0) {
      /* The element should be refined */
      if (ts->t8_element_level (el_buffer[0]) < forest->maxlevel) {
        /* only refine, if we do not exceed the maximum allowed level */
        /* Create the children and add them to the list */
        ts->t8_element_new (num_children - 1, el_buffer + 1);
        ts->t8_element_children (el_buffer[0], num_children, el_buffer);
        for (ci = num_children - 1; ci >= 0; ci--) {
          (void) sc_list_prepend (elem_list, el_buffer[ci]);
        }
      }
    }
    else {
      /* This element should not get refined,
       * we remove it from the buffer and add it to the array of new elements. */
      insert_el = t8_element_array_push (telements);
      ts->t8_element_copy (el_buffer[0], insert_el);
      ts->t8_element_destroy (1, el_buffer);
      (*num_inserted)++;
    }
  }
}

/* TODO: optimize this when we own forest_from */
void
t8_forest_adapt (t8_forest_t forest)
{
  t8_forest_t         forest_from;
  sc_list_t          *refine_list = NULL;       /* This is only needed when we adapt recursively */
  t8_element_array_t *telements, *telements_from;
  t8_locidx_t         ltree_id, num_trees;
  t8_locidx_t         el_considered;
  t8_locidx_t         el_inserted;
  t8_locidx_t         el_coarsen;
  t8_locidx_t         num_el_from;
  t8_locidx_t         el_offset;
  size_t              num_children, zz;
  t8_tree_t           tree, tree_from;
  t8_eclass_scheme_c *tscheme;
  t8_element_t      **elements, **elements_from;
  int                 refine;
  int                 ci;
  int                 num_elements;
#ifdef T8_ENABLE_DEBUG
  int                 is_family;
#endif

  T8_ASSERT (forest != NULL);
  T8_ASSERT (forest->set_from != NULL);
  T8_ASSERT (forest->set_adapt_recursive == 0
             || forest->set_adapt_recursive == 1);

  /* if profiling is enabled, measure runtime */
  if (forest->profile != NULL) {
    forest->profile->adapt_runtime = -sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occured on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("Start adapt %f %f\n", sc_MPI_Wtime (),
                           forest->profile->adapt_runtime);
  }

  forest_from = forest->set_from;
  t8_global_productionf ("Into t8_forest_adapt from %lld total elements\n",
                         (long long) forest_from->global_num_elements);

  /* TODO: Allocate memory for the trees of forest.
   * Will we do this here or in an extra function? */
  T8_ASSERT (forest->trees->elem_count == forest_from->trees->elem_count);

  if (forest->set_adapt_recursive) {
    refine_list = sc_list_new (NULL);
  }
  forest->local_num_elements = 0;
  el_offset = 0;
  num_trees = t8_forest_get_num_local_trees (forest);
  /* Iterate over the trees and build the new element arrays for each one. */
  for (ltree_id = 0; ltree_id < num_trees; ltree_id++) {
    /* Get the new and old tree and the new and old element arrays */
    tree = t8_forest_get_tree (forest, ltree_id);
    tree_from = t8_forest_get_tree (forest_from, ltree_id);
    telements = &tree->elements;
    telements_from = &tree_from->elements;
    /* Number of elements in the old tree */
    num_el_from = (t8_locidx_t) t8_element_array_get_count (telements_from);
    T8_ASSERT (num_el_from ==
               t8_forest_get_tree_num_elements (forest_from, ltree_id));
    /* Get the element scheme for this tree */
    tscheme = t8_forest_get_eclass_scheme (forest_from, tree->eclass);
    /* Index of the element we currently consider for refinement/coarsening. */
    el_considered = 0;
    /* Index into the newly inserted elements */
    el_inserted = 0;
    /* el_coarsen is the index of the first element in the new element
     * array which could be coarsened recursively. */
    el_coarsen = 0;
    /* TODO: this will generate problems with pyramidal elements */
    num_children =
      tscheme->t8_element_num_children (t8_element_array_index_locidx
                                        (telements_from, 0));
    /* Buffer for a family of new elements */
    elements = T8_ALLOC (t8_element_t *, num_children);
    /* Buffer for a family of old elements */
    elements_from = T8_ALLOC (t8_element_t *, num_children);
    /* We now iterate over all elements in this tree and check them for refinement/coarsening. */
    while (el_considered < num_el_from) {
#ifdef T8_ENABLE_DEBUG
      /* Will get set to 0 later if this is not a family */
      is_family = 1;
#endif
      /* Load the current element and at most num_children-1 many others into
       * the elements_from buffer. Stop when we are certain that they cannot from
       * a family.
       * At the end is_family will be true, if these elements form a family.
       */
      num_elements = num_children;
      for (zz = 0; zz < num_children &&
           el_considered + (t8_locidx_t) zz < num_el_from; zz++) {
        elements_from[zz] = t8_element_array_index_locidx (telements_from,
                                                           el_considered +
                                                           zz);
        if ((size_t) tscheme->t8_element_child_id (elements_from[zz]) != zz) {
          break;
        }
      }
      if (zz != num_children) {
        /* We are certain that the elements do not form a family.
         * So we will only pass the first element to the adapt callback. */
        num_elements = 1;
#ifdef T8_ENABLE_DEBUG
        is_family = 0;
#endif
      }
      T8_ASSERT (!is_family
                 || tscheme->t8_element_is_family ((const t8_element_t **)
                                                   elements_from));
      /* Pass the element, or the family to the adapt callback.
       * The output will be > 0 if the element should be refined
       *                    = 0 if the element should remain as is
       *                    < 0 if we passed a family and it should get coarsened.
       */
      refine =
        forest->set_adapt_fn (forest, forest->set_from, ltree_id,
                              el_considered, tscheme, num_elements,
                              (const t8_element_t **) elements_from);
      T8_ASSERT (is_family || refine >= 0);
      if (refine > 0 && tscheme->t8_element_level (elements_from[0]) >=
          forest->maxlevel) {
        /* Only refine an element if it does not exceed the maximum level */
        refine = 0;
      }
      if (refine > 0) {
        /* The first element is to be refined */
        if (forest->set_adapt_recursive) {
          /* Create the children of this element */
          tscheme->t8_element_new (num_children, elements);
          tscheme->t8_element_children (elements_from[0], num_children,
                                        elements);
          for (ci = num_children - 1; ci >= 0; ci--) {
            /* Prepend the children to the refine_list.
             * These should now be the only elements in the list.
             */
            (void) sc_list_prepend (refine_list, elements[ci]);
          }
          /* We now recursively check the newly created elements for refinement. */
          t8_forest_adapt_refine_recursive (forest, ltree_id, el_considered,
                                            tscheme, refine_list, telements,
                                            &el_inserted, elements);
          /* el_coarsen is the index of the first element in the new element
           * array which could be coarsened recursively.
           * We can set this here to the next element after the current family, since a family that emerges from a refinement will never be coarsened */
          el_coarsen = el_inserted + num_children;
        }
        else {
          /* The element should be refined and refinement is not recursive. */
          /* We add the children to the element array of the current tree. */
          (void) t8_element_array_push_count (telements, num_children);
          for (zz = 0; zz < num_children; zz++) {
            elements[zz] =
              t8_element_array_index_locidx (telements, el_inserted + zz);
          }
          tscheme->t8_element_children (elements_from[0], num_children,
                                        elements);
          el_inserted += num_children;
        }
        el_considered++;
      }
      else if (refine < 0) {
        /* The elements form a family and are to be coarsened. */
        /* Make room for one more new element. */
        elements[0] = t8_element_array_push (telements);
        /* Compute the parent of the current family.
         * This parent is now inserted in telements. */
        tscheme->t8_element_parent (elements_from[0], elements[0]);
        el_inserted++;
        if (forest->set_adapt_recursive) {
          /* Adaptation is recursive.
           * We check whether the just generated parent is the last in its
           * family.
           * If so, we check this family for recursive coarsening. */
          if ((size_t) tscheme->t8_element_child_id (elements[0])
              == num_children - 1) {
            t8_forest_adapt_coarsen_recursive (forest, ltree_id,
                                               el_considered, tscheme,
                                               telements, el_coarsen,
                                               &el_inserted, elements);
          }
        }
        el_considered += num_children;
      }
      else {
        /* The considered elements are neither to be coarsened nor is the first
         * one to be refined.
         * We copy the element to the new element array. */
        T8_ASSERT (refine == 0);
        elements[0] = t8_element_array_push (telements);
        tscheme->t8_element_copy (elements_from[0], elements[0]);
        el_inserted++;
        if (forest->set_adapt_recursive &&
            (size_t) tscheme->t8_element_child_id (elements[0])
            == num_children - 1) {
          /* If adaptation is recursive and this was the last element in its
           * family, we need to check for recursive coarsening. */
          t8_forest_adapt_coarsen_recursive (forest, ltree_id, el_considered,
                                             tscheme, telements, el_coarsen,
                                             &el_inserted, elements);
        }
        el_considered++;
      }
    }
    /* Check that if we had recursive adaptation, the refine list is now empty. */
    T8_ASSERT (!forest->set_adapt_recursive || refine_list->elem_count == 0);
    /* Set the new element offset of this tree */
    tree->elements_offset = el_offset;
    el_offset += el_inserted;
    /* Add to the new number of local elements. */
    forest->local_num_elements += el_inserted;
    /* Possibly shrink the telements array to the correct size */
    t8_element_array_resize (telements, el_inserted);

    /* clean up */
    T8_FREE (elements);
    T8_FREE (elements_from);
  }
  if (forest->set_adapt_recursive) {
    /* clean up */
    sc_list_destroy (refine_list);
  }

  /* We now adapted all local trees */
  /* Compute the new global number of elements */
  t8_forest_comm_global_num_elements (forest);
  t8_global_productionf ("Done t8_forest_adapt with %lld total elements\n",
                         (long long) forest->global_num_elements);

  /* if profiling is enabled, measure runtime */
  if (forest->profile != NULL) {
    forest->profile->adapt_runtime += sc_MPI_Wtime ();
    /* DO NOT DELETE THE FOLLOWING line.
     * even if you do not want this output. It fixes a bug that occured on JUQUEEN, where the
     * runtimes were computed to 0.
     * Only delete the line, if you know what you are doing. */
    t8_global_productionf ("End adapt %f %f\n", sc_MPI_Wtime (),
                           forest->profile->adapt_runtime);
  }
}

/* A refinement callback function that works with a marker array.
 * We expect the forest user_data pointer to point to an sc_array of short ints.
 * We expect one int per element, which is interpreted as follows:
 *    0 - this element should neither be refined nor coarsened
 *    1 - refine this element
 *   -1 - coarsen this elements (must be set for the whole family)
 * other values are invalid.
 * The array may be longer (to include ghosts for example), but the additional
 * entries are ignored.
 */
int
t8_forest_adapt_marker_array_callback (t8_forest_t forest,
                                       t8_forest_t forest_from,
                                       t8_locidx_t which_tree,
                                       t8_locidx_t lelement_id,
                                       t8_eclass_scheme_c * ts,
                                       int num_elements,
                                       const t8_element_t * elements[])
{
  T8_ASSERT (t8_forest_is_committed (forest_from));
  /* Get a pointer to the user data. */
  sc_array_t         *markers =
    (sc_array_t *) t8_forest_get_user_data (forest);
  /* The element has to have (at least) as many entries as elements.
   * We allow additional entries, since they may be used outside this 
   * function. For example for marker values for the ghost elements.
   */
  T8_ASSERT ((t8_locidx_t) markers->elem_count >=
             t8_forest_get_num_element (forest_from));
  T8_ASSERT (markers->elem_size == sizeof (short));

  /* Get the (process local) index of the current element by adding the tree offset
   * to the tree local index. */
  t8_locidx_t         element_index =
    t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;
  /* Get the marker value for this element */
  short               marker_value =
    *(short *) (t8_sc_array_index_locidx (markers, element_index));
  T8_ASSERT (marker_value == 0 || marker_value == 1 || marker_value == -1);

#ifdef T8_ENABLE_DEBUG
  /* In debugging mode we check that
   * if one member of a family is marked -1, then all are
   */
  if (num_elements > 1) {
    int                 ielem;
    int                 oneisminusone = 0;
    int                 oneisnotminusone = 0;
    short               temp_marker_value;
    /* Iterate over elements and store if any element is -1
     * and if any element is not -1. */
    for (ielem = 0; ielem < num_elements; ++ielem) {
      temp_marker_value =
        *(short *) (t8_sc_array_index_locidx (markers, element_index));
      if (temp_marker_value == -1) {
        oneisminusone = 1;
      }
      else {
        oneisnotminusone = 1;
      }
    }
    /* If any element is -1, then all elements must be -1, thus
     * no element is allowed to be not -1. */
    T8_ASSERT (!(oneisminusone && oneisnotminusone));
  }
#endif
  /* Return the marker of this element */
  return marker_value;
}

/* Given a forest that is to be adapted, we fill an array with refinement
 * markers. Thus, for each element we store either 0, 1, or -1, depending
 * on what will happen with the element during refinement.
 *  0 - nothing
 *  1 - refine this element
 * -1 - coarsen this element and all its siblings.
 * 
 * Thus, the value of the markers corresponds to the difference in refinement
 * levels.
 */
void
t8_forest_adapt_build_marker_array (t8_forest_t forest, sc_array_t * markers,
                                    t8_locidx_list_t *
                                    elements_that_do_not_refine)
{
  t8_forest_t         forest_from = forest->set_from;
  T8_ASSERT (t8_forest_is_committed (forest_from));
  T8_ASSERT (forest->set_adapt_fn != NULL);
  T8_ASSERT (forest->set_adapt_recursive == 0);
  const t8_locidx_t   num_trees = t8_forest_get_num_local_trees (forest_from);
  t8_locidx_t         ltreeid, ielement;
  const t8_locidx_t   num_elements = t8_forest_get_num_element (forest_from);

  /* check correct markers array size */
  T8_ASSERT (markers != NULL);
  T8_ASSERT ((t8_locidx_t) markers->elem_count >= num_elements);
  T8_ASSERT (markers->elem_size == sizeof (short));
  /* Check correct list size */
  T8_ASSERT (t8_locidx_list_is_initialized (elements_that_do_not_refine));
  T8_ASSERT (t8_locidx_list_count (elements_that_do_not_refine) == 0);

  for (ltreeid = 0; ltreeid < num_trees; ltreeid++) {
    /* Get the tree's class, number of elements and scheme */
    const t8_eclass_t   tree_class =
      t8_forest_get_tree_class (forest_from, ltreeid);
    const t8_locidx_t   elements_in_tree =
      t8_forest_get_tree_num_elements (forest_from, ltreeid);
    const t8_eclass_scheme_c *ts =
      t8_forest_get_eclass_scheme (forest_from, tree_class);
    /* Iterate over all elements of this tree and call the refinement function */
    for (ielement = 0; ielement < elements_in_tree; ++ielement) {
      /* Get a pointer to the element */
      const t8_element_t *element =
        t8_forest_get_element_in_tree (forest_from, ltreeid, ielement);
      t8_element_t      **siblings;
      /* Compute the number of siblings */
      const int           num_siblings =
        ts->t8_element_num_siblings (element);
      int                 isib;
      int                 is_family = 0;
      /* Compute whether this element and its next num_siblings followers are a family */
      if (ielement + num_siblings <
          elements_in_tree && ts->t8_element_child_id (element) == 0) {
        siblings = T8_ALLOC (t8_element_t *, num_siblings);
        /* Gather the siblings in an array */
        for (isib = 0; isib < num_siblings; ++isib) {
          siblings[isib] =
            t8_forest_get_element_in_tree (forest_from, ltreeid,
                                           ielement + isib);
        }
        is_family =
          ts->t8_element_is_family ((const t8_element_t **) siblings);
        T8_FREE (siblings);
      }
      /* If this is a family we pass the element and all siblings to the adapt_fn,
       * otherwise only the element. */
      const int           num_elements_to_adapt =
        is_family ? num_siblings : 1;
      const int           adapt_value =
        forest->set_adapt_fn (forest, forest->set_from, ltreeid,
                              ielement, (t8_eclass_scheme_c *) ts,
                              num_elements_to_adapt, &element);
      /* TODO: Use const in parameters of adapt_fn and get rid of the type casts here. */
      const t8_locidx_t   element_index =
        ielement + t8_forest_get_tree_element_offset (forest_from, ltreeid);
      /* Set the marker to 1, 0, or -1 depending on whether adapt_value is >0, 0, <0 */
      *(short *) t8_sc_array_index_locidx (markers, element_index) =
        adapt_value > 0 ? 1 : adapt_value == 0 ? 0 : -1;
      if (adapt_value <= 0) {
        /* This element does not get refined, we add it to the list of unrefined element */
        t8_locidx_list_append (elements_that_do_not_refine, element_index);
      }
      if (adapt_value < 0) {
        /* This is a family that will get coarsened. At the marker value to all siblings.
         * and add all siblings to the list of not refined elements. */

        T8_ASSERT (is_family);
        for (isib = 1; isib < num_siblings; ++isib) {
          *(short *) t8_sc_array_index_locidx (markers,
                                               element_index + isib) = -1;
          t8_locidx_list_append (elements_that_do_not_refine,
                                 element_index + isib);
        }
        /* We have to skip the siblings from the adapt check, since we know that they
         * will be coarsened. However, calling the adapt function for the seperately
         * will return 0 or 1 and result in false markers. */
        ielement += num_siblings - 1;
      }
    }                           /* End of element loop */
  }                             /* End of tree loop */
  /* We expect that no element was put twice into the list. */
  T8_ASSERT (!t8_locidx_list_has_duplicate_entries
             (elements_that_do_not_refine));
}

T8_EXTERN_C_END ();
