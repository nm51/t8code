/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element types in parallel.

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

#include <sc_refcount.h>
#include <t8_default_cxx.hxx>
#include <t8_default/t8_dprism.h>
#include <t8_default/t8_dtri.h>
#include <t8_default/t8_dtet.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest.h>
#include <sc_flops.h>
#include <sc_statistics.h>
#include <sc_options.h>



static int
t8_basic_adapt_refine_type (t8_forest_t forest, t8_locidx_t which_tree,
                            t8_eclass_scheme_c * ts,
                            int num_elements, t8_element_t * elements[])
{
  int                 level;
  int                 type;

  T8_ASSERT (num_elements == 1 || num_elements ==
             ts->t8_element_num_children (elements[0]));

  level = ts->t8_element_level(elements[0]);
  if (level >= *(int*)t8_forest_get_user_data(forest)) {
    return 0;
  }
  /* get the type of the current element */
  type = ((t8_dprism_t *) elements[0])->tri.type;
  /* refine type 0 */
  if (type == 0) {
    return 1;
  }
  return 0;
}

static void
t8_time_refine(int start_level, int end_level, int create_forest, int cube){
  t8_forest_t         forest, forest_adapt;
  sc_flopinfo_t       fi, snapshot;
  sc_statinfo_t       stats[1];
  char                vtuname[BUFSIZ];

  t8_forest_init (&forest);

  if(cube == 0){
  t8_forest_set_cmesh (forest,
                       t8_cmesh_new_bigmesh (T8_ECLASS_PRISM, 512, sc_MPI_COMM_WORLD),
                       sc_MPI_COMM_WORLD);
  }
  else{
      t8_forest_set_cmesh (forest,
                           t8_cmesh_new_hypercube (T8_ECLASS_PRISM, sc_MPI_COMM_WORLD,0,0),
                           sc_MPI_COMM_WORLD);
  }
  t8_forest_set_scheme (forest, t8_scheme_new_default_cxx ());
  t8_forest_set_level (forest, start_level);
  t8_forest_set_profiling (forest, 1);
  t8_forest_commit (forest);

  t8_forest_print_profile (forest);
  t8_cmesh_print_profile (t8_forest_get_cmesh (forest));
  snprintf (vtuname, BUFSIZ, "forest_hypercube_%s",
            t8_eclass_to_string[T8_ECLASS_PRISM]);
  t8_forest_write_vtk (forest, vtuname);
  t8_debugf ("Output to %s\n", vtuname);

  t8_forest_init (&forest_adapt);
  t8_forest_set_user_data(forest_adapt, &end_level);

  t8_forest_set_adapt (forest_adapt, forest,
                       t8_basic_adapt_refine_type, NULL, 1);

  sc_flops_start (&fi);
  sc_flops_snap (&fi, &snapshot);

  t8_forest_commit (forest_adapt);

  snprintf (vtuname, BUFSIZ, "forest_hypercube_adapt_%s",
            t8_eclass_to_string[T8_ECLASS_PRISM]);
  t8_forest_write_vtk (forest_adapt, vtuname);
  t8_debugf ("Output to %s\n", vtuname);

  sc_flops_shot (&fi, &snapshot);
  sc_stats_set1 (&stats[0], snapshot.iwtime, "New");


  t8_forest_unref (&forest_adapt);

  sc_stats_compute (sc_MPI_COMM_WORLD, 1, stats);
  sc_stats_print (t8_get_package_id (), SC_LP_STATISTICS, 1, stats, 1, 1);

}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  char                usage[BUFSIZ];
  char                help[BUFSIZ];
  int                 create_forest;
  int                 start_level = 0,end_level = 1, cube = 0;
  int                 parsed, helpme;

  /* brief help message */
  snprintf (usage, BUFSIZ, "Usage:\t%s <OPTIONS>\n\t%s -h\t"
            "for a brief overview of all options.",
            basename (argv[0]), basename (argv[0]));

  /* long help message */
  snprintf (help, BUFSIZ, "This program constructs a prism mesh of 512 prisms. "
            "\nThe user can choose the initial refinement level and the final\n"
            "refinement level of the mesh. If not set, the initial level is 0,\n"
            "the final level is 1. The program has no visual output, if desired,\n"
            "the user can switch to a hpyercube mesh.\n\n%s\n", usage);

  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  t8_init (SC_LP_DEFAULT);

  /* initialize command line argument parser */
  opt = sc_options_new (argv[0]);
  sc_options_add_switch (opt, 'h', "help", &helpme,
                         "Display a short help message.");
  sc_options_add_int (opt, 's', "slevel", &start_level, 0,
                      "initial refine level");
  sc_options_add_int (opt, 'e', "elevel", &end_level, 0,
                      "Final refine level: greater or equal to initial refine level");
  sc_options_add_int (opt, 'c', "cube", &cube, 0,
                      "cube = 1 -> use the hypercube mesh and visual output.");


  parsed =
    sc_options_parse (t8_get_package_id (), SC_LP_ERROR, opt, argc, argv);
  if(end_level < start_level){
      t8_debugf("Wrong usage of end and start level, end level set to start level\n");
      end_level = start_level;
  }
  if (helpme) {
    /* display help message and usage */
    t8_global_productionf ("%s\n", help);
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }
  else if (parsed >= 0 && 0 <= start_level && start_level <= end_level) {
    create_forest = 1;
    t8_time_refine(start_level, end_level, create_forest, cube);
  }
  else {
    /* wrong usage */
    t8_global_productionf ("\n\t ERROR: Wrong usage.\n\n");
    sc_options_print_usage (t8_get_package_id (), SC_LP_ERROR, opt, NULL);
  }

  sc_options_destroy (opt);
  sc_finalize ();

  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return 0;
}
