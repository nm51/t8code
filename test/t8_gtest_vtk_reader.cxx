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

/** \file t8_gtest_vtk_reader.cxx
* Provide tests to check the functionality of the vtk-reader. Tests each supported
* file-format for each reader. 
*/

#include <gtest/gtest.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh.h>
#include <t8_cmesh_vtk_writer.h>

#if T8_WITH_VTK
#include <t8_cmesh/t8_cmesh_vtk_reader.hxx>
#endif

/* Number of supported file-types that store poly-data */
/* TODO: Find .g example*/
#ifndef T8_POLY_SUPPORT
#define T8_POLY_SUPPORT 4
#endif

/* Number of total supported file-types (poly-data + unstructured)*/
#ifndef T8_TOTAL_SUPPORT
#define T8_TOTAL_SUPPORT 6
#endif

/* The number of the file-type that refers to .vtu files*/
#ifndef T8_VTU_TYPE
#define T8_VTU_TYPE 5
#endif

const char          t8_poly_support[T8_POLY_SUPPORT][4] =
  { "ply", "vtp", "obj", "stl" };
const char          t8_unstructured_support[2][4] = { "vtu", "vtk" };

/** A function to read from a file. The filetype is specified via file_type
 * Needs a communicator for future parallel file I/O
 * \param [in] filename     The path to the file that should be read
 * \param [in] file_type    The type of the file given by \a filename. The order of the
 *                          filetypes is defined via the arrays t8_poly_support and t8_unstructured_support
 * \param [in] comm         The mpi-communicator to use
 * \returns                 A commited cmesh 
*/
t8_cmesh_t
t8_read_from_vtk (const char *filename, int file_type, sc_MPI_Comm comm)
{
/* Check if file exists. */
  SC_CHECK_ABORTF (access (filename, R_OK) == 0, "Could not open file %s\n",
                   filename);
#if T8_WITH_VTK
  if (file_type < T8_POLY_SUPPORT) {
    return t8_cmesh_read_from_vtk_poly (filename, 1, 0, comm);
  }
  else {
    return t8_cmesh_read_from_vtk_unstructured (filename, 1, 0, comm);
  }
#else
  return NULL;
#endif /* T8_WITH_VTK */
}

/* *INDENT-OFF* */
class vtk_reader:public testing::TestWithParam < int > {
protected:
    void SetUp () override {
#if T8_WITH_VTK
        /* we iterate over the test files. */
        file_type = GetParam ();
        if(file_type < T8_POLY_SUPPORT)
        {
            snprintf (filename, BUFSIZ, "test/testfiles/vtk_reader/simple_test.%s", 
                t8_poly_support[file_type]);
        }
        else
        {
            snprintf (filename, BUFSIZ, "test/testfiles/vtk_reader/simple_test.%s", 
                t8_unstructured_support[file_type - T8_POLY_SUPPORT]);
        }
#else
        /* Do nothing if we do not link with VTK */
        GTEST_SKIP();
#endif /* T8_WITH_VTK */
    }   
    void TearDown () override {
#if T8_WITH_VTK
        /*Nothing to tear down.*/
#else
        GTEST_SKIP();
#endif /* T8_WITH_VTK */ 
    }
    char        filename[BUFSIZ];
    int         file_type;
}; 

/* Read the cmesh and check if the cmesh is committed */
TEST_P(vtk_reader, simple_test)
{
    t8_cmesh_t          cmesh = t8_read_from_vtk(filename, file_type, sc_MPI_COMM_WORLD);
    t8_debugf("[D] %s\n", filename);
    EXPECT_TRUE(t8_cmesh_is_committed(cmesh));
    t8_cmesh_destroy(&cmesh);
}

/* Construct the hypercube_hybrid mesh and write it. Read the file and compare
 * the two cmeshes. Currently we can check this only for vtu files, as the writer
 * only supports vtu files.
 * TODO: This check fails currently for vtu files. Fix this.*/
TEST_P(vtk_reader, compare_constructed)
{
    /* Run the test only for vtu*/
    if(file_type != T8_VTU_TYPE){
        GTEST_SKIP();
    }
    else{
        /* Construct the hybrid mesh*/
        t8_cmesh_t          cmesh = t8_cmesh_new_hypercube_hybrid(sc_MPI_COMM_WORLD, 0, 0, 0);
        t8_cmesh_t          cmesh_from_reader;
        t8_locidx_t         num_trees;
        t8_locidx_t         itree;
        t8_eclass_t         eclass;
        /* Write the hybrid mesh*/
        t8_cmesh_vtk_write_file(cmesh, "test/testfiles/vtk_reader/hybrid", 1.0);
        /* Currently we can only read vtu files. 
         * TODO: Change this, as soon as we can read pvtu files */
        cmesh_from_reader = t8_read_from_vtk("test/testfiles/vtk_reader/hybrid_0000.vtu", file_type, sc_MPI_COMM_WORLD); 
        /* Check if the number of trees is the same */
        num_trees = t8_cmesh_get_num_local_trees(cmesh);
        EXPECT_EQ(num_trees, t8_cmesh_get_num_local_trees(cmesh_from_reader));
        /* Check if the class of every tree is the same. */
        for(itree = 0; itree < num_trees; itree++){
            eclass = t8_cmesh_get_tree_class(cmesh, itree);
            EXPECT_EQ(eclass, t8_cmesh_get_tree_class(cmesh_from_reader, itree));
        }

        /* This test fails currently*/
         EXPECT_TRUE(t8_cmesh_is_equal(cmesh, cmesh_from_reader));

        t8_cmesh_destroy(&cmesh);
        t8_cmesh_destroy(&cmesh_from_reader);
    }
    
}

INSTANTIATE_TEST_SUITE_P (t8_gtest_poly, vtk_reader,
                        testing::Range(0, T8_TOTAL_SUPPORT));
/* *INDENT-ON* */
