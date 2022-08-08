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

/** \file t8_cmesh_vtk_reader.cxx
* Implementation of a Reader for vtk/vtu files using the vtk-library.
* The functions can only be used when t8code is linked with the vtk-library.
*/

#include <t8_cmesh.h>
#include <t8_data/t8_shmem.h>
#include <t8_cmesh_vtk_writer.h>
#include <t8_cmesh_vtk_reader.hxx>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_vtk_helper.hxx>

#if T8_WITH_VTK
#include <vtkUnstructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkTriangleFilter.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformation.h>
#endif
T8_EXTERN_C_BEGIN ();

t8_cmesh_t
t8_cmesh_parallel_read_from_vtk_unstructured (const char *filename,
                                              sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  t8_cmesh_t          cmesh;
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkSmartPointer < vtkCellData > cellData;
  vtkNew < vtkXMLPUnstructuredGridReader > p_reader;
  int                 mpiret;
  int                 mpisize;
  int                 mpirank;
  int                 main_proc = 0;    /*TODO: Should be an argument */
  t8_gloidx_t         local_num_trees;
  t8_gloidx_t         global_num_trees;
  t8_gloidx_t         first_tree = 0;
  t8_gloidx_t         last_tree;
  int                 rank_it;
  t8_shmem_init (comm);
  p_reader->SetFileName (filename);
  p_reader->UpdateInformation ();
  mpiret = sc_MPI_Comm_size (comm, &mpisize);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (comm, &mpirank);
  SC_CHECK_MPI (mpiret);
  vtkInformation     *Info = p_reader->GetOutputInformation (0);
  Info->Set (vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER (),
             mpirank);
  Info->Set (vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES (),
             mpisize);
  p_reader->Update ();

  unstructuredGrid = p_reader->GetOutput ();
  cellData = unstructuredGrid->GetCellData ();
  t8_cmesh_init (&cmesh);
  local_num_trees =
    (t8_gloidx_t) t8_vtk_iterate_cells (unstructuredGrid, cellData, comm,
                                        &cmesh);
  t8_shmem_array_t    tree_offsets = t8_cmesh_alloc_offsets (mpisize, comm);

  //t8_shmem_array_allgather(&local_num_trees, 1, T8_MPI_GLOIDX, tree_offsets, 1, T8_MPI_GLOIDX);

  t8_shmem_array_prefix (&local_num_trees, tree_offsets, 1, T8_MPI_GLOIDX,
                         sc_MPI_SUM, comm);

  for (rank_it = 0; rank_it <= mpisize; rank_it++) {
    t8_debugf ("[D] %i, offset: %li\n", rank_it,
               t8_shmem_array_get_gloidx (tree_offsets, rank_it));
  }
  first_tree = t8_shmem_array_get_gloidx (tree_offsets, mpirank);
  last_tree = first_tree + local_num_trees;

  t8_debugf ("[D] first_tree: %li, last_tree: %li\n", first_tree, last_tree);

  t8_cmesh_set_partition_range (cmesh, 3, first_tree, last_tree);

  t8_cmesh_commit (cmesh, comm);
  return cmesh;

#else
  /*Return empty cmesh if not linked against vtk */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return t8_cmesh_new_empty (comm, 0, 0);
#endif
}

/*Construct a cmesh given a filename and a*/
t8_cmesh_t
t8_cmesh_read_from_vtk_unstructured (const char *filename,
                                     const int num_files,
                                     const int compute_face_neigh,
                                     sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  t8_cmesh_t          cmesh;
  t8_cmesh_init (&cmesh);
  /*The Incoming data must be an unstructured Grid */
  vtkSmartPointer < vtkUnstructuredGrid > unstructuredGrid;
  vtkSmartPointer < vtkCellData > cellData;
  /*Prepare grid for translation */
  unstructuredGrid = t8_read_unstructured (filename);

  /* Get the Data of the all cells */
  cellData = unstructuredGrid->GetCellData ();

  /*Actual translation */
  t8_vtk_iterate_cells (unstructuredGrid, cellData, comm, &cmesh);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;

#else
  /*Return empty cmesh if not linked against vtk */
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return t8_cmesh_new_empty (comm, 0, 0);
#endif
}

t8_cmesh_t
t8_cmesh_read_from_vtk_poly (const char *filename, const int num_files,
                             const int compute_face_neigh, sc_MPI_Comm comm)
{
#if T8_WITH_VTK
  vtkSmartPointer < vtkPolyData > poly_data;
  vtkSmartPointer < vtkCellArray > cells;
  vtkSmartPointer < vtkCellData > cell_data;
  vtkSmartPointer < vtkPolyData > triangulated;
  vtkNew < vtkTriangleFilter > tri_filter;

  t8_cmesh_t          cmesh;
  t8_cmesh_init (&cmesh);

  /* Prepare the poly-data for the translation from vtk to t8code.
   * We split all polygons (which are not supported by t8code) to
   * triangles, vertices and lines. */
  poly_data = t8_read_poly (filename);
  tri_filter->SetInputData (poly_data);
  /*PolyVertex to vertex */
  tri_filter->PassVertsOn ();
  /*PolyLines to lines */
  tri_filter->PassLinesOn ();
  tri_filter->Update ();
  triangulated = tri_filter->GetOutput ();

  cell_data = triangulated->GetCellData ();

  t8_vtk_iterate_cells (triangulated, cell_data, comm, &cmesh);
  t8_cmesh_commit (cmesh, comm);
  return cmesh;
#else
  t8_global_errorf
    ("WARNING: t8code is not linked against the vtk library. Without proper linking t8code cannot use the vtk-reader\n");
  return t8_cmesh_new_empty (comm, 0, 0);
#endif
}

T8_EXTERN_C_END ();
