/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "coupler.H"
#include <iostream>
#include <cmath>

constexpr double STRIDE = 1.0;
constexpr int X_EXTENT = 20;
constexpr int Y_EXTENT = 20;
constexpr int VCOUNT = X_EXTENT * Y_EXTENT;
constexpr int QCOUNT = (X_EXTENT - 1) * (Y_EXTENT - 1);

int GLOBAL_X_OFFSET;
int GLOBAL_Y_OFFSET;

int stop( int, void * )
{
  MPI_Abort( MPI_COMM_WORLD, 1 );
  return 0;
}

void setNodePos( int i, int j, vec3 & pos )
{
  const int x_offset = j + X_EXTENT * GLOBAL_X_OFFSET;
  const int y_offset = i + Y_EXTENT * GLOBAL_Y_OFFSET;

  pos.e0 = STRIDE * x_offset;
  pos.e1 = STRIDE * y_offset;
  pos.e2 = 0.0;
}

bool checkNodePos( int i, int j, const vec3 & pos )
{
  vec3 temp;
  setNodePos( i, j, temp );
  return pos.e0 == temp.e0 && pos.e1 == temp.e1 && pos.e2 == temp.e2;
}

void setNodeVel( int i, int j, vec3 & vel )
{
  const int x_offset = j + X_EXTENT * GLOBAL_X_OFFSET;
  const int y_offset = i + Y_EXTENT * GLOBAL_Y_OFFSET;

  vel.e0 = 1.0 / x_offset;
  vel.e1 = 1.0 / y_offset;
  vel.e2 = 0.0;
}

bool checkNodeVel( int i, int j, const vec3 & vel )
{
  vec3 temp;
  setNodeVel( i, j, temp );
  return vel.e0 == temp.e0 && vel.e1 == temp.e1 && vel.e2 == temp.e2;
}

void invalidateVec3( vec3 & v )
{
  v.e0 = -1.0;
  v.e1 = -1.0;
  v.e2 = -1.0;
}

bool onBoundary( int i, int j )
{
  return !(i == 0 || j == 0 || i == Y_EXTENT - 2 || j == X_EXTENT - 2);
}

int main( int argc, char * argv[] )
{
  MPI_Init( &argc, &argv );
  int rank, N;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &N );

  // for debugging, tell HDF5 to abort on error
  H5Eset_auto ( H5E_DEFAULT, stop, NULL );

  const int GLOBAL_DIM = std::sqrt( N );
  if( GLOBAL_DIM * GLOBAL_DIM != N )
  {
    std::cout << "Number of ranks must be a perfect square." << std::endl;
    MPI_Abort( MPI_COMM_WORLD, 1 );
  }

  GLOBAL_X_OFFSET = rank % GLOBAL_DIM;
  GLOBAL_Y_OFFSET = rank / GLOBAL_DIM;

  const int GLOBAL_NODE_OFFSET = rank * VCOUNT;
  const int GLOBAL_FACE_OFFSET = rank * QCOUNT;

  // GEOS would replace this logic here.
  vec3 * node_pos = new vec3[VCOUNT];
  vec3 * node_vel = new vec3[VCOUNT];
  int * quads = new int[4 * QCOUNT];
  bool * on_boundary = new bool[QCOUNT];
  double * pressure = new double[QCOUNT];

  /* Initialize the nodes. */
  int node_index = 0;
  for( int i = 0; i < Y_EXTENT; ++i )
  {
    for( int j = 0; j < X_EXTENT; ++j )
    {
      setNodePos( i, j, node_pos[node_index] );
      setNodeVel( i, j, node_vel[node_index] );
      node_index++;
    }
  }

  /* Initialize the quads. */
  int quad_index = 0;
  for( int i = 0; i < Y_EXTENT - 1; ++i )
  {
    for( int j = 0; j < X_EXTENT - 1; ++j )
    {
      quads[4 * quad_index] = X_EXTENT * i + j;
      quads[4 * quad_index + 1] = X_EXTENT * i + j + 1;
      quads[4 * quad_index + 2] = X_EXTENT * (i + 1) + j + 1;
      quads[4 * quad_index + 3] = X_EXTENT * (i + 1) + j;

      pressure[quad_index] = 1.0 / (quad_index + GLOBAL_FACE_OFFSET * rank);

      if( onBoundary( i, j ))
      {
        on_boundary[quad_index] = true;
      }
      else
      {
        on_boundary[quad_index] = false;
      }

      quad_index++;
    }
  }

  /* Write out the data on the boundary. */
  double dt = 2.5;
  FieldMap face_fields;
  face_fields["pressure"] = std::make_tuple( H5T_NATIVE_DOUBLE, 1, pressure );

  FieldMap node_fields;
  node_fields["position"] = std::make_tuple( H5T_NATIVE_DOUBLE, 3, node_pos );
  node_fields["velocity"] = std::make_tuple( H5T_NATIVE_DOUBLE, 3, node_vel );

  int face_offset, n_faces_written, node_offet, n_nodes_written;
  writeBoundaryFile( MPI_COMM_WORLD, "GEOSboundary.hdf5", face_offset,
                     n_faces_written, node_offet, n_nodes_written, dt, QCOUNT,
                     VCOUNT, quads, on_boundary, face_fields, node_fields );

  /* Invalidate the data on the boundary. */
  quad_index = 0;
  for( int i = 0; i < Y_EXTENT - 1; ++i )
  {
    for( int j = 0; j < X_EXTENT - 1; ++j )
    {
      if( onBoundary( i, j ))
      {
        pressure[quad_index] = -1.0;
        for( int k = 0; k < 4; ++k )
        {
          int node = quads[4 * quad_index + k];
          invalidateVec3( node_pos[node] );
          invalidateVec3( node_vel[node] );
        }
      }

      quad_index++;
    }
  }

  /* Read in the boundary file. */
  readBoundaryFile( MPI_COMM_WORLD, "GEOSboundary.hdf5", face_offset,
                    n_faces_written, QCOUNT, node_offet, n_nodes_written, VCOUNT,
                    face_fields, node_fields );

  /* Check that the data is back to normal. */
  bool success = true;
  node_index = 0;
  for( int i = 0; i < Y_EXTENT; ++i )
  {
    for( int j = 0; j < X_EXTENT; ++j )
    {
      bool result = true;
      result &= checkNodePos( i, j, node_pos[node_index] );
      result &= checkNodeVel( i, j, node_vel[node_index] );
      success &= result;

      if( !result )
      {
        std::cout << "Rank " << rank << ": Node error at (" << i << ", " << j << ")." << std::endl;
      }

      node_index++;
    }
  }

  quad_index = 0;
  for( int i = 0; i < Y_EXTENT - 1; ++i )
  {
    for( int j = 0; j < X_EXTENT - 1; ++j )
    {
      if( pressure[quad_index] != 1.0 / (quad_index + GLOBAL_FACE_OFFSET * rank))
      {
        std::cout << "Rank " << rank << ": Face error at (" << i << ", " << j
                  << ")." << std::endl;
        success = false;
      }

      quad_index++;
    }
  }

  if( success )
  {
    std::cout << "Rank " << rank <<  ": TESTS PASS!!!" << std::endl;
  }
  else
  {
    std::cout << "Rank " << rank <<  ": TESTS FAIL!!!" << std::endl;
  }

  delete[] node_pos;
  delete[] node_vel;
  delete[] quads;
  delete[] on_boundary;
  delete[] pressure;

  MPI_Finalize();

  return 0;
}
