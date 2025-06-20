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

#ifndef _COUPLER_HPP_
#define _COUPLER_HPP_

#include <mpi.h>
#include <hdf5.h>
#include <vector>
#include <string>
#include <cstdio>
#include <tuple>
#include <map>
#include <cstdint>

/* Map from field name to the HDF5 data type, the number of components per
 * object and a pointer to the data. */
using FieldMap_out = std::map< std::string, std::tuple< hid_t, std::int64_t, void * > >;
using FieldMap_in = std::map< std::string, std::tuple< hid_t, std::int64_t, const void * > >;

void waitForFileExistence( MPI_Comm comm, const char * filename );

/*!
 * \brief Write out a boundary file with the given data.
 *
 * \param [in] comm the communicator used in writing the file.
 * \param [in] filename the name of the file to write out to.
 * \param [in] dt the current time step.
 * \param [in] on_boundary true iff the respective face is to be written out.
 * \param [out] face_offset the global offset at which this ranks face data is
 *  written.
 * \param [out] n_faces_to_write the number of faces this rank will write out.
 * \param [in] n_faces the number of faces in the entire local mesh.
 * \param [in] faces connectivity of the faces which are assumed to be quads.
 * \param [in] face_fields map from face field names to fields.
 * \param [out] node_offset the global offset at which this ranks node data is
 *  written out.
 * \param [out] n_nodes_to_write the number of nodes this rank will write out.
 * \param [in] n_nodes the number of nodes in the entire local mesh.
 * \param [in] node_fields map from node field names to fields.
 */
void writeBoundaryFile( MPI_Comm comm, const char * filename, double dt, const bool * on_boundary,
                        std::int64_t & face_offset, std::int64_t & n_faces_to_write, std::int64_t n_faces,
                        const std::int64_t * faces, const FieldMap_in & face_fields,
                        std::int64_t & node_offset, std::int64_t & n_nodes_to_write, std::int64_t n_nodes,
                        const FieldMap_in & node_fields );

/*!
 * \brief Read in a boundary file into the provided fields.
 *
 * \param [in] comm the communicator used in reading from the file.
 * \param [in] filename the name of the file to read from.
 * \param [in] face_offset the global offset at which this rank will read
 *  face data.
 * \param [in] n_faces_to_read the number of faces this rank will read in.
 * \param [in] n_faces the number of faces in the entire local mesh.
 * \param [in] node_offset the global offset at which this rank will read node
 *  data.
 * \param [in] n_nodes_to_read the number of nodes this rank will read in.
 * \param [in] n_nodes the number of nodes in the entire local mesh.
 * \param [in/out] face_fields map from face field names to fields.
 * \param [in/out] node_fields map from node field names to fields.
 */
void readBoundaryHeader( MPI_Comm comm,
                         const char * filename,
                         double & dt,
                         std::int64_t & n_faces,
                         std::int64_t & n_nodes );

void readBoundaryFile( MPI_Comm comm, const char * filename,
                       std::int64_t face_offset, std::int64_t n_faces_to_read, std::int64_t n_faces, FieldMap_out & face_fields,
                       std::int64_t node_offset, std::int64_t n_nodes_to_read, std::int64_t n_nodes, FieldMap_out & node_fields );

#endif
