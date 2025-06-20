/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file EmebeddedSurfacesParallelSynchronization.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACESPARALLELSYNCHRONIZATION_HPP_
#define GEOS_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACESPARALLELSYNCHRONIZATION_HPP_

#include "common/DataTypes.hpp"
#include "common/StdContainerWrappers.hpp"

namespace geos
{
class MeshLevel;
class NeighborCommunicator;

struct NewObjectLists
{
  std::set< localIndex > newNodes;
  std::set< localIndex > newEdges;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newElements;

  void insert( NewObjectLists const & newObjects )
  {
    newNodes.insert( newObjects.newNodes.begin(),
                     newObjects.newNodes.end() );

    newEdges.insert( newObjects.newEdges.begin(),
                     newObjects.newEdges.end() );

    for( auto & iter : newObjects.newElements )
    {
      std::pair< localIndex, localIndex > const & key = iter.first;
      std::set< localIndex > const & values = iter.second;
      newElements[key].insert( values.begin(), values.end() );
    }

  }

};

namespace embeddedSurfacesParallelSynchronization
{

void synchronizeTopology( MeshLevel & mesh,
                          stdVector< NeighborCommunicator > & neighbors,
                          NewObjectLists & newObjects,
                          int const mpiCommOrder,
                          string const fractureRegionName );
}

} /* namespace geos */

#endif /* GEOS_PHYSICSSOLVERS_SURFACEGENERATION_EMBEDDEDSURFACESPARALLELSYNCHRONIZATION_HPP_ */
