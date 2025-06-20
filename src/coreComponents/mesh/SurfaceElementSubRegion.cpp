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
 * @file SurfaceElementSubRegion.cpp
 */


#include "SurfaceElementSubRegion.hpp"
#include "ElementRegionManager.hpp"
#include "MeshFields.hpp"

namespace geos
{

using namespace dataRepository;

SurfaceElementSubRegion::SurfaceElementSubRegion( string const & name,
                                                  dataRepository::Group * const parent ):
  ElementSubRegionBase( name, parent ),
  m_unmappedGlobalIndicesInToNodes(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_elementAperture(),
  m_elementArea(),
  m_normalVector(),
  m_tangentVector1(),
  m_tangentVector2()
{
  registerWrapper( viewKeyStruct::nodeListString(), &m_toNodesRelation ).
    setDescription( "Map to the nodes attached to each SurfaceElement." );

  registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation ).
    setDescription( "Map to the edges attached to each SurfaceElement." );

  registerField( fields::elementAperture{}, &m_elementAperture );

  registerField( fields::elementArea{}, &m_elementArea );

  registerField( fields::normalVector{}, &m_normalVector ).
    reference().resizeDimension< 1 >( 3 );

  registerField( fields::tangentVector1{}, &m_tangentVector1 ).
    reference().resizeDimension< 1 >( 3 );

  registerField( fields::tangentVector2{}, &m_tangentVector2 ).
    reference().resizeDimension< 1 >( 3 );

  excludeWrappersFromPacking( { viewKeyStruct::nodeListString(),
                                viewKeyStruct::edgeListString(),
                                viewKeyStruct::surfaceElementsToCellRegionsString(),
                                viewKeyStruct::surfaceElementsToCellSubRegionsString(),
                                viewKeyStruct::surfaceElementsToCellIndexString() } );


}

SurfaceElementSubRegion::~SurfaceElementSubRegion()
{}

} /* namespace geos */
