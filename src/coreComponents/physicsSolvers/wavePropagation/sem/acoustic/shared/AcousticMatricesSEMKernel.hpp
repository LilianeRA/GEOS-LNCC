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
 * @file AcousticMatricesSEMKernel.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_
#define GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_

namespace geos
{

struct AcousticMatricesSEM
{

  template< typename FE_TYPE >
  struct MassMatrix
  {

    MassMatrix( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}
    /**
     * @brief Launches the precomputation of the mass matrices
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @tparam FE_TYPE the type of discretization
     * @param[in] finiteElement The finite element discretization used
     * @param[in] size the number of cells in the subRegion
     * @param[in] numFacesPerElem number of faces per element
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[out] mass diagonal of the mass matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeMassMatrix( localIndex const size,
                       arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                       arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                       arrayView1d< real32 const > const velocity,
                       arrayView1d< real32 const > const density,
                       arrayView1d< real32 > const mass )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {

        real32 const invC2 = 1.0 / ( density[e] * pow( velocity[e], 2 ) );
        // only the eight corners of the mesh cell are needed to compute the Jacobian
        real64 xLocal[ 8 ][ 3 ];
        for( localIndex a = 0; a < 8; ++a )
        {
          localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
          for( localIndex i = 0; i < 3; ++i )
          {
            xLocal[a][i] = nodeCoords( nodeIndex, i );
          }
        }
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
        for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
        {
          real32 const localIncrement = invC2 * m_finiteElement.computeMassTerm( q, xLocal );
          RAJA::atomicAdd< ATOMIC_POLICY >( &mass[elemsToNodes( e, q )], localIncrement );
        }
      } ); // end loop over element
    }

    FE_TYPE const & m_finiteElement;

  };
  template< typename FE_TYPE >
  struct DampingMatrix
  {

    DampingMatrix( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launches the precomputation of the damping matrices
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToFaces map from elements to faces
     * @param[in] facesToNodes map from face to nodes
     * @param[in] facesDomainBoundaryIndicator flag equal to 1 if the face is on the boundary, and to 0 otherwise
     * @param[in] freeSurfaceFaceIndicator flag equal to 1 if the face is on the free surface, and to 0 otherwise
     * @param[in] velocity cell-wise velocity
     * @param[in] density cell-wise density
     * @param[out] damping diagonal of the damping matrix
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeDampingMatrix( localIndex const size,
                          arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                          arrayView2d< localIndex const > const elemsToFaces,
                          ArrayOfArraysView< localIndex const > const facesToNodes,
                          arrayView1d< integer const > const facesDomainBoundaryIndicator,
                          arrayView1d< localIndex const > const freeSurfaceFaceIndicator,
                          arrayView1d< real32 const > const velocity,
                          arrayView1d< real32 const > const density,
                          arrayView1d< real32 > const damping )
    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        for( localIndex i = 0; i < elemsToFaces.size( 1 ); ++i )
        {
          localIndex const f = elemsToFaces( e, i );
          // face on the domain boundary and not on free surface
          if( facesDomainBoundaryIndicator[f] == 1 && freeSurfaceFaceIndicator[f] != 1 )
          {
            // only the four corners of the mesh face are needed to compute the Jacobian
            real64 xLocal[ 4 ][ 3 ];
            for( localIndex a = 0; a < 4; ++a )
            {
              localIndex const nodeIndex = facesToNodes( f, FE_TYPE::meshIndexToLinearIndex2D( a ) );
              for( localIndex d = 0; d < 3; ++d )
              {
                xLocal[a][d] = nodeCoords( nodeIndex, d );
              }
            }
            real32 const alpha = 1.0 / (density[e] * velocity[e]);
            constexpr localIndex numNodesPerFace = FE_TYPE::numNodesPerFace;
            for( localIndex q = 0; q < numNodesPerFace; ++q )
            {
              real32 const localIncrement = alpha * m_finiteElement.computeDampingTerm( q, xLocal );
              RAJA::atomicAdd< ATOMIC_POLICY >( &damping[facesToNodes( f, q )], localIncrement );
            }
          }
        }
      } );
    }

    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;

  };

  template< typename FE_TYPE >
  struct GradientKappaBuoyancy
  {

    GradientKappaBuoyancy( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launch the computation of the 2 gradients relative to the coeff of the wave equation K=1/rho*c2 and b=1/rho
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] q_dt2 second order derivative in time of backward
     * @param[in] q_n current time step of backward
     * @param[in] p_n current time step of forward
     * @param[out] grad first part of gradient vector with respect to K=1/rho*c2
     * @param[out] grad2 second part of gradient vector with respact to b=1/rho
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeGradient( localIndex const size,
                     arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                     arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                     arrayView1d< integer const > const elemGhostRank,
                     arrayView1d< real32 const > const q_dt2,
                     arrayView1d< real32 const > const q_n,
                     arrayView1d< real32 const > const p_n,
                     arrayView1d< real32 > const grad,
                     arrayView1d< real32 > const grad2 )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        if( elemGhostRank[e]<0 )
        {
          // only the eight corners of the mesh cell are needed to compute the Jacobian
          real64 xLocal[ 8 ][ 3 ];
          for( localIndex a = 0; a < 8; ++a )
          {
            localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = nodeCoords( nodeIndex, i );
            }
          }
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            localIndex nodeIdx = elemsToNodes( e, q );
            grad[e] += q_dt2[nodeIdx] * p_n[nodeIdx] * m_finiteElement.computeMassTerm( q, xLocal );
            m_finiteElement.template computeStiffnessTerm<>( q, xLocal, [&] ( const int i, const int j, const real64 val )
            {
              grad2[e] += val* q_n[elemsToNodes( e, j )] * p_n[elemsToNodes( e, i )];
            } );
          }
        }
      } ); // end loop over element
    }
    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;
  };

  template< typename FE_TYPE >
  struct ImagingCondition
  {

    ImagingCondition( FE_TYPE const & finiteElement )
      : m_finiteElement( finiteElement )
    {}

    /**
     * @brief Launch the computation of the imaging condition for RTM
     * @tparam EXEC_POLICY the execution policy
     * @tparam ATOMIC_POLICY the atomic policy
     * @param[in] size the number of cells in the subRegion
     * @param[in] nodeCoords coordinates of the nodes
     * @param[in] elemsToNodes map from element to nodes
     * @param[in] q_n current time step of backward
     * @param[in] p_n current time step of forward
     * @param[out] imag will contain the imaging condition value
     */
    template< typename EXEC_POLICY, typename ATOMIC_POLICY >
    void
    computeImagingCondition( localIndex const size,
                             arrayView2d< WaveSolverBase::wsCoordType const, nodes::REFERENCE_POSITION_USD > const nodeCoords,
                             arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                             arrayView1d< integer const > const elemGhostRank,
                             arrayView1d< real32 const > const q_n,
                             arrayView1d< real32 const > const p_n,
                             arrayView1d< real32 > const imag )

    {
      forAll< EXEC_POLICY >( size, [=] GEOS_HOST_DEVICE ( localIndex const e )
      {
        if( elemGhostRank[e]<0 )
        {
          // only the eight corners of the mesh cell are needed to compute the Jacobian
          real64 xLocal[ 8 ][ 3 ];
          for( localIndex a = 0; a < 8; ++a )
          {
            localIndex const nodeIndex = elemsToNodes( e, FE_TYPE::meshIndexToLinearIndex3D( a ) );
            for( localIndex i = 0; i < 3; ++i )
            {
              xLocal[a][i] = nodeCoords( nodeIndex, i );
            }
          }
          constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            localIndex nodeIdx = elemsToNodes( e, q );
            imag[e] += q_n[nodeIdx] * p_n[nodeIdx] * m_finiteElement.computeMassTerm( q, xLocal );
          }
        }
      } ); // end loop over element
    }
    /// The finite element space/discretization object for the element type in the subRegion
    FE_TYPE const & m_finiteElement;
  };

};

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICMATRICESSEMKERNEL_HPP_
