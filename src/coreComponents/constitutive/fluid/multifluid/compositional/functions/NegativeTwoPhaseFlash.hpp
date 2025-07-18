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
 * @file NegativeTwoPhaseFlash.hpp
 */

#ifndef GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
#define GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_

#include "RachfordRice.hpp"
#include "KValueInitialization.hpp"
#include "FugacityCalculator.hpp"
#include "FlashData.hpp"
#include "constitutive/fluid/multifluid/MultiFluidConstants.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/ComponentProperties.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"

namespace geos
{

namespace constitutive
{
namespace compositional
{

struct NegativeTwoPhaseFlash
{
  using Deriv = constitutive::multifluid::DerivativeOffset;

public:
  /**
   * @brief Perform negative two-phase EOS flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] flashData The parameters required for the flash
   * @param[in/out] kValues The phase equilibrium ratios
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @return an indicator of success of the flash
   */
  template< int USD1, int USD2 >
  GEOS_HOST_DEVICE
  static bool compute( integer const numComps,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition,
                       ComponentProperties::KernelWrapper const & componentProperties,
                       FlashData const & flashData,
                       arraySlice2d< real64, USD1 > const & kValues,
                       real64 & vapourPhaseMoleFraction,
                       arraySlice1d< real64, USD2 > const & liquidComposition,
                       arraySlice1d< real64, USD2 > const & vapourComposition );

  /**
   * @brief Calculate derivatives from the two-phase negative flash
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] flashData The parameters required for the flash
   * @param[in] vapourFraction the calculated vapour (gas) mole fraction
   * @param[in] liquidComposition the calculated liquid phase composition
   * @param[in] vapourComposition the calculated vapour phase composition
   * @param[out] vapourFractionDerivs derivatives of the calculated vapour (gas) mole fraction
   * @param[out] liquidCompositionDerivs derivatives of the calculated liquid phase composition
   * @param[out] vapourCompositionDerivs derivatives of the calculated vapour phase composition
   */
  template< integer USD1, integer USD2, integer USD3 >
  GEOS_HOST_DEVICE
  static void computeDerivatives( integer const numComps,
                                  real64 const pressure,
                                  real64 const temperature,
                                  arraySlice1d< real64 const > const & composition,
                                  ComponentProperties::KernelWrapper const & componentProperties,
                                  FlashData const & flashData,
                                  real64 const & vapourFraction,
                                  arraySlice1d< real64 const, USD1 > const & liquidComposition,
                                  arraySlice1d< real64 const, USD1 > const & vapourComposition,
                                  arraySlice1d< real64, USD2 > const & vapourFractionDerivs,
                                  arraySlice2d< real64, USD3 > const & liquidCompositionDerivs,
                                  arraySlice2d< real64, USD3 > const & vapourCompositionDerivs );

private:
  /**
   * @brief Calculate which components are present.
   * @details Creates a list of indices whose components have non-zero mole fraction.
   * @param[in] numComps number of components
   * @param[in] composition the composition of the fluid
   * @param[out] presentComponents the list of present components
   * @return the number of present components
   */
  template< typename ARRAY >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static integer calculatePresentComponents( integer const numComps,
                                             arraySlice1d< real64 const > const & composition,
                                             ARRAY & presentComponents )
  {
    // Check for machine-zero feed values
    integer presentCount = 0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      if( MultiFluidConstants::minForSpeciesPresence < composition[ic] )
      {
        presentComponents[presentCount++] = ic;
      }
    }
    presentComponents.resize( presentCount );
    return presentCount;
  }

  /**
   * @brief Normalise a composition in place to ensure that the components add up to unity
   * @param[in] numComps number of components
   * @param[in/out] composition composition to be normalized
   * @return the sum of the given values
   */
  template< integer USD >
  GEOS_HOST_DEVICE
  GEOS_FORCE_INLINE
  static real64 normalizeComposition( integer const numComps,
                                      arraySlice1d< real64, USD > const & composition )
  {
    real64 totalMoles = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      totalMoles += composition[ic];
    }
    real64 const oneOverTotalMoles = 1.0 / (totalMoles + MultiFluidConstants::epsilon);
    for( integer ic = 0; ic < numComps; ++ic )
    {
      composition[ic] *= oneOverTotalMoles;
    }
    return totalMoles;
  }

  /**
   * @brief Calculate the logarithms of the fugacity ratios
   * @param[in] numComps number of components
   * @param[in] pressure pressure
   * @param[in] temperature temperature
   * @param[in] composition composition of the mixture
   * @param[in] componentProperties The compositional component properties
   * @param[in] flashData The parameters required for the flash
   * @param[in] kValues The k-values
   * @param[in] presentComponents The indices of the present components
   * @param[out] vapourPhaseMoleFraction the calculated vapour (gas) mole fraction
   * @param[out] liquidComposition the calculated liquid phase composition
   * @param[out] vapourComposition the calculated vapour phase composition
   * @param[out] logLiquidFugacity the calculated log fugacity ratios for the liquid phase
   * @param[out] logVapourFugacity the calculated log fugacity ratios for the vapour phase
   * @param[out] fugacityRatios the fugacity rations
   * @return The error
   */
  template< integer USD1, integer USD2 >
  GEOS_HOST_DEVICE
  static real64 computeFugacityRatio(
    integer const numComps,
    real64 const pressure,
    real64 const temperature,
    arraySlice1d< real64 const > const & composition,
    ComponentProperties::KernelWrapper const & componentProperties,
    FlashData const & flashData,
    arraySlice1d< real64 const, USD1 > const & kValues,
    arraySlice1d< integer const > const & presentComponents,
    real64 & vapourPhaseMoleFraction,
    arraySlice1d< real64, USD2 > const & liquidComposition,
    arraySlice1d< real64, USD2 > const & vapourComposition,
    arraySlice1d< real64 > const & logLiquidFugacity,
    arraySlice1d< real64 > const & logVapourFugacity,
    arraySlice1d< real64 > const & fugacityRatios );

  /**
   * @brief Solve the lineat system for the derivatives of the flash
   * @param[in/out] A the coefficient matrix. Destroyed after call
   * @param[in/out] X the rhs and solution
   * @return @c true if the problem is well solved @c false otherwise
   */
  template< int USD >
  GEOS_HOST_DEVICE
  static bool solveLinearSystem( arraySlice2d< real64, USD > const & A,
                                 arraySlice2d< real64, USD > const & X )
  {
#if defined(GEOS_DEVICE_COMPILE)
    GEOS_UNUSED_VAR( A );
    GEOS_UNUSED_VAR( X );
    return false;
#else
    BlasLapackLA::solveLinearSystem( A, X );
    return true;
#endif
  }

};

template< int USD1, int USD2 >
GEOS_HOST_DEVICE
bool NegativeTwoPhaseFlash::compute( integer const numComps,
                                     real64 const pressure,
                                     real64 const temperature,
                                     arraySlice1d< real64 const > const & composition,
                                     ComponentProperties::KernelWrapper const & componentProperties,
                                     FlashData const & flashData,
                                     arraySlice2d< real64, USD1 > const & kValues,
                                     real64 & vapourPhaseMoleFraction,
                                     arraySlice1d< real64, USD2 > const & liquidComposition,
                                     arraySlice1d< real64, USD2 > const & vapourComposition )
{
  constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  stackArray1d< real64, maxNumComps > logLiquidFugacity( numComps );
  stackArray1d< real64, maxNumComps > logVapourFugacity( numComps );
  stackArray1d< real64, maxNumComps > fugacityRatios( numComps );
  stackArray1d< integer, maxNumComps > componentIndices( numComps );
  auto const & kVapourLiquid = kValues[0];

  calculatePresentComponents( numComps, composition, componentIndices );

  auto const presentComponents = componentIndices.toSliceConst();

  // Initialise compositions to feed composition
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic];
    vapourComposition[ic] = composition[ic];
  }

  // Check if k-Values need to be initialised
  bool needInitialisation = false;
  for( integer ic = 0; ic < numComps; ++ic )
  {
    if( kVapourLiquid[ic] < MultiFluidConstants::epsilon )
    {
      needInitialisation = true;
      break;
    }
  }

  if( needInitialisation )
  {
    if( flashData.liquidEos == EquationOfStateType::SoreideWhitson )
    {
      KValueInitialization::computeSoreideWhitsonKvalue( numComps,
                                                         pressure,
                                                         temperature,
                                                         componentProperties,
                                                         presentComponents,
                                                         kVapourLiquid );
    }
    else
    {
      KValueInitialization::computeWilsonGasLiquidKvalue( numComps,
                                                          pressure,
                                                          temperature,
                                                          componentProperties,
                                                          kVapourLiquid );
    }
  }

  bool converged = false;
  for( localIndex iterationCount = 0; iterationCount < MultiFluidConstants::maxSSIIterations; ++iterationCount )
  {
    real64 const error = computeFugacityRatio( numComps,
                                               pressure,
                                               temperature,
                                               composition,
                                               componentProperties,
                                               flashData,
                                               kVapourLiquid.toSliceConst(),
                                               presentComponents,
                                               vapourPhaseMoleFraction,
                                               liquidComposition,
                                               vapourComposition,
                                               logLiquidFugacity.toSlice(),
                                               logVapourFugacity.toSlice(),
                                               fugacityRatios.toSlice() );

    // Compute fugacity ratios and check convergence
    converged = (error < MultiFluidConstants::fugacityTolerance);

    if( converged )
    {
      break;
    }

    // Update K-values
    for( integer ic = 0; ic < numComps; ++ic )
    {
      kVapourLiquid[ic] *= exp( fugacityRatios[ic] );
    }
  }

  // Retrieve physical bounds from negative flash values
  if( vapourPhaseMoleFraction < MultiFluidConstants::epsilon )
  {
    vapourPhaseMoleFraction = 0.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }
  }
  else if( 1.0 - vapourPhaseMoleFraction < MultiFluidConstants::epsilon )
  {
    vapourPhaseMoleFraction = 1.0;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidComposition[ic] = composition[ic];
      vapourComposition[ic] = composition[ic];
    }
  }

  return converged;
}

template< integer USD1, integer USD2, integer USD3 >
GEOS_HOST_DEVICE
void NegativeTwoPhaseFlash::computeDerivatives(
  integer const numComps,
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const > const & composition,
  ComponentProperties::KernelWrapper const & componentProperties,
  FlashData const & flashData,
  real64 const & vapourFraction,
  arraySlice1d< real64 const, USD1 > const & liquidComposition,
  arraySlice1d< real64 const, USD1 > const & vapourComposition,
  arraySlice1d< real64, USD2 > const & vapourFractionDerivs,
  arraySlice2d< real64, USD3 > const & liquidCompositionDerivs,
  arraySlice2d< real64, USD3 > const & vapourCompositionDerivs )
{
  constexpr integer maxNumComps = MultiFluidConstants::MAX_NUM_COMPONENTS;
  constexpr integer maxNumDofs = MultiFluidConstants::MAX_NUM_COMPONENTS + 2;

  integer const numDofs = numComps + 2;

  auto const setZero = []( real64 & val ) { val = 0.0; };
  LvArray::forValuesInSlice( vapourFractionDerivs, setZero );
  LvArray::forValuesInSlice( liquidCompositionDerivs, setZero );
  LvArray::forValuesInSlice( vapourCompositionDerivs, setZero );

  // Check if we are single or 2-phase
  if( vapourFraction < MultiFluidConstants::epsilon )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
    }
  }
  else if( 1.0 - vapourFraction < MultiFluidConstants::epsilon )
  {
    for( integer ic = 0; ic < numComps; ++ic )
    {
      liquidCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
      vapourCompositionDerivs( ic, Deriv::dC + ic ) = 1.0;
    }
  }
  else
  {
    // Calculate the liquid and vapour fugacities and derivatives
    stackArray1d< real64, maxNumComps > logLiquidFugacity( numComps );
    stackArray1d< real64, maxNumComps > logVapourFugacity( numComps );
    stackArray2d< real64, maxNumComps * maxNumDofs > logLiquidFugacityDerivs( numComps, numDofs );
    stackArray2d< real64, maxNumComps * maxNumDofs > logVapourFugacityDerivs( numComps, numDofs );

    FugacityCalculator::computeLogFugacity( numComps,
                                            pressure,
                                            temperature,
                                            liquidComposition,
                                            componentProperties,
                                            flashData.liquidEos,
                                            flashData,
                                            logLiquidFugacity );
    FugacityCalculator::computeLogFugacity( numComps,
                                            pressure,
                                            temperature,
                                            vapourComposition,
                                            componentProperties,
                                            flashData.vapourEos,
                                            flashData,
                                            logVapourFugacity );

    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       liquidComposition,
                                                       componentProperties,
                                                       flashData.liquidEos,
                                                       flashData,
                                                       logLiquidFugacity.toSliceConst(),
                                                       logLiquidFugacityDerivs.toSlice() );
    FugacityCalculator::computeLogFugacityDerivatives( numComps,
                                                       pressure,
                                                       temperature,
                                                       vapourComposition,
                                                       componentProperties,
                                                       flashData.vapourEos,
                                                       flashData,
                                                       logVapourFugacity.toSliceConst(),
                                                       logVapourFugacityDerivs.toSlice() );

    constexpr integer maxNumVals = 2*MultiFluidConstants::MAX_NUM_COMPONENTS+1;
    integer const numVals = 2*numComps;
    StackArray< real64, 2, maxNumVals * maxNumVals, MatrixLayout::COL_MAJOR_PERM > A( numVals + 1, numVals + 1 );
    StackArray< real64, 2, maxNumVals * maxNumVals, MatrixLayout::COL_MAJOR_PERM > X( numVals + 1, numVals + 1 );

    LvArray::forValuesInSlice( A.toSlice(), setZero );
    LvArray::forValuesInSlice( X.toSlice(), setZero );

    for( integer ic = 0; ic < numComps; ++ic )
    {
      integer const xi = ic;
      integer const yi = ic + numComps;
      integer const vi = numVals;

      integer e = ic;
      A( e, xi ) = 1.0 - vapourFraction;
      A( e, yi ) = vapourFraction;
      A( e, vi ) = vapourComposition[ic] - liquidComposition[ic];

      e = ic + numComps;
      real64 const phiL = exp( logLiquidFugacity( ic ) );
      real64 const phiV = exp( logVapourFugacity( ic ) );
      for( integer jc = 0; jc < numComps; ++jc )
      {
        integer const xj = jc;
        integer const yj = jc + numComps;
        real64 const dPhiLdx = logLiquidFugacityDerivs( ic, Deriv::dC+jc );
        real64 const dPhiVdy = logVapourFugacityDerivs( ic, Deriv::dC+jc );
        A( e, xj ) =  liquidComposition[ic] * phiL * dPhiLdx;
        A( e, yj ) = -vapourComposition[ic] * phiV * dPhiVdy;
      }
      A( e, xi ) += phiL;
      A( e, yi ) -= phiV;

      e = numVals;
      A( e, xi ) = -1.0;
      A( e, yi ) =  1.0;
    }

    // Pressure and temperature derivatives
    for( integer ic = 0; ic < numComps; ++ic )
    {
      real64 const phiL = -liquidComposition[ic] * exp( logLiquidFugacity( ic ) );
      real64 const phiV = vapourComposition[ic] * exp( logVapourFugacity( ic ) );
      X( ic + numComps, Deriv::dP ) = phiL * logLiquidFugacityDerivs( ic, Deriv::dP ) +
                                      phiV * logVapourFugacityDerivs( ic, Deriv::dP );
      X( ic + numComps, Deriv::dT ) = phiL * logLiquidFugacityDerivs( ic, Deriv::dT ) +
                                      phiV * logVapourFugacityDerivs( ic, Deriv::dT );
    }

    // Composition derivatives
    for( integer kc = 0; kc < numComps; ++kc )
    {
      integer const idof = Deriv::dC + kc;

      for( integer ic = 0; ic < numComps; ++ic )
      {
        X( ic, idof ) = -composition[ic];
      }
      X( kc, idof ) += 1.0;
    }

    // Solve linear system
    solveLinearSystem( A.toSlice(), X.toSlice() );

    // Fill in the derivatives
    for( integer idof = 0; idof < numDofs; ++idof )
    {
      for( integer ic = 0; ic < numComps; ++ic )
      {
        liquidCompositionDerivs( ic, idof ) = X( ic, idof );
        vapourCompositionDerivs( ic, idof ) = X( ic + numComps, idof );
      }
      vapourFractionDerivs( idof ) = X( numVals, idof );
    }
  }
}

template< integer USD1, integer USD2 >
GEOS_HOST_DEVICE
real64 NegativeTwoPhaseFlash::computeFugacityRatio(
  integer const numComps,
  real64 const pressure,
  real64 const temperature,
  arraySlice1d< real64 const > const & composition,
  ComponentProperties::KernelWrapper const & componentProperties,
  FlashData const & flashData,
  arraySlice1d< real64 const, USD1 > const & kValues,
  arraySlice1d< integer const > const & presentComponents,
  real64 & vapourPhaseMoleFraction,
  arraySlice1d< real64, USD2 > const & liquidComposition,
  arraySlice1d< real64, USD2 > const & vapourComposition,
  arraySlice1d< real64 > const & logLiquidFugacity,
  arraySlice1d< real64 > const & logVapourFugacity,
  arraySlice1d< real64 > const & fugacityRatios )
{
  // Solve Rachford-Rice Equation
  vapourPhaseMoleFraction = RachfordRice::solve( kValues, composition, presentComponents );

  // Assign phase compositions
  for( integer ic = 0; ic < numComps; ++ic )
  {
    liquidComposition[ic] = composition[ic] / ( 1.0 + vapourPhaseMoleFraction * ( kValues[ic] - 1.0 ) );
    vapourComposition[ic] = kValues[ic] * liquidComposition[ic];
  }

  normalizeComposition( numComps, liquidComposition );
  normalizeComposition( numComps, vapourComposition );

  FugacityCalculator::computeLogFugacity( numComps,
                                          pressure,
                                          temperature,
                                          liquidComposition.toSliceConst(),
                                          componentProperties,
                                          flashData.liquidEos,
                                          flashData,
                                          logLiquidFugacity );
  FugacityCalculator::computeLogFugacity( numComps,
                                          pressure,
                                          temperature,
                                          vapourComposition.toSliceConst(),
                                          componentProperties,
                                          flashData.vapourEos,
                                          flashData,
                                          logVapourFugacity );

  // Compute fugacity ratios and calculate the error
  real64 error = 0.0;
  for( integer const ic : presentComponents )
  {
    fugacityRatios[ic] = ( logLiquidFugacity[ic] - logVapourFugacity[ic] ) + log( liquidComposition[ic] ) - log( vapourComposition[ic] );
    error += (fugacityRatios[ic]*fugacityRatios[ic]);
  }
  return LvArray::math::sqrt( error );
}

} // namespace compositional

} // namespace constitutive

} // namespace geos

#endif //GEOS_CONSTITUTIVE_FLUID_MULTIFLUID_COMPOSITIONAL_FUNCTIONS_NEGATIVETWOPHASEFLASH_HPP_
