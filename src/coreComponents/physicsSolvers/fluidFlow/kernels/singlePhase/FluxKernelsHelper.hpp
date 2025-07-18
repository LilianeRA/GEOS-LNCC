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
 * @file FluxKernelsHelper.hpp
 */

#ifndef GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUXKERNELSHELPER_HPP
#define GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUXKERNELSHELPER_HPP

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "constitutive/fluid/singlefluid/SingleFluidLayouts.hpp"

namespace geos
{

namespace singlePhaseFluxKernelsHelper
{

/**
 * @brief The type for element-based data. Consists entirely of ArrayView's.
 *
 * Can be converted from ElementRegionManager::ElementViewConstAccessor
 * by calling .toView() or .toViewConst() on an accessor instance
 */
template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

GEOS_HOST_DEVICE
inline
void computeSinglePhaseFlux( localIndex const ( &seri )[2],
                             localIndex const ( &sesri )[2],
                             localIndex const ( &sei )[2],
                             real64 const ( &transmissibility )[2],
                             real64 const ( &dTrans_dPres )[2],
                             ElementViewConst< arrayView1d< real64 const > > const & pres,
                             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                             ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & dens,
                             ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dDens,
                             ElementViewConst< arrayView1d< real64 const > > const & mob,
                             ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & dMob,
                             real64 & alpha,
                             real64 & mobility,
                             real64 & potGrad,
                             real64 & fluxVal,
                             real64 ( & dFlux_dP )[2],
                             real64 & dFlux_dTrans )
{
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 0 >;
  // average density
  real64 densMean = 0.0;
  real64 dDensMean_dP[2];

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    densMean        += 0.5 * dens[seri[ke]][sesri[ke]][sei[ke]][0];
    dDensMean_dP[ke] = 0.5 * dDens[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dP];
  }

  // compute potential difference
  real64 dpotGrad_dTrans = 0.0;
  real64 sumWeightGrav = 0.0;
  real64 potScale = 0.0;
  int signpotGradf[2] = {1, -1};

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    localIndex const er  = seri[ke];
    localIndex const esr = sesri[ke];
    localIndex const ei  = sei[ke];

    real64 const pressure = pres[er][esr][ei];
    real64 const gravD = gravCoef[er][esr][ei];
    real64 const pot = transmissibility[ke] * ( pressure - densMean * gravD );

    potGrad += pot;
    dpotGrad_dTrans += signpotGradf[ke] * ( pressure - densMean * gravD );
    sumWeightGrav += transmissibility[ke] * gravD;

    potScale = fmax( potScale, fabs( pot ) );
  }

  // compute upwinding tolerance
  real64 constexpr upwRelTol = 1e-8;
  real64 const upwAbsTol = fmax( potScale * upwRelTol, LvArray::NumericLimits< real64 >::epsilon );

  // decide mobility coefficients - smooth variation in [-upwAbsTol; upwAbsTol]
  alpha = ( potGrad + upwAbsTol ) / ( 2 * upwAbsTol );

  real64 dMobility_dP[2]{};
  if( alpha <= 0.0 || alpha >= 1.0 )
  {
    // happy path: single upwind direction
    localIndex const ke = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );
    mobility = mob[seri[ke]][sesri[ke]][sei[ke]];
    dMobility_dP[ke] = dMob[seri[ke]][sesri[ke]][sei[ke]][DerivOffset::dP];
  }
  else
  {
    // sad path: weighted averaging
    real64 const mobWeights[2] = { alpha, 1.0 - alpha };
    for( localIndex ke = 0; ke < 2; ++ke )
    {
      mobility += mobWeights[ke] * mob[seri[ke]][sesri[ke]][sei[ke]];
      dMobility_dP[ke] = mobWeights[ke] * dMob[seri[ke]][sesri[ke]][sei[ke]][DerivOffset::dP];
    }
  }

  // compute the final flux and derivative w.r.t transmissibility
  fluxVal = mobility * potGrad;

  dFlux_dTrans = mobility * dpotGrad_dTrans;

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    dFlux_dP[ke] = mobility * ( transmissibility[ke] - dDensMean_dP[ke] * sumWeightGrav )
                   + dMobility_dP[ke] * potGrad + dFlux_dTrans * dTrans_dPres[ke];
  }

}


template< typename ENERGYFLUX_DERIVATIVE_TYPE >
GEOS_HOST_DEVICE
void computeEnthalpyFlux( localIndex const ( &seri )[2],
                          localIndex const ( &sesri )[2],
                          localIndex const ( &sei )[2],
                          real64 const ( &transmissibility )[2],
                          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & enthalpy,
                          ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dEnthalpy,
                          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                          ElementViewConst< arrayView3d< real64 const, constitutive::singlefluid::USD_FLUID_DER > > const & dDens,
                          ElementViewConst< arrayView2d< real64 const, constitutive::singlefluid::USD_FLUID > > const & dMob,
                          real64 const & alpha,
                          real64 const & mobility,
                          real64 const & potGrad,
                          real64 const & massFlux,
                          real64 const & dMassFlux_dTrans,
                          real64 const ( &dMassFlux_dP )[2],
                          real64 ( & dMassFlux_dT )[2],
                          real64 & energyFlux,
                          real64 & dEnergyFlux_dTrans,
                          ENERGYFLUX_DERIVATIVE_TYPE & dEnergyFlux_dP,
                          ENERGYFLUX_DERIVATIVE_TYPE & dEnergyFlux_dT )
{
  // Step 1: compute the derivatives of the mean density at the interface wrt temperature
  using DerivOffset = constitutive::singlefluid::DerivativeOffsetC< 1 >;
  real64 dDensMean_dT[2]{0.0, 0.0};

  for( integer ke = 0; ke < 2; ++ke )
  {
    real64 const dDens_dT = dDens[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dT];
    dDensMean_dT[ke] = 0.5 * dDens_dT;
  }

  // Step 2: compute the derivatives of the potential difference wrt temperature
  //***** calculation of flux *****

  real64 dGravHead_dT[2]{0.0, 0.0};

  // compute potential difference
  for( integer ke = 0; ke < 2; ++ke )
  {
    localIndex const er  = seri[ke];
    localIndex const esr = sesri[ke];
    localIndex const ei  = sei[ke];

    // compute derivative of gravity potential difference wrt temperature
    real64 const gravD = transmissibility[ke] * gravCoef[er][esr][ei];

    for( integer i = 0; i < 2; ++i )
    {
      dGravHead_dT[i] += dDensMean_dT[i] * gravD;
    }
  }

  // Step 3: compute the derivatives of the (upwinded) compFlux wrt temperature
  // *** upwinding ***

  // Step 3.1: compute the derivative of the mass flux wrt temperature
  for( integer ke = 0; ke < 2; ++ke )
  {
    dMassFlux_dT[ke] -= dGravHead_dT[ke];
  }

  for( integer ke = 0; ke < 2; ++ke )
  {
    dMassFlux_dT[ke] *= mobility;
  }

  real64 dMob_dT[2]{};

  if( alpha <= 0.0 || alpha >= 1.0 )
  {
    localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );

    dMob_dT[k_up] = dMob[seri[k_up]][sesri[k_up]][sei[k_up]][DerivOffset::dT];
  }
  else
  {
    real64 const mobWeights[2] = { alpha, 1.0 - alpha };
    for( integer ke = 0; ke < 2; ++ke )
    {
      dMob_dT[ke] = mobWeights[ke] * dMob[seri[ke]][sesri[ke]][sei[ke]][DerivOffset::dT];
    }
  }

  // add contribution from upstream cell mobility derivatives
  for( integer ke = 0; ke < 2; ++ke )
  {
    dMassFlux_dT[ke] += dMob_dT[ke] * potGrad;
  }

  // Step 4: compute the enthalpy flux
  real64 enthalpyTimesMobWeight = 0.0;
  real64 dEnthalpy_dP[2]{0.0, 0.0};
  real64 dEnthalpy_dT[2]{0.0, 0.0};

  if( alpha <= 0.0 || alpha >= 1.0 )
  {
    localIndex const k_up = 1 - localIndex( fmax( fmin( alpha, 1.0 ), 0.0 ) );

    enthalpyTimesMobWeight = enthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0];
    dEnthalpy_dP[k_up] = dEnthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dP];
    dEnthalpy_dT[k_up] = dEnthalpy[seri[k_up]][sesri[k_up]][sei[k_up]][0][DerivOffset::dT];
  }
  else
  {
    real64 const mobWeights[2] = { alpha, 1.0 - alpha };
    for( integer ke = 0; ke < 2; ++ke )
    {
      enthalpyTimesMobWeight += mobWeights[ke] * enthalpy[seri[ke]][sesri[ke]][sei[ke]][0];
      dEnthalpy_dP[ke] = mobWeights[ke] * dEnthalpy[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dP];
      dEnthalpy_dT[ke] = mobWeights[ke] * dEnthalpy[seri[ke]][sesri[ke]][sei[ke]][0][DerivOffset::dT];
    }
  }

  energyFlux += massFlux * enthalpyTimesMobWeight;
  dEnergyFlux_dTrans = enthalpyTimesMobWeight * dMassFlux_dTrans;

  for( integer ke = 0; ke < 2; ++ke )
  {
    dEnergyFlux_dP[ke] += dMassFlux_dP[ke] * enthalpyTimesMobWeight;
    dEnergyFlux_dT[ke] += dMassFlux_dT[ke] * enthalpyTimesMobWeight;
  }

  for( integer ke = 0; ke < 2; ++ke )
  {
    dEnergyFlux_dP[ke] += massFlux * dEnthalpy_dP[ke];
    dEnergyFlux_dT[ke] += massFlux * dEnthalpy_dT[ke];
  }
}


template< typename ENERGYFLUX_DERIVATIVE_TYPE >
GEOS_HOST_DEVICE
void computeConductiveFlux( localIndex const ( &seri )[2],
                            localIndex const ( &sesri )[2],
                            localIndex const ( &sei )[2],
                            ElementViewConst< arrayView1d< real64 const > > const & temperature,
                            real64 const ( &thermalTrans )[2],
                            real64 & energyFlux,
                            ENERGYFLUX_DERIVATIVE_TYPE & dEnergyFlux_dT )
{
  for( integer ke = 0; ke < 2; ++ke )
  {
    localIndex const er  = seri[ke];
    localIndex const esr = sesri[ke];
    localIndex const ei  = sei[ke];

    energyFlux += thermalTrans[ke] * temperature[er][esr][ei];
    dEnergyFlux_dT[ke] += thermalTrans[ke];
  }
}

} // namespace singlePhaseFluxKernelsHelper

} // namespace geos

#endif //GEOS_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASE_FLUXKERNELSHELPER_HPP
