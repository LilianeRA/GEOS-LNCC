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
 * @file Damage.hpp
 * @brief This class overrides the SSLE constitutive updates to account for a Damage field
 *
 * In a phase-field for fracture model, the damage variable affects the Elasticity equation
 * with the degradation of the stresses. Instead of sigma = C : epsilon, we have sigma = g(d)*C:epsilon,
 * where g(d) is the degradation function. This degradation function can either be a quadratic one
 * (set LORENTZ 0) or a quasi-quadratic one (set LORENTZ 1). In general, the quadratic one will give you
 * brittle fracture behaviour. The quasi-quadratic one, combined with linear dissipation, will give you
 * cohesive fracture behaviour, with a user-defined critical stress. If you use quadratic dissipation in
 * your damage solver, set QUADRATIC_DISSIPATION to 1.
 *
 * References:
 *
 * Miehe, Christian; Hofacker, Martina; Welschinger, Fabian. A phase field model for rate-independent crack
 * propagation: Robust algorithmic implementation based on operator splits.
 * Computer Methods in Applied Mechianics and Engineering, v. 199, n. 45-48, p. 2765-2778, 2010
 *
 * Borden, Micheal J., et al. A phase-field description of dynamic brittle fracture.
 * Computer Methods in Applied Mechanics and Engineering, v. 217, p. 77-95, 2012
 *
 * Bourdin, Blaise; Francfort, Gille A.; Marigo, Jean-Jacques. The variational approach to fracture.
 * Journal of Elasticity, v. 91, n. 1-3, p. 5-148, 2008.
 *
 */

#ifndef GEOS_CONSTITUTIVE_SOLID_DAMAGE_HPP_
#define GEOS_CONSTITUTIVE_SOLID_DAMAGE_HPP_

#include "constitutive/solid/SolidBase.hpp"
#include "InvariantDecompositions.hpp"
#include "ElasticIsotropic.hpp"

namespace geos
{
namespace constitutive
{

// DAMAGE MODEL UPDATES
//
// NOTE: This model uses the m_newStress array to represent the stress in an
//       elastic, "undamaged" configuration.  We then scale the results
//       by the damage factor whenever the true stress is requested through an update
//       function.  The developer should be very cautious if accessing the stress
//       directly through an arrayView, as it does not represent the true stress.
//
// NOTE: This model is designed to work with phase field implementation where m_damage
//       is updated externally to the material model routine.  A modified implementation
//       would be required to use this directly within a SolidMechanics-only solver, to
//       internally update the damage variable and consistently linearize the system.

template< typename UPDATE_BASE >
class DamageUpdates : public UPDATE_BASE
{
public:
  template< typename ... PARAMS >
  DamageUpdates( arrayView2d< real64 > const & inputNewDamage,
                 arrayView2d< real64 > const & inputOldDamage,
                 arrayView3d< real64 > const & inputDamageGrad,
                 arrayView2d< real64 > const & inputStrainEnergyDensity,
                 arrayView2d< real64 > const & inputVolumetricStrain,
                 arrayView2d< real64 > const & inputExtDrivingForce,
                 real64 const & inputLengthScale,
                 arrayView1d< real64 > const & inputCriticalFractureEnergy,
                 real64 const & inputcriticalStrainEnergy,
                 real64 const & inputDegradationLowerLimit,
                 integer const & inputExtDrivingForceFlag,
                 arrayView1d< real64 > const & inputTensileStrength,
                 arrayView1d< real64 > const & inputCompressStrength,
                 arrayView1d< real64 > const & inputDeltaCoefficient,
                 arrayView1d< real64 > const & inputBiotCoefficient,
                 PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_newDamage( inputNewDamage ),
    m_oldDamage( inputOldDamage ),
    m_damageGrad( inputDamageGrad ),
    m_strainEnergyDensity( inputStrainEnergyDensity ),
    m_volStrain( inputVolumetricStrain ),
    m_extDrivingForce ( inputExtDrivingForce ),
    m_lengthScale( inputLengthScale ),
    m_criticalFractureEnergy( inputCriticalFractureEnergy ),
    m_criticalStrainEnergy( inputcriticalStrainEnergy ),
    m_degradationLowerLimit( inputDegradationLowerLimit ),
    m_extDrivingForceFlag( inputExtDrivingForceFlag ),
    m_tensileStrength( inputTensileStrength ),
    m_compressStrength( inputCompressStrength ),
    m_deltaCoefficient( inputDeltaCoefficient ),
    m_biotCoefficient( inputBiotCoefficient )
  {}

  using DiscretizationOps = typename UPDATE_BASE::DiscretizationOps;

  using UPDATE_BASE::smallStrainNoStateUpdate;
  using UPDATE_BASE::smallStrainUpdate;
  using UPDATE_BASE::smallStrainNoStateUpdate_StressOnly;
  using UPDATE_BASE::smallStrainUpdate_StressOnly;
  using UPDATE_BASE::saveConvergedState;

  using UPDATE_BASE::m_disableInelasticity;

  //Standard quadratic degradation functions

  inline
  GEOS_HOST_DEVICE
  virtual real64 getDegradationValue( localIndex const k,
                                      localIndex const q ) const
  {
    real64 pf;

    if( m_extDrivingForceFlag )
    {
      pf = LvArray::math::max( LvArray::math::min( 1.0, m_newDamage( k, q )), 0.0 );
    }
    else
    {
      pf = m_newDamage( k, q );
    }

    // Set a lower bound tolerance for the degradation
    real64 const eps = m_degradationLowerLimit;

    return ((1 - eps)*(1 - pf)*(1 - pf) + eps);
  }


  inline
  GEOS_HOST_DEVICE
  virtual real64 getDegradationDerivative( localIndex const k, real64 const d ) const
  {
    GEOS_UNUSED_VAR( k );

    return -2*(1 - d);
  }


  inline
  GEOS_HOST_DEVICE
  virtual real64 getDegradationSecondDerivative( localIndex const k, real64 const d ) const
  {
    GEOS_UNUSED_VAR( k, d );

    return 2.0;
  }

  //Damage dependence function on fluid pressure terms and its derivatives

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual real64 pressureDamageFunction( localIndex const k,
                                         localIndex const q ) const
  {
    real64 pf = fmax( fmin( 1.0, m_newDamage( k, q )), 0.0 );

    return 0.5*(1 + LvArray::math::cos( M_PI*pf ));
  }


  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual real64 pressureDamageFunctionDerivative( real64 const d ) const
  {
    return -0.5*M_PI*LvArray::math::sin( M_PI*d );
  }


  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual real64 pressureDamageFunctionSecondDerivative( real64 const d ) const
  {
    return -0.5*M_PI*M_PI*LvArray::math::cos( M_PI*d );
  }

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual real64 getDamage( localIndex const k,
                            localIndex const q ) const
  {
    return m_newDamage( k, q );
  }

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual real64 getOldDamage( localIndex const k,
                               localIndex const q ) const
  {
    return m_oldDamage( k, q );
  }

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual void getDamageGrad( localIndex const k,
                              localIndex const q,
                              real64 ( & damageGrad )[3] ) const
  {
    for( int dim=0; dim < 3; ++dim )
    {
      damageGrad[dim] = m_damageGrad[k][q][dim];
    }

  }

  GEOS_FORCE_INLINE
  GEOS_HOST_DEVICE
  virtual void updateBiotCoefficient( localIndex const k,
                                      real64 const biotCoefficient ) const
  {
    m_biotCoefficient[k] = biotCoefficient;
  }

  GEOS_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const & timeIncrement,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const override
  {
    UPDATE_BASE::smallStrainUpdate( k, q, timeIncrement, strainIncrement, stress, stiffness );

    if( m_disableInelasticity )
    {
      return;
    }

    real64 factor = getDegradationValue( k, q );

    // compute volumetric and deviatoric strain invariants
    real64 strain[6] = {0};

    UPDATE_BASE::getElasticStrain( k, q, strain );

    real64 traceOfStrain = strain[0] + strain[1] + strain[2];

    m_volStrain( k, q ) = traceOfStrain;

    if( m_extDrivingForceFlag )
    {
      real64 stressP;
      real64 stressQ;
      real64 deviator[6];

      twoInvariant::stressDecomposition( stress,
                                         stressP,
                                         stressQ,
                                         deviator );

      real64 const mu    = stiffness.m_shearModulus;
      real64 const kappa = stiffness.m_bulkModulus;

      // compute invariants of degraded stress
      real64 I1 = factor * stressP * 3.;
      real64 sqrt_J2 = factor * stressQ / sqrt( 3. );

      real64 const criticalFractureEnergy = m_criticalFractureEnergy[k];
      real64 const tensileStrength = m_tensileStrength[k];
      real64 const compressStrength = m_compressStrength[k];
      real64 const deltaCoefficient = m_deltaCoefficient[k];

      // Calculate the external driving force according to Kumar et al.
      real64 beta0 = deltaCoefficient * 0.375 * criticalFractureEnergy / m_lengthScale;

      real64 beta1 = -0.375 * criticalFractureEnergy / m_lengthScale * ((1 + deltaCoefficient)*(compressStrength - tensileStrength)/2./compressStrength/tensileStrength)
                     - (8*mu + 24*kappa - 27*tensileStrength) * (compressStrength - tensileStrength) / 144. / mu / kappa
                     - m_lengthScale / criticalFractureEnergy * ((mu + 3*kappa)*(pow( compressStrength, 3 ) - pow( tensileStrength, 3 ))*tensileStrength/18/(mu*mu)/(kappa*kappa));

      real64 beta2 = -0.375 * criticalFractureEnergy / m_lengthScale * (sqrt( 3. )*(1 + deltaCoefficient)*(compressStrength + tensileStrength)/2./compressStrength/tensileStrength)
                     + (8*mu + 24*kappa - 27*tensileStrength)*(compressStrength + tensileStrength) / 48. / sqrt( 3. ) / mu / kappa
                     + m_lengthScale / criticalFractureEnergy * ((mu + 3*kappa)*(pow( compressStrength, 3 ) + pow( tensileStrength, 3 ))*tensileStrength/6./sqrt( 3. )/(mu*mu)/(kappa*kappa));

      real64 beta3 = m_lengthScale * (tensileStrength/mu/kappa) / criticalFractureEnergy;

      m_extDrivingForce( k, q ) = 1. / (1 + beta3*I1*I1) * (beta2 * sqrt_J2 + beta1*I1 + beta0);
    }

    LvArray::tensorOps::scale< 6 >( stress, factor );

    stiffness.scaleParams( factor );
  }


  // TODO: The code below assumes the strain energy density will never be
  //       evaluated in a non-converged / garbage configuration.

  GEOS_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override
  {
    real64 const sed = SolidBaseUpdates::getStrainEnergyDensity( k, q );

    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }

    return m_strainEnergyDensity( k, q );
  }

  GEOS_HOST_DEVICE
  virtual real64 getVolStrain( localIndex const k,
                               localIndex const q ) const
  {
    return m_volStrain( k, q );
  }


  GEOS_HOST_DEVICE
  real64 getRegularizationLength() const
  {
    return m_lengthScale;
  }

  GEOS_HOST_DEVICE
  real64 getCriticalFractureEnergy( localIndex const k ) const
  {
    return m_criticalFractureEnergy[k];
  }

  GEOS_HOST_DEVICE
  virtual real64 getEnergyThreshold( localIndex const k,
                                     localIndex const q ) const
  {
    #if LORENTZ
    return m_criticalStrainEnergy;
    #else
    if( m_extDrivingForceFlag )
      return 3*m_criticalFractureEnergy[k]/(16 * m_lengthScale) + 0.5 * m_extDrivingForce( k, q );
    else
      return 3*m_criticalFractureEnergy[k]/(16 * m_lengthScale);

    #endif


  }

  GEOS_HOST_DEVICE
  virtual real64 getBiotCoefficient( localIndex const k ) const
  {
    return m_biotCoefficient[k];
  }

  GEOS_HOST_DEVICE
  virtual void saveConvergedState( localIndex const k,
                                   localIndex const q ) const override final
  {
    ElasticIsotropicUpdates::saveConvergedState( k, q );
    m_oldDamage[k][q] = m_newDamage[k][q];
  }

  /// The new damage value on all quadrature points
  arrayView2d< real64 > const m_newDamage;

  /// The old damage value on all quadrature points
  arrayView2d< real64 > const m_oldDamage;

  /// The damage gradient on all quadrature points
  arrayView3d< real64 > const m_damageGrad;

  /// The strain energy density to drive fracture at the quadrature point
  arrayView2d< real64 > const m_strainEnergyDensity;

  /// The volumetric strain at the quadrature point
  arrayView2d< real64 > const m_volStrain;

  /// The external driving force that accounts for the nucleation of fracture
  arrayView2d< real64 > const m_extDrivingForce;

  /// The phase-field regularization length
  real64 const m_lengthScale;

  /// A reference view to the critical fracture energy for each element
  arrayView1d< real64 > const m_criticalFractureEnergy;

  /// The value of the critical strain energy above which crack/damage initiates
  real64 const m_criticalStrainEnergy;

  /// The lower limit of the degradation function
  real64 const m_degradationLowerLimit;

  /// The flag to indicate if the external driving force is used for fracture nucleation
  integer const m_extDrivingForceFlag;

  /// A reference view to the tensile strength for each element
  arrayView1d< real64 > const m_tensileStrength;

  /// A reference view to the compressive strength for each element
  arrayView1d< real64 > const m_compressStrength;

  /// A reference view to the delta coefficient for each element
  arrayView1d< real64 > const m_deltaCoefficient;

  /// A reference view to the Biot coefficient for each element
  arrayView1d< real64 > const m_biotCoefficient;
};



class DamageBase : public SolidBase
{};

template< typename BASE >
class Damage : public BASE
{
public:

  /// @typedef Alias for LinearElasticIsotropicUpdates
  using KernelWrapper = DamageUpdates< typename BASE::KernelWrapper >;

  Damage( string const & name, dataRepository::Group * const parent );
  virtual ~Damage() override = default;

  static string catalogName() { return string( "Damage" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }

  virtual void postInputInitialization() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void saveConvergedState() const override;
  /// *** The interface to get member variables
  arrayView2d< real64 const > getNewDamage() const { return m_newDamage; }
  arrayView2d< real64 const > getOldDamage() const { return m_oldDamage; }

  arrayView2d< real64 const > getExtDrivingForce() const { return m_extDrivingForce; }


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_newDamage.toView(),
                                                                       m_oldDamage.toView(),
                                                                       m_damageGrad.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_volStrain.toView(),
                                                                       m_extDrivingForce.toView(),
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy.toView(),
                                                                       m_criticalStrainEnergy,
                                                                       m_degradationLowerLimit,
                                                                       m_extDrivingForceFlag,
                                                                       m_tensileStrength.toView(),
                                                                       m_compressStrength.toView(),
                                                                       m_deltaCoefficient.toView(),
                                                                       m_biotCoefficient.toView() );
  }

  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr char const * newDamageString() { return "newDamage"; }
    static constexpr char const * oldDamageString() { return "oldDamage"; }
    static constexpr char const * damageGradString() { return "damageGrad"; }
    static constexpr char const * strainEnergyDensityString() { return "strainEnergyDensity"; }
    static constexpr char const * volumetricStrainString() { return "volumetricStrain"; }
    static constexpr char const * extDrivingForceString() { return "extDrivingForce"; }
    /// string/key for regularization length
    static constexpr char const * lengthScaleString() { return "lengthScale"; }
    /// string/key for default Gc
    static constexpr char const * defaultCriticalFractureEnergyString() { return "defaultCriticalFractureEnergy"; }
    /// string/key for Gc
    static constexpr char const * criticalFractureEnergyString() { return "criticalFractureEnergy"; }
    /// string/key for sigma_c
    static constexpr char const * criticalStrainEnergyString() { return "criticalStrainEnergy"; }
    /// string/key for degradation lower limit
    static constexpr char const * degradationLowerLimitString() { return "degradationLowerLimit"; }
    // string/key for c_e switch
    static constexpr char const * extDrivingForceFlagString() { return "extDrivingForceFlag"; }
    /// string/key for the default tensile strength
    static constexpr char const * defaultTensileStrengthString() { return "defaultTensileStrength"; }
    /// string/key for the uniaxial tensile strength
    static constexpr char const * tensileStrengthString() { return "tensileStrength"; }
    /// string/key for the default compressive strength
    static constexpr char const * defaultCompressStrengthString() { return "defaultCompressiveStrength"; }
    /// string/key for the uniaxial compressive strength
    static constexpr char const * compressStrengthString() { return "compressiveStrength"; }
    /// string/key for the default delta coefficient in computing the external driving force
    static constexpr char const * defaultDeltaCoefficientString() { return "defaultDeltaCoefficient"; }
    /// string/key for a delta coefficient in computing the external driving force
    static constexpr char const * deltaCoefficientString() { return "deltaCoefficient"; }
    /// string/key for the Biot coefficient
    static constexpr char const * biotCoefficientString() { return "biotCoefficient"; }
  };


protected:

  /// The new damage value on all quadrature points
  array2d< real64 > m_newDamage;

  /// The old damage value on all quadrature points
  array2d< real64 > m_oldDamage;

  /// The damage gradient on all quadrature points
  array3d< real64 > m_damageGrad;

  /// The strain energy density to drive fracture at the quadrature point
  array2d< real64 > m_strainEnergyDensity;

  /// The volumetric strain at the quadrature point
  array2d< real64 > m_volStrain;

  /// The external driving force that accounts for the nucleation of fracture
  array2d< real64 > m_extDrivingForce;

  /// The phase-field regularization length
  real64 m_lengthScale;

  /// The default value of critical fracture energy
  real64 m_defaultCriticalFractureEnergy;

  /// The value of the critical strain energy above which crack/damage initiates
  real64 m_criticalStrainEnergy;

  /// The lower limit of the degradation function
  real64 m_degradationLowerLimit;

  /// The flag to indicate if the external driving force is used for fracture nucleation
  integer m_extDrivingForceFlag;

  /// The default value of the tensile strength
  real64 m_defaultTensileStrength;

  /// The default value of the compressive strength
  real64 m_defaultCompressStrength;

  /// The default value of the delta coefficient (should be calibrated)
  real64 m_defaultDeltaCoefficient;

  /// The value of Biot coefficient (to remove it)
  array1d< real64 > m_biotCoefficient;

  /// The critical fracture energy for each cell
  array1d< real64 > m_criticalFractureEnergy;

  /// The tensile strength for each cell
  array1d< real64 > m_tensileStrength;

  /// The compressive strength for each cell
  array1d< real64 > m_compressStrength;

  /// The delta coefficient for each cell
  array1d< real64 > m_deltaCoefficient;
};

}
} /* namespace geos */

#endif /* GEOS_CONSTITUTIVE_SOLID_DAMAGE_HPP_ */
