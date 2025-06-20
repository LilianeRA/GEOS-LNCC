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
 * @file SpringSlider.cpp
 */

#include "SpringSlider.hpp"

#include "dataRepository/InputFlags.hpp"
#include "mesh/DomainPartition.hpp"
#include "rateAndStateFields.hpp"
#include "physicsSolvers/solidMechanics/contact/ContactFields.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "constitutive/contact/RateAndStateFriction.hpp"
#include "ExplicitQDRateAndState.hpp"
#include "constitutive/ConstitutivePassThru.hpp"


namespace geos
{

using namespace dataRepository;
using namespace fields;
using namespace constitutive;

template< typename RSSOLVER_TYPE >
SpringSlider< RSSOLVER_TYPE >::SpringSlider( const string & name,
                                             Group * const parent ):
  RSSOLVER_TYPE( name, parent )
{}

template< typename RSSOLVER_TYPE >
SpringSlider< RSSOLVER_TYPE >::~SpringSlider()
{
  // TODO Auto-generated destructor stub
}

template< typename RSSOLVER_TYPE >
void SpringSlider< RSSOLVER_TYPE >::registerDataOnMesh( Group & meshBodies )
{
  RSSOLVER_TYPE::registerDataOnMesh( meshBodies );

  this->forDiscretizationOnMeshTargets( meshBodies, [&] ( string const &,
                                                          MeshLevel & mesh,
                                                          string_array const & regionNames )
  {
    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                 [&]( localIndex const,
                                                                      SurfaceElementSubRegion & subRegion )
    {
      // 3-component functions on fault
      string const labels3Comp[3] = { "normal", "tangent1", "tangent2" };
      subRegion.registerField< contact::dispJump >( this->getName() ).
        setDimLabels( 1, labels3Comp ).
        reference().template resizeDimension< 1 >( 3 );

      subRegion.registerField< contact::dispJump_n >( this->getName() ).
        setDimLabels( 1, labels3Comp ).
        reference().template resizeDimension< 1 >( 3 );

      string const labels2Comp[2] = { "tangent1", "tangent2" };

      subRegion.registerField< contact::deltaSlip >( this->getName() ).
        setDimLabels( 1, labels2Comp ).
        reference().template resizeDimension< 1 >( 2 );

      subRegion.registerField< contact::deltaSlip_n >( this->getName() ).
        setDimLabels( 1, labels2Comp ).
        reference().template resizeDimension< 1 >( 2 );

      subRegion.registerWrapper< string >( viewKeyStruct::frictionLawNameString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE ).
        setSizedFromParent( 0 );

      string & frictionLawName = subRegion.getReference< string >( viewKeyStruct::frictionLawNameString() );
      frictionLawName = PhysicsSolverBase::getConstitutiveName< FrictionBase >( subRegion );
      GEOS_ERROR_IF( frictionLawName.empty(), GEOS_FMT( "{}: FrictionBase model not found on subregion {}",
                                                        this->getDataContext(), subRegion.getDataContext() ) );
    } );
  } );
}

template< typename RSSOLVER_TYPE >
real64 SpringSlider< RSSOLVER_TYPE >::updateStresses( real64 const & time_n,
                                                      real64 const & dt,
                                                      const int cycleNumber,
                                                      DomainPartition & domain ) const

{
  GEOS_UNUSED_VAR( cycleNumber, time_n );
  // Spring-slider shear traction computation
  this->forDiscretizationOnMeshTargets( domain.getMeshBodies(), [&]( string const &,
                                                                     MeshLevel & mesh,
                                                                     string_array const & regionNames )

  {
    mesh.getElemManager().forElementSubRegions< SurfaceElementSubRegion >( regionNames,
                                                                           [&]( localIndex const,
                                                                                SurfaceElementSubRegion & subRegion )
    {



      string const & frictionLawName = subRegion.template getReference< string >( viewKeyStruct::frictionLawNameString() );
      constitutive::ConstitutiveBase & frictionLaw = subRegion.getConstitutiveModel< constitutive::ConstitutiveBase >( frictionLawName );

      constitutive::ConstitutivePassThru< RateAndStateFrictionBase >::execute( frictionLaw, [&] ( auto & castedFrictionLaw )
      {
        typename TYPEOFREF( castedFrictionLaw ) ::KernelWrapper frictionKernelWrapper = castedFrictionLaw.createKernelUpdates();

        updateShearTraction( subRegion, frictionKernelWrapper, dt );
      } );
    } );
  } );
  return dt;
}

template class SpringSlider< ImplicitQDRateAndState >;
template class SpringSlider< ExplicitQDRateAndState >;

namespace
{
typedef SpringSlider< ImplicitQDRateAndState > ImplicitSpringSlider;
typedef SpringSlider< ExplicitQDRateAndState > ExplicitSpringSlider;
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ImplicitSpringSlider, string const &, dataRepository::Group * const )
REGISTER_CATALOG_ENTRY( PhysicsSolverBase, ExplicitSpringSlider, string const &, dataRepository::Group * const )
}

} // namespace geos
