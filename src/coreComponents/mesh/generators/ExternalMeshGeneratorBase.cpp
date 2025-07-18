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

#include "ExternalMeshGeneratorBase.hpp"

namespace geos
{

using namespace dataRepository;

ExternalMeshGeneratorBase::ExternalMeshGeneratorBase( string const & name,
                                                      dataRepository::Group * const parent )
  : MeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the mesh file" );

  registerWrapper( viewKeyStruct::translateString(), &m_translate ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { 0.0, 0.0, 0.0 } ).
    setDescription( "Translate the coordinates of the vertices by a given vector (prior to scaling)" );

  registerWrapper( viewKeyStruct::scaleString(), &m_scale ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( { 1.0, 1.0, 1.0 } ).
    setDescription( "Scale the coordinates of the vertices by given scale factors (after translation)" );

  registerWrapper( viewKeyStruct::volumicFieldsToImportString(), &m_volumicFieldsToImport ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Volumic fields to be imported from the external mesh file" );

  registerWrapper( viewKeyStruct::volumicFieldsInGEOSString(), &m_volumicFieldsInGEOS ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the volumic fields in GEOS to import into" );

  registerWrapper( viewKeyStruct::surfacicFieldsToImportString(), &m_surfacicFieldsToImport ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Surfacic fields to be imported from the external mesh file" );

  registerWrapper( viewKeyStruct::surfacicFieldsInGEOSString(), &m_surfacicFieldsInGEOS ).
    setRTTypeName( rtTypes::CustomTypes::groupNameRefArray ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Names of the surfacic fields in GEOS to import into" );
}

void ExternalMeshGeneratorBase::postInputInitialization()
{
  auto const checkSizes = [this]( string_array const & from, string_array const & to,
                                  string const & fromKey, string const & toKey )
  {
    GEOS_THROW_IF_NE_MSG( from.size(), to.size(),
                          getWrapperDataContext( fromKey ) <<
                          " and " << getWrapperDataContext( toKey ) <<
                          " must contain the same number of values.",
                          InputError );
  };
  checkSizes( m_volumicFieldsToImport, m_volumicFieldsInGEOS, viewKeyStruct::volumicFieldsToImportString(), viewKeyStruct::volumicFieldsInGEOSString() );
  checkSizes( m_surfacicFieldsToImport, m_surfacicFieldsInGEOS, viewKeyStruct::surfacicFieldsToImportString(), viewKeyStruct::surfacicFieldsInGEOSString() );

  auto const checkDuplicates = [this]( string_array const & v, string const & key )
  {
    std::set< string > const tmp{ v.begin(), v.end() };
    bool const hasDuplicates = tmp.size() != LvArray::integerConversion< std::size_t >( v.size() );

    GEOS_THROW_IF( hasDuplicates,
                   getWrapperDataContext( key ) << ": '" << stringutilities::join( v, ", " ) <<
                   "' already present in list of fields to import.",
                   InputError );
  };
  checkDuplicates( m_volumicFieldsInGEOS, viewKeyStruct::volumicFieldsInGEOSString() );
  checkDuplicates( m_surfacicFieldsInGEOS, viewKeyStruct::surfacicFieldsInGEOSString() );

  // Building the fields mapping from the two separated input/output vectors.
  auto const buildMapping = [&]( string_array const & from,
                                 string_array const & to ) -> std::map< string, string >
  {
    std::map< string, string > mapping;
    for( size_t i = 0; i < from.size(); i++ )
    {
      mapping[from[i]] = to[i];
    }
    return mapping;
  };

  MeshGeneratorBase::m_volumicFields = buildMapping( m_volumicFieldsToImport, m_volumicFieldsInGEOS );
  MeshGeneratorBase::m_surfacicFields = buildMapping( m_surfacicFieldsToImport, m_surfacicFieldsInGEOS );
}

} // namespace geos
