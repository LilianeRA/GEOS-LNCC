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

// Source includes
#include "Group.hpp"
#include "ConduitRestart.hpp"
#include "common/format/StringUtilities.hpp"
#include "common/format/table/TableData.hpp"
#include "common/format/table/TableFormatter.hpp"
#include "common/format/table/TableLayout.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "GroupContext.hpp"
#if defined(GEOS_USE_PYGEOSX)
#include "python/PyGroupType.hpp"
#endif

namespace geos
{
namespace dataRepository
{

Group::Group( string const & name,
              Group * const parent ):
  Group( name, parent->getConduitNode() )
{
  GEOS_ERROR_IF( parent == nullptr, "Should not be null (for Group named " << name << ")." );
  m_parent = parent;
}

Group::Group( string const & name,
              conduit::Node & rootNode ):
  m_parent( nullptr ),
  m_sizedFromParent( 0 ),
  m_wrappers(),
  m_subGroups(),
  m_size( 0 ),
  m_capacity( 0 ),
  m_name( name ),
  m_logLevel( 0 ),
  m_restart_flags( RestartFlags::WRITE_AND_READ ),
  m_input_flags( InputFlags::INVALID ),
  m_conduitNode( rootNode[ name ] ),
  m_logLevelsRegistry( std::make_unique< LogLevelsRegistry >() ),
  m_dataContext( std::make_unique< GroupContext >( *this ) )
{}

Group::~Group()
{
// TODO enable this and fix bugs this exposes.
//  m_conduitNode.parent()->remove( m_name );
}

Group::CatalogInterface::CatalogType & Group::getCatalog()
{
  static Group::CatalogInterface::CatalogType catalog;
  return catalog;
}

WrapperBase & Group::registerWrapper( std::unique_ptr< WrapperBase > wrapper )
{
  // Extract `wrapperName` first to prevent from UB call order in the `insert` call.
  string const wrapperName = wrapper->getName();
  return *m_wrappers.insert( wrapperName, wrapper.release(), true );
}

void Group::deregisterWrapper( string const & name )
{
  GEOS_ERROR_IF( !hasWrapper( name ),
                 "Wrapper " << name << " doesn't exist in Group" << getDataContext() << '.' );
  m_wrappers.erase( name );
  m_conduitNode.remove( name );
}

void Group::resize( indexType const newSize )
{
  forWrappers( [newSize] ( WrapperBase & wrapper )
  {
    if( wrapper.sizedFromParent() == 1 )
    {
      wrapper.resize( newSize );
    }
  } );

  forSubGroups( [newSize] ( Group & subGroup )
  {
    if( subGroup.sizedFromParent() == 1 )
    {
      subGroup.resize( newSize );
    }
  } );

  m_size = newSize;
  if( m_size > m_capacity )
  {
    m_capacity = m_size;
  }
}

void Group::reserve( indexType const newSize )
{
  forWrappers( [newSize] ( WrapperBase & wrapper )
  {
    if( wrapper.sizedFromParent() == 1 )
    {
      wrapper.reserve( newSize );
    }
  } );

  forSubGroups( [newSize] ( Group & subGroup )
  {
    if( subGroup.sizedFromParent() == 1 )
    {
      subGroup.resize( newSize );
    }
  } );

  m_capacity = newSize;
}

string Group::getPath() const
{
  // In the Conduit node hierarchy everything begins with 'Problem', we should change it so that
  // the ProblemManager actually uses the root Conduit Node but that will require a full rebaseline.
  string const noProblem = getConduitNode().path().substr( stringutilities::cstrlen( dataRepository::keys::ProblemManager ) );
  return noProblem.empty() ? "/" : noProblem;
}

string Group::processInputName( xmlWrapper::xmlNode const & targetNode,
                                xmlWrapper::xmlNodePos const & targetNodePos,
                                string_view parentNodeName,
                                xmlWrapper::xmlNodePos const & parentNodePos,
                                std::set< string > & siblingNames )
{
  GEOS_THROW_IF( targetNode.type() != xmlWrapper::xmlNodeType::node_element,
                 GEOS_FMT( "Error in node named \"{}\" ({}): GEOS XML nodes cannot contain "
                           "text data nor anything but XML nodes.\nErroneous content: \"{}\"\n",
                           parentNodeName, parentNodePos.toString(),
                           stringutilities::trimSpaces( targetNode.value() ) ),
                 InputError );

  string targetNodeName;
  try
  { // read & validate the name attribute
    xmlWrapper::readAttributeAsType( targetNodeName, "name",
                                     rtTypes::getTypeRegex< string >( rtTypes::CustomTypes::groupName ),
                                     targetNode, string( "" ) );
  } catch( std::exception const & ex )
  {
    xmlWrapper::processInputException( ex, "name", targetNode, targetNodePos );
  }

  if( targetNodeName.empty() )
  { // if there is no name attribute, we use the node tag as a node name
    targetNodeName = targetNode.name();
  }
  else
  { // Make sure names are not duplicated by checking all previous siblings
    GEOS_THROW_IF( siblingNames.count( targetNodeName ) != 0,
                   GEOS_FMT( "Error at node named \"{}\" ({}): "
                             "An XML block cannot contain children with duplicated names.\n",
                             targetNodeName, targetNodePos.toString() ),
                   InputError );
    siblingNames.insert( targetNodeName );
  }

  return targetNodeName;
}

void Group::processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                       xmlWrapper::xmlNode & targetNode )
{
  xmlWrapper::xmlNodePos targetNodePos = xmlDocument.getNodePosition( targetNode );
  processInputFileRecursive( xmlDocument, targetNode, targetNodePos );
}
void Group::processInputFileRecursive( xmlWrapper::xmlDocument & xmlDocument,
                                       xmlWrapper::xmlNode & targetNode,
                                       xmlWrapper::xmlNodePos const & targetNodePos )
{
  xmlDocument.addIncludedXML( targetNode );

  if( targetNodePos.isFound() )
  {
    m_dataContext = std::make_unique< DataFileContext >( targetNode, targetNodePos );
  }

  // Handle the case where the node was imported from a different input file
  // Set the path prefix to make sure all relative Path variables are interpreted correctly
  string const oldPrefix = std::string( Path::getPathPrefix() );
  xmlWrapper::xmlAttribute filePath = targetNode.attribute( xmlWrapper::filePathString );
  if( filePath )
  {
    Path::setPathPrefix( getAbsolutePath( splitPath( filePath.value() ).first ) );
  }

  // Loop over the child nodes of the targetNode
  std::set< string > childNames;
  for( xmlWrapper::xmlNode childNode : targetNode.children() )
  {
    xmlWrapper::xmlNodePos childNodePos = xmlDocument.getNodePosition( childNode );
    string const childName = processInputName( childNode, childNodePos,
                                               getName(), targetNodePos, childNames );
    // Create children
    Group * newChild = createChild( childNode.name(), childName );
    if( newChild == nullptr )
    {
      newChild = getGroupPointer( childName );
    }
    if( newChild != nullptr )
    {
      newChild->processInputFileRecursive( xmlDocument, childNode, childNodePos );
    }
  }

  processInputFile( targetNode, targetNodePos );

  // Restore original prefix once the node is processed
  Path::setPathPrefix( oldPrefix );
}

void Group::processInputFile( xmlWrapper::xmlNode const & targetNode,
                              xmlWrapper::xmlNodePos const & nodePos )
{
  std::set< string > processedAttributes;
  for( std::pair< string const, WrapperBase * > & pair : m_wrappers )
  {
    if( pair.second->processInputFile( targetNode, nodePos ) )
    {
      processedAttributes.insert( pair.first );
    }
  }

  for( xmlWrapper::xmlAttribute attribute : targetNode.attributes() )
  {
    string const attributeName = attribute.name();
    if( !xmlWrapper::isFileMetadataAttribute( attributeName ) )
    {
      GEOS_THROW_IF( processedAttributes.count( attributeName ) == 0,
                     GEOS_FMT( "Error in {}: XML Node at '{}' contains unused attribute '{}'.\n"
                               "Valid attributes are:\n{}\nFor more details, please refer to documentation at:\n"
                               "http://geosx-geosx.readthedocs-hosted.com/en/latest/docs/sphinx/userGuide/Index.html",
                               getDataContext(), targetNode.path(), attributeName,
                               dumpInputOptions() ),
                     InputError );
    }
  }
}

void Group::postInputInitializationRecursive()
{
  m_logLevelsRegistry = nullptr;
  for( auto const & subGroupIter : m_subGroups )
  {
    subGroupIter.second->postInputInitializationRecursive();
  }
  postInputInitialization();
}

void Group::registerDataOnMeshRecursive( Group & meshBodies )
{
  registerDataOnMesh( meshBodies );
  for( auto && subGroup : m_subGroups )
  {
    subGroup.second->registerDataOnMeshRecursive( meshBodies );
  }
}

Group * Group::createChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0( "Adding Object " << childKey<<" named "<< childName<<" from Group::Catalog." );
  return &registerGroup( childName,
                         CatalogInterface::factory( childKey, getDataContext(),
                                                    childName, this ) );
}

void Group::printDataHierarchy( integer const indent ) const
{
  GEOS_LOG( string( indent, '\t' ) << getName() << " : " << LvArray::system::demangleType( *this ) );
  for( auto & view : wrappers() )
  {
    GEOS_LOG( string( indent, '\t' ) << "-> " << view.second->getName() << " : "
                                     << LvArray::system::demangleType( *view.second ) );
  }
  GEOS_LOG( string( indent, '\t' ) );

  for( auto & group : m_subGroups )
  {
    group.second->printDataHierarchy( indent + 1 );
  }
}

string Group::dumpInputOptions() const
{
  string rval;
  TableLayout const logLayout = TableLayout( "", {TableLayout::Column()
                                                    .setName( "name" )
                                                    .setValuesAlignment( TableLayout::Alignment::left ),
                                                  TableLayout::Column()
                                                    .setName( "Requirement" )
                                                    .setValuesAlignment( TableLayout::Alignment::center ),
                                                  TableLayout::Column()
                                                    .setName( "Description" )
                                                    .setValuesAlignment( TableLayout::Alignment::left ) } )
                                  .setMaxColumnWidth( 80 );
  TableData logData;
  for( auto const & wrapper : m_wrappers )
  {
    if( wrapper.second->getInputFlag() == InputFlags::OPTIONAL ||
        wrapper.second->getInputFlag() == InputFlags::REQUIRED )
    {
      logData.addRow( wrapper.second->getName(),
                      InputFlagToString( wrapper.second->getInputFlag() ),
                      wrapper.second->getDescription() );
    }
  }

  TableTextFormatter logFormatter( logLayout );
  rval.append( logFormatter.toString( logData ));

  return rval;
}

string Group::dumpSubGroupsNames() const
{
  if( numSubGroups() == 0 )
  {
    return getName() + " has no children.";
  }
  else
  {
    return "The children of " + getName() + " are: " +
           "{ " + stringutilities::join( getSubGroupsNames(), ", " ) + " }";
  }
}

string Group::dumpWrappersNames() const
{
  if( numWrappers() == 0 )
  {
    return getName() + " has no wrappers.";
  }
  else
  {
    return "The wrappers of " + getName() + " are: " +
           "{ " + stringutilities::join( getWrappersNames(), ", " ) + " }";
  }
}

void Group::deregisterGroup( string const & name )
{
  GEOS_ERROR_IF( !hasGroup( name ), "Group " << name << " doesn't exist." );
  m_subGroups.erase( name );
  m_conduitNode.remove( name );
}

void Group::initializationOrder( string_array & order )
{
  for( auto & subGroupIter : m_subGroups )
  {
    order.emplace_back( subGroupIter.first );
  }
}

void Group::initialize_postMeshGeneration()
{
  string_array initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initialize_postMeshGeneration();
  }
}


void Group::initialize()
{
  initializePreSubGroups();

  string_array initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initialize();
  }

  initializePostSubGroups();
}


void Group::initializePostInitialConditions()
{
  initializePostInitialConditionsPreSubGroups();

  string_array initOrder;
  initializationOrder( initOrder );

  for( auto const & groupName : initOrder )
  {
    getGroup( groupName ).initializePostInitialConditions();
  }

  initializePostInitialConditionsPostSubGroups();
}

template< bool DO_PACKING >
localIndex Group::packImpl( buffer_unit_type * & buffer,
                            string_array const & wrapperNames,
                            arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  localIndex packedSize = 0;
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, getName() );

  packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "Wrappers" ) );

  // `wrappers` are considered for packing if they match the size of this Group instance.
  // A way to check this is to check the sufficient (but not necessary...) condition `wrapper.sizedFromParent()`.
  stdVector< WrapperBase const * > wrappers;
  for( string const & wrapperName: wrapperNames )
  {
    if( hasWrapper( wrapperName ) )
    {
      WrapperBase const & wrapper = getWrapperBase( wrapperName );

      if( wrapper.sizedFromParent() )
      {
        wrappers.push_back( &wrapper );
      }
    }
    else
    {
      GEOS_ERROR( "Wrapper " << wrapperName << " not found in Group " << getDataContext() << "." );
    }
  }

  // Now we pack the `wrappers`.
  packedSize += bufferOps::Pack< DO_PACKING >( buffer, LvArray::integerConversion< localIndex >( wrappers.size() ) );
  for( WrapperBase const * wrapper: wrappers )
  {
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, wrapper->getName() );
    if( packList.empty() )
    {
      packedSize += wrapper->pack< DO_PACKING >( buffer, true, onDevice, events );
    }
    else
    {
      packedSize += wrapper->packByIndex< DO_PACKING >( buffer, packList, true, onDevice, events );
    }
  }

  if( recursive > 0 )
  {
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, string( "SubGroups" ) );
    packedSize += bufferOps::Pack< DO_PACKING >( buffer, m_subGroups.size() );
    for( auto const & keyGroupPair : m_subGroups )
    {
      packedSize += bufferOps::Pack< DO_PACKING >( buffer, keyGroupPair.first );
      packedSize += keyGroupPair.second->packImpl< DO_PACKING >( buffer, wrapperNames, packList, recursive, onDevice, events );
    }
  }

  return packedSize;
}

localIndex Group::packSize( string_array const & wrapperNames,
                            arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  buffer_unit_type * dummy;
  return this->packImpl< false >( dummy, wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::packSize( arrayView1d< localIndex const > const & packList,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  stdVector< string > const tmp = mapKeys( m_wrappers );
  string_array wrapperNames;
  wrapperNames.insert( wrapperNames.begin(), tmp.begin(), tmp.end() );
  return this->packSize( wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::packSize( string_array const & wrapperNames,
                            integer const recursive,
                            bool onDevice,
                            parallelDeviceEvents & events ) const
{
  arrayView1d< localIndex const > nullArray;
  return packSize( wrapperNames, nullArray, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        string_array const & wrapperNames,
                        arrayView1d< localIndex const > const & packList,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  return this->packImpl< true >( buffer, wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        arrayView1d< localIndex const > const & packList,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  stdVector< string > const tmp = mapKeys( m_wrappers );
  string_array wrapperNames;
  wrapperNames.insert( wrapperNames.begin(), tmp.begin(), tmp.end() );
  return this->pack( buffer, wrapperNames, packList, recursive, onDevice, events );
}


localIndex Group::pack( buffer_unit_type * & buffer,
                        string_array const & wrapperNames,
                        integer const recursive,
                        bool onDevice,
                        parallelDeviceEvents & events ) const
{
  arrayView1d< localIndex const > nullArray;
  return pack( buffer, wrapperNames, nullArray, recursive, onDevice, events );
}

localIndex Group::unpack( buffer_unit_type const * & buffer,
                          arrayView1d< localIndex > & packList,
                          integer const recursive,
                          bool onDevice,
                          parallelDeviceEvents & events,
                          MPI_Op GEOS_UNUSED_PARAM( op ) )
{
  localIndex unpackedSize = 0;
  string groupName;
  unpackedSize += bufferOps::Unpack( buffer, groupName );
  GEOS_ERROR_IF( groupName != getName(), "Group::unpack(): group names do not match" );

  string wrappersLabel;
  unpackedSize += bufferOps::Unpack( buffer, wrappersLabel );
  GEOS_ERROR_IF( wrappersLabel != "Wrappers", "Group::unpack(): wrapper label incorrect" );

  localIndex numWrappers;
  unpackedSize += bufferOps::Unpack( buffer, numWrappers );
  for( localIndex a=0; a<numWrappers; ++a )
  {
    string wrapperName;
    unpackedSize += bufferOps::Unpack( buffer, wrapperName );
    getWrapperBase( wrapperName ).unpackByIndex( buffer, packList, true, onDevice, events );
  }


  if( recursive > 0 )
  {
    string subGroups;
    unpackedSize += bufferOps::Unpack( buffer, subGroups );
    GEOS_ERROR_IF( subGroups != "SubGroups", "Group::unpack(): group names do not match" );

    decltype( m_subGroups.size()) numSubGroups;
    unpackedSize += bufferOps::Unpack( buffer, numSubGroups );
    GEOS_ERROR_IF( numSubGroups != m_subGroups.size(), "Group::unpack(): incorrect number of subGroups" );

    for( auto const & index : m_subGroups )
    {
      GEOS_UNUSED_VAR( index );
      string subGroupName;
      unpackedSize += bufferOps::Unpack( buffer, subGroupName );
      unpackedSize += getGroup( subGroupName ).unpack( buffer, packList, recursive, onDevice, events );
    }
  }

  return unpackedSize;
}


void Group::prepareToWrite()
{
  if( getRestartFlags() == RestartFlags::NO_WRITE )
  {
    return;
  }

  forWrappers( [] ( WrapperBase & wrapper )
  {
    wrapper.registerToWrite();
  } );

  m_conduitNode[ "__size__" ].set( m_size );

  forSubGroups( []( Group & subGroup )
  {
    subGroup.prepareToWrite();
  } );
}


void Group::finishWriting()
{
  if( getRestartFlags() == RestartFlags::NO_WRITE )
  {
    return;
  }

  forWrappers( [] ( WrapperBase & wrapper )
  {
    wrapper.finishWriting();
  } );

  forSubGroups( []( Group & subGroup )
  {
    subGroup.finishWriting();
  } );
}


void Group::loadFromConduit()
{
  if( getRestartFlags() != RestartFlags::WRITE_AND_READ )
  {
    return;
  }

  m_size = m_conduitNode.fetch_existing( "__size__" ).value();
  localIndex const groupSize = m_size;

  forWrappers( [&]( WrapperBase & wrapper )
  {
    if( !( wrapper.loadFromConduit()) )
    {
      if( wrapper.sizedFromParent() == 1 )
      {
        wrapper.resize( groupSize );
      }
    }
  } );

  forSubGroups( []( Group & subGroup )
  {
    subGroup.loadFromConduit();
  } );
}

void Group::postRestartInitializationRecursive()
{
  forSubGroups( [&]( Group & subGroup )
  {
    subGroup.postRestartInitializationRecursive();
  } );

  postRestartInitialization();
}

Group const & Group::getBaseGroupByPath( string const & path ) const
{
  Group const * currentGroup = this;
  string::size_type previousPosition = 0;

  if( path[ 0 ] == '/' )
  {
    bool foundTarget = false;
    for( int i=0; i<1000; ++i )
    {
      if( currentGroup->m_parent != nullptr )
      {
        currentGroup = currentGroup->m_parent;
      }
      else
      {
        foundTarget = true;
        previousPosition = 1;
        break;
      }
    }
    GEOS_THROW_IF( !foundTarget,
                   "Could not find the specified path start.\n"<<
                   "Specified path is " << path,
                   std::domain_error );
  }

  string::size_type currentPosition;
  do
  {
    currentPosition = path.find( '/', previousPosition );
    string const curGroupName = path.substr( previousPosition, currentPosition - previousPosition );

    previousPosition = currentPosition + 1;

    if( curGroupName == "" || curGroupName == "." || curGroupName==currentGroup->m_name )
    {
      continue;
    }
    else if( curGroupName == ".." )
    {
      currentGroup = &this->getParent();
    }
    else
    {
      currentGroup = &currentGroup->getGroup( curGroupName );
    }
  }
  while( currentPosition != string::npos );

  return *currentGroup;
}

localIndex Group::getSubGroupIndex( keyType const & key ) const
{
  return getSubGroups().getIndex( key );
}

#if defined(GEOS_USE_PYGEOSX)
PyTypeObject * Group::getPythonType() const
{ return geos::python::getPyGroupType(); }
#endif

stdVector< string > Group::getSubGroupsNames() const
{
  stdVector< string > childrenNames;
  childrenNames.reserve( numSubGroups() );
  forSubGroups( [&]( Group const & subGroup ){ childrenNames.push_back( subGroup.getName() ); } );
  return childrenNames;
}

stdVector< string > Group::getWrappersNames() const
{
  stdVector< string > wrappersNames;
  wrappersNames.reserve( numWrappers() );
  forWrappers( [&]( WrapperBase const & wrapper ){ wrappersNames.push_back( wrapper.getName() ); } );
  return wrappersNames;
}

} /* end namespace dataRepository */
} /* end namespace geos  */
