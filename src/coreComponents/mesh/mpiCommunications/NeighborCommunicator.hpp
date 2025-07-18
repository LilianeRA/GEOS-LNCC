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

#ifndef GEOS_MESH_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_
#define GEOS_MESH_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_

#include "common/MpiWrapper.hpp"
#include "mesh/FieldIdentifiers.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "LvArray/src/limits.hpp"
#include "../ElementRegionManager.hpp"

namespace geos
{
inline int CommTag( int const GEOS_UNUSED_PARAM( senderRank ),
                    int const GEOS_UNUSED_PARAM( receiverRank ),
                    int const comm )
{
//  int m_size;
//  MPI_Comm_size( MPI_COMM_GEOS, &m_size );
//  return senderRank * m_size + receiverRank + m_size * m_size * comm;
  return comm;
}

class MeshLevel;
class MPI_iCommData;

class NeighborCommunicator
{
public:
  using ElemAdjListViewType = ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > >;
  using ElemAdjListRefWrapType = ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > >;
  using ElemAdjListRefType = ElementRegionManager::ElementReferenceAccessor< localIndex_array >;

  explicit NeighborCommunicator( int rank );

  NeighborCommunicator():
    NeighborCommunicator( -1 ){};

  void mpiISendReceive( buffer_unit_type const * const sendBuffer,
                        int const sendSize,
                        MPI_Request & sendRequest,
                        buffer_unit_type * const receiveBuffer,
                        int const receiveSize,
                        MPI_Request & receiveRequest,
                        int const commID,
                        MPI_Comm mpiComm );


  void mpiISendReceiveBufferSizes( int const commID,
                                   MPI_Request & mpiSendRequest,
                                   MPI_Request & mpiRecvRequest,
                                   MPI_Comm mpiComm );


  void mpiISendReceiveBuffers( int const commID,
                               MPI_Request & mpiSendRequest,
                               MPI_Request & mpiRecvRequest,
                               MPI_Comm mpiComm );

  template< typename T >
  void mpiISendReceive( T const * const sendBuffer,
                        int const sendSize,
                        MPI_Request & sendReq,
                        T * const recvBuffer,
                        int const recvSize,
                        MPI_Request & recvReq,
                        int const commID,
                        MPI_Comm mpiComm )
  {
    mpiISendReceive( reinterpret_cast< buffer_unit_type const * >( sendBuffer ),
                     sendSize * sizeof(T),
                     sendReq,
                     reinterpret_cast< buffer_unit_type * >( recvBuffer ),
                     recvSize * sizeof(T),
                     recvReq,
                     commID,
                     mpiComm );
  }


  template< typename T >
  void mpiISendReceiveSizes( array1d< T > const & sendBuffer,
                             MPI_Request & sendReq,
                             MPI_Request & recvReq,
                             int const commID,
                             MPI_Comm mpiComm );

  template< typename T >
  void mpiISendReceiveData( array1d< T > const & sendBuffer,
                            MPI_Request & sendReq,
                            array1d< T > & recvBuffer,
                            MPI_Request & recvReq,
                            int const commID,
                            MPI_Comm mpiComm );

  void mpiWaitAll( int const commID,
                   MPI_Request & mpiSendRequest,
                   MPI_Status & mpiSendStatus,
                   MPI_Request & mpiRecvRequest,
                   MPI_Status & mpiReceiveStatus );

  void mpiWaitAll( int const commID );

  /**
   * @brief Post a non-blocking recieve for the buffer size
   *        information for a later call to PostRecv.
   * @note Get the receive MPI_Request by calling
   *       GetSizeRecvRequest with the same CommID
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in
   * @return The return code from the interal call to
   *         MPI_iRecv.
   */
  int postSizeRecv( int const commID,
                    MPI_Request & mpiSizeRecvRequest );

  /**
   * @brief Post a non-blocking send for the buffer
   *        size information so the neighbor knows
   *        how much data it will recieve from a
   *        subsequent call to PostSend.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The return code from the internal call to
   *         MPI_iSend.
   */
  int postSizeSend( int const commID,
                    MPI_Request & mpiSizeSendRequest );

  /**
   * @brief Post a non-blocking receive for the actual
   *        information being recieved from the neighbor.
   * @note Get the receive request by calling GetRecvRequest
   *       with the same commID.
   * @note Make sure the request for any required size recvs
   *       has completed prior to calling this function, as
   *       it uses the size recv information to resize the
   *       internal recv buffer.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The return code from the internal call to
   *         MPI_iRecv.
   */
  int postRecv( int const commID,
                MPI_Request & mpRecvRequest );

  /**
   * @brief Post a non-blocking send for the actual information
   *        being sent to the neighboring rank.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The return code from the internal call to
   *         MPI_iSend.
   */
  int postSend( int const commID,
                MPI_Request & mpiSendRequest );

  /**
   * Posts non-blocking sends to m_neighborRank for
   *  both the size and regular communication buffers
   *  to exchange ghost information.
   *  Additionally posts a non-blocking recv for size
   *  information from m_neighborRank, this size recv
   *  must be completed before PostRecv is called in order
   *  to correctly resize the receive buffer.
   */
  void prepareAndSendGhosts( bool const contactActive,
                             integer const depth,
                             MeshLevel & mesh,
                             int const commID,
                             MPI_Request & mpiRecvSizeRequest,
                             MPI_Request & mpiSendSizeRequest,
                             MPI_Request & mpiSendRequest );

  /**
   * Unpack the receive buffer and process ghosting
   *  information recieved from m_neighborRank.
   *  This must be called after PostRecv is called, and
   *  the request associated with that recv has
   *  completed (retrieve the request using GetRecvRequest)
   */
  void unpackGhosts( MeshLevel & meshLevel,
                     int const commID );

  void unpackGhostsData( MeshLevel & meshLevel,
                         int const commID );

  /**
   * Posts non-blocking sends to m_neighborRank for
   *  both the size and regular communication buffers
   *  to exchange synchronization lists.
   *  Additionally posts a non-blocking recv for size
   *  information from m_neighborRank, this size recv
   *  must be completed before PostRecv is called in order
   *  to correctly resize the receive buffer.
   */
  void prepareAndSendSyncLists( MeshLevel const & mesh,
                                int const commID,
                                MPI_Request & mpiSizeRecvRequest,
                                MPI_Request & mpiSizeSendRequest,
                                MPI_Request & mpiSendRequest );

  /**
   * Unpack the receive buffer and process synchronization
   *  list information received from m_neighborRank.
   *  This must be called after PostRecv is called, and
   *  the request associated with that recv has
   *  completed (retrieve the request using GetRecvRequest)
   */
  void unpackAndRebuildSyncLists( MeshLevel & meshLevel,
                                  int const CommID );

  void packCommBufferForSync( FieldIdentifiers const & fieldsToBeSync,
                              MeshLevel const & meshLevel,
                              int const commID,
                              bool onDevice,
                              parallelDeviceEvents & events );

  int packCommSizeForSync( FieldIdentifiers const & fieldsToBeSync,
                           MeshLevel const & meshLevel,
                           int const commID,
                           bool onDevice,
                           parallelDeviceEvents & events );

  void unpackBufferForSync( FieldIdentifiers const & fieldsToBeSync,
                            MeshLevel & meshLevel,
                            int const commID,
                            bool onDevice,
                            parallelDeviceEvents & events,
                            MPI_Op op=MPI_REPLACE );

  int neighborRank() const { return m_neighborRank; }

  void clear();

  static int constexpr maxComm = 100;

  buffer_type const & receiveBuffer( int commID ) const
  {
    return m_receiveBuffer[commID];
  }
  buffer_type & receiveBuffer( int commID )
  {
    return m_receiveBuffer[commID];
  }

  int const & receiveBufferSize( int commID ) const
  {
    return m_receiveBufferSize[commID];
  }
  int & receiveBufferSize( int commID )
  {
    return m_receiveBufferSize[commID];
  }

  buffer_type const & sendBuffer( int commID ) const
  {
    return m_sendBuffer[commID];
  }
  buffer_type & sendBuffer( int commID )
  {
    return m_sendBuffer[commID];
  }

  void resizeSendBuffer( int const commID, int const newSize )
  {
    m_sendBufferSize[commID] = newSize;
    m_sendBuffer[commID].resize( newSize );
  }

  void resizeRecvBuffer( int const commID, int const newSize )
  {
    m_receiveBufferSize[commID] = newSize;
    m_receiveBuffer[commID].resize( newSize );
  }

  void addNeighborGroupToMesh( MeshLevel & mesh ) const;

private:

  int m_neighborRank;

  int m_sendBufferSize[maxComm];
  int m_receiveBufferSize[maxComm];

  stdVector< buffer_type > m_sendBuffer;
  stdVector< buffer_type > m_receiveBuffer;

  localIndex_array m_nodeUnpackList;
  localIndex_array m_edgeUnpackList;
  localIndex_array m_faceUnpackList;
  ElemAdjListRefType m_elementAdjacencyReceiveListArray;
  buffer_type::size_type m_unpackedSize = 0;
  buffer_unit_type const * m_receiveBufferPtr = nullptr;



};

template< typename T >
void NeighborCommunicator::mpiISendReceiveSizes( array1d< T > const & sendBuffer,
                                                 MPI_Request & sendReq,
                                                 MPI_Request & recvReq,
                                                 int const commID,
                                                 MPI_Comm mpiComm )
{
  m_sendBufferSize[commID] = LvArray::integerConversion< int >( sendBuffer.size());

  mpiISendReceive( &m_sendBufferSize[commID],
                   1,
                   sendReq,
                   &m_receiveBufferSize[commID],
                   1,
                   recvReq,
                   commID,
                   mpiComm );
}

template< typename T >
void NeighborCommunicator::mpiISendReceiveData( array1d< T > const & sendBuffer,
                                                MPI_Request & sendReq,
                                                array1d< T > & recvBuffer,
                                                MPI_Request & recvReq,
                                                int const commID,
                                                MPI_Comm mpiComm )
{
  recvBuffer.resize( m_receiveBufferSize[commID] );

  mpiISendReceive( sendBuffer.data(),
                   m_sendBufferSize[commID],
                   sendReq,
                   recvBuffer.data(),
                   m_receiveBufferSize[commID],
                   recvReq,
                   commID,
                   mpiComm );
}

} /* namespace geos */

#endif /* GEOS_MESH_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_ */
