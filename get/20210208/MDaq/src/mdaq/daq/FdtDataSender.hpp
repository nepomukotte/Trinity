/**
 * @file FdtDataSender.hpp
 * @date 25 mars 2013
 * @author sizun
 * @note
 * SVN tag: $Id: FdtDataSender.hpp 1885 2017-05-19 06:49:05Z anvar $
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef mdaq_daq_FdtDataSender_hpp_INCLUDED
#define mdaq_daq_FdtDataSender_hpp_INCLUDED

#include "mdaq/daq/DataSender.h"
#include "mdaq/daq/TcpDataTransmit.h"
#include "mdaq/utl/net/Socket.h"
#ifdef VXWORKS
#include <mdaq/utl/net/ZBufConnection.h>
#endif
#include "utl/net/SocketAddress.h"
#include "utl/Logging.h"

namespace mdaq
{
namespace daq
{
//__________________________________________________________________________________________________
/**
 * A data sender implementing the agreed-upon interface with GANIL Narval actors.
 * It includes:
 *  - periodic emission of heart beats 0x00000000 when there are not data to send.
 *  - message 0xFFFFFFFF sent before disconnection.
 */
#ifdef VXWORKS
template < class T = mdaq::utl::net::ZBufConnection >
#else
template < class T >
#endif
class FdtDataSender : public mdaq::daq::DataSender, public mdaq::daq::TcpDataTransmit
{
public:
	FdtDataSender(const ::utl::net::SocketAddress & sockAddr,
			const ::utl::net::SocketAddress & bindAddr = ::utl::net::SocketAddress(INADDR_ANY, 0));
	virtual ~FdtDataSender();
	virtual void start() { connect(); DataSender::start(); };
	virtual void stop()  { DataSender::stop(); disconnect(); };
	void connect();
	void disconnect();
	virtual void flushData(mdaq::daq::CircularBuffer & circularBuffer);
	virtual void sendData(mdaq::daq::CircularBuffer & circularBuffer);
	virtual void flushChunk(const void* data, uint64_t nBytes);
	virtual void sendHeartBeat();
	virtual std::string senderType() const { return "FDT"; }
public:
	typedef T TcpConnectType; ///< TCP connection type (vanilla TCP or zero-copy buffer).
protected:
	TcpConnection & tcpConnection() { return *tcpConnection_; } ///< Returns TCP connection.
	void sendGoodbye();
protected:
	TcpConnectionPtr tcpConnection_; ///< TCP connection.
	::utl::net::SocketAddress receiverAddress_; ///< Address of the server to connect to.
	::utl::net::SocketAddress senderAddress_; ///< Address of the client to send from.
public:
	static const uint32_t MSG_IMALIVE; ///< Contents of the heart beat message.
	static const uint32_t MSG_GOODBYE; ///< Contents of the disconnection warning message.
};
//__________________________________________________________________________________________________
template < class T > const uint32_t FdtDataSender< T >::MSG_IMALIVE = 0x0ul;
template < class T > const uint32_t FdtDataSender< T >::MSG_GOODBYE = 0xFFFFFFFFul;
//__________________________________________________________________________________________________
template < class T >
FdtDataSender< T >::FdtDataSender(const ::utl::net::SocketAddress & sockAddr,
		const ::utl::net::SocketAddress & bindAddr)
	: receiverAddress_(sockAddr), senderAddress_(bindAddr)
{
}
//__________________________________________________________________________________________________
/**
 * Destroys data sender and all its tasks.
 */
template < class T >
FdtDataSender< T >::~FdtDataSender()
{
	try
	{
		disconnect();
	}
	catch (const PeerShutdownException & e)
	{
	}
}
//__________________________________________________________________________________________________
/**
 * Establishes connection.
 */
template < class T >
void FdtDataSender< T >::connect()
{
	disconnect();
	LOG_INFO() << "Connecting to " << receiverAddress_ << " ...";
	tcpConnection_.reset(new TcpConnectType(receiverAddress_, senderAddress_));
	sendHeartBeat();
}
//__________________________________________________________________________________________________
/**
 * Closes TCP connection after sending GOODBYE message.
 */
template < class T >
void FdtDataSender< T >::disconnect()
{
	if (tcpConnection_.get())
	{
		try
		{
			sendGoodbye();
			tcpConnection_.reset();
		}
		catch (const Socket::Exception & e)
		{
			LOG_ERROR() << "Error closing FDT connection: " << e.what();
		}
	}
}
//__________________________________________________________________________________________________
/**
 * Sends all contents of the circular buffer.
 */
template < class T >
void FdtDataSender<T>::sendData(mdaq::daq::CircularBuffer & circularBuffer)
{
	if (not tcpConnection_.get()) return;

	while (true)
	{
		size_t chunkSize_B = circularBuffer.nextChunkSize_B();
		if (chunkSize_B > 0)
		{
			while (chunkSize_B > 0u)
			{
				size_t bytesSent = tcpConnection().trySend(circularBuffer.readPtr.value(), chunkSize_B);
				chunkSize_B -= bytesSent;
				circularBuffer.readPtr.incr(bytesSent);
			}
		}
		else
		{
			break;
		}
	}
}
//__________________________________________________________________________________________________
/**
 * Sends a sequence of bytes.
 */
template < class T >
void FdtDataSender<T>::flushChunk(const void* data, uint64_t nBytes)
{
	const uint8_t* ptr = reinterpret_cast<const uint8_t*>(data);
	while (nBytes > 0u)
	{
		size_t bytesSent = tcpConnection().trySend(ptr, nBytes);
		nBytes -= bytesSent;
		ptr += bytesSent;
	}
}
//__________________________________________________________________________________________________
/**
 * Sends all contents of the circular buffer.
 */
template < class T >
void FdtDataSender<T>::flushData(mdaq::daq::CircularBuffer & circularBuffer)
{
	sendData(circularBuffer);
}
//__________________________________________________________________________________________________
/**
 * Sends an IMALIVE 32 bits message.
 */
template < class T >
void FdtDataSender< T >::sendHeartBeat()
{
	if (not tcpConnection_.get()) return;

	try
	{
		LOG_DEBUG() << "Sending heart beat";
		tcpConnection().doSend(&MSG_IMALIVE, 4u);
	}
	catch (const Socket::Exception & e)
	{
		LOG_INFO() << "Attempting to reconnect...";
		tcpConnection_.reset();
		connect();
	}
}
//__________________________________________________________________________________________________
/**
 * Sends a GOODBYE 32 bits message.
 */
template < class T >
void FdtDataSender< T >::sendGoodbye()
{
	LOG_INFO() << "Sending GOODBYE";
	tcpConnection().doSend(&MSG_GOODBYE, 4u);
}
//__________________________________________________________________________________________________
} /* namespace daq */
} /* namespace get */
#endif /* mdaq_daq_FdtDataSender_hpp_INCLUDED */
