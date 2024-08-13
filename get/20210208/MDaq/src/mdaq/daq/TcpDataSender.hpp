/* =====================================================================================================================
 * TcpDataSender.hpp
 * ---------------------------------------------------------------------------------------------------------------------
 * template class mdaq::daq::TcpDataSender
 * Created on: Apr 4, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * (c) Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * (c) Grand Accelerateur National d'Ions Lourds (Ganil)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Frederic Saillant (saillant@ganil.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mdaq_daq_TcpDataSender_hpp_INCLUDED
#define mdaq_daq_TcpDataSender_hpp_INCLUDED

#include "DataSender.h"
#include "mdaq/daq/TcpDataTransmit.h"

#ifdef VXWORKS
#include <mdaq/utl/net/ZBufConnection.h>
#endif

#include <utl/Logging.h>
#include <utl/BinIO.hpp>
/*
#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
*/
#include <string>
#include <vector>

namespace mdaq
{
namespace daq
{

template < class T >
class TcpDataSender : public DataSender, public TcpDataTransmit
{
public:
	typedef T TcpConnectType;

public:
	TcpDataSender(const ::utl::net::SocketAddress & sockAddr,
			const ::utl::net::SocketAddress & bindAddr = ::utl::net::SocketAddress(INADDR_ANY, 0));
	virtual ~TcpDataSender();
	virtual void sendData(mdaq::daq::CircularBuffer & circularBuffer); ///< @override abstract DataSender::sendData
	virtual void flushData(mdaq::daq::CircularBuffer & circularBuffer);
	virtual void flushChunk(const void* data, uint64_t nBytes);
	TcpConnection& tcpConnection() { return *tcpConnection_; }
	size_t set_sendBufSize(size_t bufSize) { return tcpConnection().set_sendBufSize(bufSize); }
	size_t sendBufSize() { return tcpConnection().sendBufSize(); }
	virtual std::string senderType() const { return "TCP"; }
protected:
	TcpConnectionPtr tcpConnection_;
	const char* nextDataChunkBegin;
	size_t nextDataChunkSize_B;
};

template < class T >
TcpDataSender<T>::TcpDataSender(const ::utl::net::SocketAddress& sockAddr,
		const ::utl::net::SocketAddress & bindAddr)
	: nextDataChunkBegin(0), nextDataChunkSize_B(0u)
{
	tcpConnection_.reset(new TcpConnectType(sockAddr, bindAddr));
}

template < class T >
TcpDataSender<T>::~TcpDataSender()
{
}

template < class T >
void TcpDataSender<T>::sendData(mdaq::daq::CircularBuffer & circularBuffer)
{
	if (alwaysFlushData_)
	{
		flushData(circularBuffer);
	}
	else
	{
		bool bufferWrap = false;
		for (size_t chunkSize_B = circularBuffer.nextChunkSize_B(bufferWrap);
				chunkSize_B > 0;
				chunkSize_B = circularBuffer.nextChunkSize_B(bufferWrap))
		{
			if (chunkSize_B >= sendBufSize())
			{
				tcpConnection().doSend(circularBuffer.readPtr.value(), sendBufSize());
				circularBuffer.readPtr.incr(sendBufSize());
			}
			else if (bufferWrap)
			{
				tcpConnection().doSend(circularBuffer.readPtr.value(), chunkSize_B);
				circularBuffer.readPtr.incr(chunkSize_B);
			}
			else
			{
				break;
			}
		}
	}
}

template < class T >
void TcpDataSender<T>::flushData(mdaq::daq::CircularBuffer & circularBuffer)
{
	while (true)
	{
		size_t chunkSize_B = circularBuffer.nextChunkSize_B();
		if (chunkSize_B > 0)
		{
			/*{
				// Dump data to console
				const char* data = reinterpret_cast<const char*>(circularBuffer.readPtr.value());
				LOG_DEBUG() << "Sending chunk of " << chunkSize_B << " B starting at address " << (void*) data;
				std::ostringstream dump;
				static const size_t maxDump_B = 256; //557568;
				dump << "First " << std::min(chunkSize_B, maxDump_B) << " byte(s) are:";
				dump << std::hex << std::right << std::setfill('0');
				for (size_t i = 0; i < chunkSize_B and i < maxDump_B; ++i)
				{
					if (i % 32 == 0) dump << "\n\t" << std::setw(9) << i << ": ";
					if (i % 8 == 0) dump << ' ';
					dump << std::setw(2) << (int16_t) data[i] << ' ' ;
				}
				LOG_DEBUG() << dump.str();
			}*/
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

template < class T >
void TcpDataSender<T>::flushChunk(const void* data, uint64_t nBytes)
{
	const uint8_t* ptr = reinterpret_cast<const uint8_t*>(data);
	while (nBytes > 0u)
	{
		size_t bytesSent = tcpConnection().trySend(ptr, nBytes);
		nBytes -= bytesSent;
		ptr += bytesSent;
	}
}

typedef TcpDataSender < ::mdaq::utl::net::TcpConnection  >  StdDataSender;

#ifdef VXWORKS
typedef TcpDataSender < ::mdaq::utl::net::ZBufConnection > ZBufDataSender;
#endif

} // namespace daq
} // namespace mdaq

#endif /* mdaq_daq_TcpDataSender_h_INCLUDED */
