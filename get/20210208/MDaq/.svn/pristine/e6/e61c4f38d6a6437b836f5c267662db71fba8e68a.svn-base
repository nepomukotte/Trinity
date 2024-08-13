/* =============================================================================
 * UdpTransceiver.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::utl::net::UdpTransceiver
 * Created on: 5 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the MDaq2 project.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#include "UdpTransceiver.h"
#include <utl/Logging.h>

namespace mdaq
{
namespace utl
{
namespace net
{

UdpTransceiver::UdpTransceiver()
	: Socket(SOCK_DGRAM)
{
	open();
}

UdpTransceiver::UdpTransceiver(const SocketAddress& inboundAddr)
	: Socket(SOCK_DGRAM)
{
	open();
	bind(inboundAddr);
}

UdpTransceiver::~UdpTransceiver()
{
}

#ifdef VXWORKS_KERNEL
UdpTransceiver::UdpTransceiver(const UdpTransceiver&) {}
UdpTransceiver& UdpTransceiver::operator=(const UdpTransceiver&) {}
#endif

void UdpTransceiver::doSend(const SocketAddress& dest, const void* data, size_t nBytes)
{
	const char* dataPtr = reinterpret_cast<const char*>(data);
	while(nBytes > 0u)
	{
		size_t bytesSent = trySend(dest, dataPtr, nBytes);
		nBytes -= bytesSent;
		dataPtr += bytesSent;
	}
}

Milliseconds UdpTransceiver::doSend(const SocketAddress& dest, const void* data, size_t nBytes, Milliseconds timeout_ms)
{
	const char* dataPtr = reinterpret_cast<const char*>(data);
	while(nBytes > 0u)
	{
		size_t bytesSent = trySend(dest, dataPtr, nBytes, timeout_ms);
		nBytes -= bytesSent;
		dataPtr += bytesSent;
	}
	return timeout_ms;
}

Milliseconds UdpTransceiver::doRcve(SocketAddress& source, void* data, size_t nBytes, Milliseconds timeout_ms)
{
	char* dataPtr = reinterpret_cast<char*>(data);
	while(nBytes > 0u)
	{
		size_t bytesRcved = tryRcve(source, dataPtr, nBytes, timeout_ms);
		nBytes -= bytesRcved;
		dataPtr += bytesRcved;
	}
	return timeout_ms;
}

size_t UdpTransceiver::trySend(const SocketAddress& dest, const void* data, size_t nBytes)
{
	sockaddr_in dest_sockaddr;
	dest_sockaddr.sin_family = AF_INET;
	dest_sockaddr.sin_port = htons(dest.port());
	dest_sockaddr.sin_addr.s_addr = dest.ipAddress().in_addr();

	ssize_t bytesSent =
		::sendto(bsd_socket_,
#ifdef VXWORKS
			const_cast<char*>(reinterpret_cast<DataInType>(data)),
#else
			data,
#endif
			nBytes,
			0,
			reinterpret_cast<sockaddr*>(&dest_sockaddr),
			sizeof(dest_sockaddr));
	if(bytesSent < 0)
	{
		throw SendException(dest, nBytes);
	}
	return bytesSent;
}

size_t UdpTransceiver::trySend(const SocketAddress& dest, const void* data, size_t nBytes, Milliseconds& timeout_ms)
{
	timeout_ms = waitForSend(timeout_ms);
	return trySend(dest, data, nBytes);
}

size_t UdpTransceiver::tryRcve(SocketAddress& source, void* data, size_t maxBytes)
{

	sockaddr_in src_sockaddr;
	SocklenType sockLen = sizeof(src_sockaddr);
	ssize_t bytesRcved =
		::recvfrom(bsd_socket_,
			reinterpret_cast<DataOutType>(data),
			maxBytes,
			0,
			reinterpret_cast<sockaddr*>(&src_sockaddr),
			&sockLen);
	source.from_sockaddr_in(src_sockaddr);
	if(bytesRcved == -1)
	{
		throw RcveException(getSocketAddress());
	}
	else if(bytesRcved == 0)
	{
		throw PeerShutdownException(source);
	}
	return bytesRcved;
}

size_t UdpTransceiver::tryRcve(SocketAddress& source, void* data, size_t maxBytes, Milliseconds& timeout_ms)
{
	timeout_ms = waitForRcve(timeout_ms);
	return tryRcve(source, data, maxBytes);
}

// =============================================================================
// mdaq::utl::net::UdpTransceiver::SendException class methods body
// -----------------------------------------------------------------------------

UdpTransceiver::SendException::~SendException() throw()
{
}

void UdpTransceiver::SendException::writeContext(std::ostream& stream) const
{
	stream << "Failed to send "<< std::dec << nBytes << " bytes of data to UDP endpoint " << sockAddr.toString();
}

// =============================================================================
// mdaq::utl::net::UdpTransceiver::RcveException class methods body
// -----------------------------------------------------------------------------

UdpTransceiver::RcveException::~RcveException() throw()
{
}

void UdpTransceiver::RcveException::writeContext(std::ostream& stream) const
{
	stream << "Failed to receive data from UDP socket " << sockAddr.toString();
}

} /* namespace net */
} /* namespace utl */
} /* namespace mdaq */
