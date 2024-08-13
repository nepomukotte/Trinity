/**
 * @file TcpConnection.cpp
 * @created Aug 8, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::net::TcpConnection
 * Created on: Aug 8, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq software project.
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
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "TcpConnection.h"
#include "utl/Logging.h"
#ifdef VXWORKS

#include <vxWorks.h>
#include <sockLib.h>

#else

#include <sys/socket.h>

#endif

#include <string.h>
#include <stdio.h>
#include <errno.h>

using ::utl::net::SocketAddress;

namespace mdaq
{
namespace utl
{
namespace net
{
/**
 * Creates a TCP connection to a remote server.
 * @param sockAddr Remote address to connect to.
 * @param bindAddr Local address to bind to.
 */
TcpConnection::TcpConnection(const SocketAddress & sockAddr, const ::utl::net::SocketAddress & bindAddr)
	: Socket(SOCK_STREAM)
{
	// Bind a specific local interface (and let the system choose a port)
	Socket::bind(bindAddr);

	// Connect
	sockaddr_in sockaddr;
	::memset(reinterpret_cast<char*>(&sockaddr) , 0, sizeof(sockaddr)); // set all bytes to 0

	sockaddr.sin_family = AF_INET;
	sockaddr.sin_addr.s_addr = sockAddr.isValid() ? sockAddr.ipAddress().in_addr() : INADDR_ANY;
	sockaddr.sin_port = htons(sockAddr.port());

	MSG_IF(verbosity > 0) << "Connecting to " << sockAddr.toString();

	if (::connect(bsd_socket_, (struct sockaddr *) &sockaddr, sizeof(sockaddr)) < 0)
	{
		throw ConnectException(sockAddr);
	}

	setOption<int>(SOL_SOCKET, SO_KEEPALIVE, 1);
}

TcpConnection::TcpConnection(int a_bsd_socket)
	: Socket(SOCK_STREAM, a_bsd_socket)
{
	setOption<int>(SOL_SOCKET, SO_KEEPALIVE, 1);
}

TcpConnection::~TcpConnection()
{
}

#ifdef VXWORKS_KERNEL
TcpConnection::TcpConnection(const TcpConnection&) {}
TcpConnection& TcpConnection::operator=(const TcpConnection&) {}
#endif

size_t TcpConnection::trySend(const void* data, size_t nBytes)
{
	if (nBytes > sendBufSize()) nBytes = sendBufSize(); // Limiting sent data to socket buffer size
	ssize_t bytesSent = ::send(bsd_socket_, reinterpret_cast<DataInType>(data), nBytes, 0);
	if (bytesSent == -1)
	{
		throw SendException(getPeerAddress(), nBytes);
	}
	return static_cast<size_t>(bytesSent);
}

size_t TcpConnection::trySend(const void* data, size_t nBytes, Milliseconds& timeout_ms)
{
	timeout_ms = waitForSend(timeout_ms);
	return trySend(data, nBytes);
}

size_t TcpConnection::tryRcve(void* data, size_t maxBytes)
{
	ssize_t bytesRcved = ::recv(bsd_socket_, reinterpret_cast<DataOutType>(data), maxBytes, 0);
	if (bytesRcved == -1)
	{
		throw RcveException(getSocketAddress());
	}
	else if (bytesRcved == 0)
	{
		throw PeerShutdownException(getSocketAddress());
	}
	return static_cast<size_t>(bytesRcved);
}

size_t TcpConnection::tryRcve(void* data, size_t maxBytes, Milliseconds& timeout_ms)
{
	timeout_ms = waitForRcve(timeout_ms);
	return tryRcve(data, maxBytes);
}

void TcpConnection::doSend(const void* data, size_t nBytes)
{
	const char* dataPtr = reinterpret_cast<const char*>(data);
	while (nBytes > 0u)
	{
		size_t bytesSent = trySend(dataPtr, nBytes);
		nBytes -= bytesSent;
		dataPtr += bytesSent;
	}
}

Milliseconds TcpConnection::doSend(const void* data, size_t nBytes, Milliseconds timeout_ms)
{
	const char* dataPtr = reinterpret_cast<const char*>(data);
	while (nBytes > 0u)
	{
		size_t bytesSent = trySend(dataPtr, nBytes, timeout_ms);
		nBytes -= bytesSent;
		dataPtr += bytesSent;
	}
	return timeout_ms;
}

Milliseconds TcpConnection::doRcve(void* data, size_t nBytes, Milliseconds timeout_ms)
{
	char* dataPtr = reinterpret_cast<char*>(data);
	while(nBytes > 0u)
	{
		size_t bytesRcved = tryRcve(dataPtr, nBytes, timeout_ms);
		nBytes -= bytesRcved;
		dataPtr += bytesRcved;
	}
	return timeout_ms;
}

std::auto_ptr<TcpConnection> TcpConnection::from_bsd_socket(int a_bsd_socket)
{
	std::auto_ptr<TcpConnection> tcpConnection(new TcpConnection(a_bsd_socket));
	return tcpConnection;
}

// =============================================================================
// mdaq::utl::net::TcpConnection::ConnectException class methods body
// -----------------------------------------------------------------------------

TcpConnection::ConnectException::~ConnectException() throw()
{
}

void TcpConnection::ConnectException::writeContext(std::ostream& stream) const
{
	stream << "TCP connection to " << sockAddr.toString() << " failed";
}

// =============================================================================
// mdaq::utl::net::TcpConnection::SendException class methods body
// -----------------------------------------------------------------------------

TcpConnection::SendException::~SendException() throw()
{
}

void TcpConnection::SendException::writeContext(std::ostream& stream) const
{
	stream << "Failed to send "<< std::dec << nBytes << " bytes of data to TCP connection " << sockAddr.toString();
}

// =============================================================================
// mdaq::utl::net::TcpConnection::RcveException class methods body
// -----------------------------------------------------------------------------

TcpConnection::RcveException::~RcveException() throw()
{
}

void TcpConnection::RcveException::writeContext(std::ostream& stream) const
{
	stream << "Failed to receive data from TCP connection " << sockAddr.toString();
}

} // namespace net
} // namespace utl
} // namespace mdaq
