/**
 * @file Socket.cpp
 * @created Aug 8, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::net::Socket
 * Created on: Aug 8, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of Mordicus Real-Time Software Development Framework (MDaq).
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

#include "Socket.h"
#include <utl/Logging.h>

#if defined(VXWORKS)
#   include <inetLib.h>
#   include <selectLib.h>
#   include <timers.h>
#elif defined(WIN32)

#else
#   include <unistd.h>
#   include <sys/socket.h>
#   include <sys/time.h>
#   include <netinet/in.h>
#   include <time.h>
#endif

#if defined(__APPLE__)
#   include <mach/mach_time.h>
#endif

#if defined(WIN32)
inline int close(SOCKET s)
{
	return closesocket(s);
}
#endif

#include <errno.h>

namespace mdaq
{
namespace utl
{
namespace net
{
// =============================================================================
// mdaq::utl::net::Socket class methods body
// -----------------------------------------------------------------------------

int Socket::verbosity = 0;

Socket::Socket(int bsd_socketType)
	: bsd_socketType_(bsd_socketType), bsd_socket_(INVALID_SOCKET)
{
	open();
}

Socket::Socket(int bsd_socketType, int a_bsd_socket)
	: bsd_socketType_(bsd_socketType), bsd_socket_(a_bsd_socket)
{
	cache_refresh();
}

Socket::~Socket()
{
	close();
}

#ifdef VXWORKS_KERNEL
Socket::Socket(const Socket&) {}
Socket& Socket::operator=(const Socket&) {}
#endif

void Socket::cache_refresh()
{
	SocklenType bufSize_;
	getOption(SOL_SOCKET, SO_SNDBUF, bufSize_);
	cache_sendBufSize = static_cast<size_t>(bufSize_);
	MSG_IF(verbosity > 0) << "Socket send size: " << std::dec << cache_sendBufSize;

	getOption(SOL_SOCKET, SO_RCVBUF, bufSize_);
	cache_rcveBufSize = static_cast<size_t>(bufSize_);
	MSG_IF(verbosity > 0) << "Socket rcve size: " << std::dec << cache_rcveBufSize;

	getOption(SOL_SOCKET, SO_TYPE, bsd_socketType_);
}

size_t Socket::set_sendBufSize(size_t bufSize)
{
	SocklenType bufSize_ = static_cast<SocklenType>(bufSize);
	setOption(SOL_SOCKET, SO_SNDBUF, bufSize_);
	cache_refresh();
	return sendBufSize();
}

//size_t Socket::sendBufSize() const
//{
//	size_t bufSize;
//	getOption(SOL_SOCKET, SO_SNDBUF, bufSize);
//	return bufSize;
//}

size_t Socket::set_rcveBufSize(size_t bufSize)
{
	SocklenType bufSize_ = static_cast<SocklenType>(bufSize);
	setOption(SOL_SOCKET, SO_RCVBUF, bufSize_);
	cache_refresh();
	return rcveBufSize();
}

//size_t Socket::rcveBufSize() const
//{
//	size_t bufSize;
//	getOption(SOL_SOCKET, SO_RCVBUF, bufSize);
//	return bufSize;
//}

/**
 * Returns the current address to which the socket is bound.
 */
SocketAddress Socket::getSocketAddress() const
{
	SocketAddress address;
	sockaddr_in bsd_sockaddr;
	SocklenType size = sizeof(bsd_sockaddr);
	if(getsockname(bsd_socket_, reinterpret_cast<sockaddr*>(&bsd_sockaddr), &size) == 0)
	{
		address.from_sockaddr_in(bsd_sockaddr);
	}
	return address;
}

/**
 * Returns the address of the peer connected to the socket.
 */
SocketAddress Socket::getPeerAddress() const
{
	SocketAddress address;
	sockaddr_in bsd_sockaddr;
	SocklenType size = sizeof(bsd_sockaddr);
	if(getpeername(bsd_socket_, reinterpret_cast<sockaddr*>(&bsd_sockaddr), &size) == 0)
	{
		address.from_sockaddr_in(bsd_sockaddr);
	}
	return address;
}

void Socket::open()
{
	if (not isValid())
	{
		if ((bsd_socket_ = ::socket(AF_INET, bsd_socketType_, 0)) == SOCKET_ERROR)
		{
			throw CreateException(bsd_socketType_);
		}
		setOption(SOL_SOCKET, SO_REUSEADDR, 1);
	}
	cache_refresh();
}

void Socket::close()
{
	if (isValid())
	{
		::shutdown(bsd_socket_, 2); // shut down both receives and sends
		::close(bsd_socket_);
		bsd_socket_ = INVALID_SOCKET;
	}
}

void Socket::bind(const SocketAddress& inboundAddr)
{
	sockaddr_in bsd_sockaddr;
	::memset(reinterpret_cast<char*>(&bsd_sockaddr) , 0, sizeof(bsd_sockaddr)); // set all bytes to 0

	bsd_sockaddr.sin_family = AF_INET;
	bsd_sockaddr.sin_addr.s_addr = inboundAddr.isValid() ? inboundAddr.ipAddress().in_addr() : INADDR_ANY;
	bsd_sockaddr.sin_port = htons(inboundAddr.port());

    if(::bind(bsd_socket_, reinterpret_cast<sockaddr*>(&bsd_sockaddr), sizeof(bsd_sockaddr)) < 0)
    {
		throw BindException(bsd_socket_, inboundAddr);
    }
}

#if defined(__APPLE__)
/**
 * Returns the conversion factor to convert Mach absolute time into real world time in milliseconds.
 */
double Socket::mach2ms()
{
	static double factor = 0.0;

	if (factor == 0.0)
	{
		mach_timebase_info_data_t info;
		kern_return_t err = mach_timebase_info(&info);

		// Convert the timebase into milliseconds
		if (err == 0)
			factor = 1e-6 * (double) info.numer / (double) info.denom;
	}

	return factor;
}
#endif

Milliseconds Socket::waitForComm(bool isSending, const Milliseconds& timeout_ms)
{
	if (isValid())
	{
		if (timeout_ms == 0) return timeout_ms;

		fd_set descrSet, *readSet = NULL, *writeSet = NULL, errorSet;
		FD_ZERO(&descrSet);
		FD_ZERO(&errorSet);
		FD_SET(bsd_socket_, &descrSet);
		if (isSending) writeSet = &descrSet;
		else readSet = &descrSet;

		int nReady;
		Milliseconds remaining_ms = timeout_ms;
		if (timeout_ms == INDEFINITE)
		{
			nReady = ::select(bsd_socket_+1, readSet, writeSet, &errorSet, NULL);
		}
		else
		{
			Timeval timeout = {0, 0};
			timeout.tv_sec = timeout_ms / 1000;
			timeout.tv_usec = 1000 * (timeout_ms - timeout.tv_sec * 1000);

#if defined(__APPLE__)
			uint64_t startTime, stopTime;
			startTime = mach_absolute_time();
#else
			struct timespec startTime, stopTime;
			clock_gettime(CLOCK_REALTIME, &startTime);
#endif

			nReady = ::select(bsd_socket_+1, readSet, writeSet, NULL, &timeout);

#if defined(__APPLE__)
			stopTime = mach_absolute_time();
			remaining_ms -= (stopTime - startTime)*mach2ms();
#else
			clock_gettime(CLOCK_REALTIME, &stopTime);
			remaining_ms -= (stopTime.tv_sec - startTime.tv_sec)*1000 + (stopTime.tv_nsec - startTime.tv_nsec) /1000000;
#endif
			if (remaining_ms < 0) remaining_ms = 0;
		}
		if (nReady == 0 or not FD_ISSET(bsd_socket_, &descrSet))
		{
			throw TimeoutException(timeout_ms);
		}
		if (FD_ISSET(bsd_socket_, &errorSet))
		{
			throw ErrorException(bsd_socket_);
		}
		return remaining_ms;
	}
	else throw InvalidException(bsd_socket_);
	return timeout_ms;
}

// =============================================================================
// mdaq::utl::net::Socket::Exception class methods body
// -----------------------------------------------------------------------------

Socket::Exception::Exception()
	: excep_errno(errno)
{
}

Socket::Exception::~Exception() throw()
{
}

void Socket::Exception::writeMessage(std::ostream& stream) const
{
	writeContext(stream);
	stream << ": " << strerror(excep_errno) << " (errno=" << excep_errno << ").";
}

// =============================================================================
// mdaq::utl::net::Socket::CreateException class methods body
// -----------------------------------------------------------------------------

Socket::CreateException::~CreateException() throw()
{
}

void Socket::CreateException::writeContext(std::ostream& stream) const
{
	stream << "Failed to create ";
	switch(bsd_socketType)
	{
	case SOCK_STREAM: stream << "TCP"    ; break;
	case SOCK_DGRAM : stream << "UDP"    ; break;
	default         : stream << "UNKNOWN"; break;
	}
	stream << "-type socket";
}

// =============================================================================
// mdaq::utl::net::Socket::OptionException class methods body
// -----------------------------------------------------------------------------

Socket::OptionException::~OptionException() throw()
{
}

void Socket::OptionException::writeContext(std::ostream& stream) const
{
	stream << "Failed to " << (whileSetting ? "set" : "get") << " socket option: ";
	switch(optName)
	{
	case SO_SNDBUF    : stream << "send buffer size"    ; break;
	case SO_RCVBUF    : stream << "receive buffer size" ; break;
	case SO_REUSEADDR : stream << "reuse address"    ; break;
	default           : stream << "UNKNOWN"; break;
	}
}

// =============================================================================
// mdaq::utl::net::Socket::BindException class methods body
// -----------------------------------------------------------------------------

Socket::BindException::~BindException() throw()
{
}

void Socket::BindException::writeContext(std::ostream& stream) const
{
	stream << "Binding socket " << socket << " to " << sockAddr.toString() << " failed";
}

// =============================================================================
// mdaq::utl::net::Socket::TimeoutException class methods body
// -----------------------------------------------------------------------------

Socket::TimeoutException::~TimeoutException() throw()
{
}

void Socket::TimeoutException::writeContext(std::ostream& stream) const
{
	stream << "Timeout of " << timeout_ms << " milliseconds";
}

// =============================================================================
// mdaq::utl::net::Socket::Socket::InvalidException class methods body
// -----------------------------------------------------------------------------

Socket::InvalidException::~InvalidException() throw()
{
}

void Socket::InvalidException::writeContext(std::ostream& stream) const
{
	stream << "Invalid socket: " << bsd_socket;
}

// =============================================================================
// mdaq::utl::net::Socket::ErrorException class methods body
// -----------------------------------------------------------------------------

Socket::ErrorException::~ErrorException() throw()
{
}

void Socket::ErrorException::writeContext(std::ostream& stream) const
{
	stream << "Socket error in socket: " << bsd_socket;
}

// =============================================================================
// mdaq::utl::net::TcpAcceptor class methods body
// -----------------------------------------------------------------------------

Socket::PeerShutdownException::~PeerShutdownException() throw()
{
}

void Socket::PeerShutdownException::writeContext(std::ostream& stream) const
{
	stream << "Connection shutdown by peer of socket " << sockAddr.toString();
}

} // namespace net
} // namespace utl
} // namespace mdaq
