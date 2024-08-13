/**
 * @file Socket.h
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

#ifndef mdaq_utl_net_Socket_h_INCLUDED
#define mdaq_utl_net_Socket_h_INCLUDED

#include <utl/net/SocketAddress.h>
#include <utl/Exception.h>

#if defined(VXWORKS)
#   include <vxWorks.h>
#   include <sockLib.h>
#endif

#if defined(WIN32)
#   include <windows.h>
#   include <winsock2.h>
typedef int socklen_t;
#else
typedef int SOCKET;
#   define INVALID_SOCKET -1
#   define SOCKET_ERROR -1
#endif

#if defined(VXWORKS)
typedef int SocklenType;
typedef char* DataOutType;
typedef const char* DataInType;
#else
typedef socklen_t SocklenType;
typedef void* DataOutType;
typedef const void* DataInType;
#endif

#include <cstring>

namespace mdaq
{
namespace utl
{
namespace net
{
typedef long Milliseconds;
static const Milliseconds INDEFINITE = -1;

using ::utl::net::SocketAddress;

class Socket
{
public:
	class Exception;
	class CreateException;
	class BindException;
	class OptionException;
	class TimeoutException;
	class InvalidException;
	class ErrorException;
	class PeerShutdownException;

	typedef struct timeval Timeval;

public:
	Socket(int bsd_socketType);
	virtual ~Socket();

	bool isValid()  const { return bsd_socket_ != INVALID_SOCKET; }
	SOCKET fd() const { return bsd_socket_; }
	SocketAddress getSocketAddress() const;
	SocketAddress getPeerAddress() const;

	size_t set_sendBufSize(size_t bufSize);
	size_t sendBufSize() const { return cache_sendBufSize; }
	size_t set_rcveBufSize(size_t bufSize);
	size_t rcveBufSize() const { return cache_rcveBufSize; }

	void open();
	void close();

	Milliseconds waitForComm(bool isSending, const Milliseconds& timeout_ms);
	Milliseconds waitForSend(Milliseconds timeout_ms = INDEFINITE) { return waitForComm(true, timeout_ms); }
	Milliseconds waitForRcve(Milliseconds timeout_ms = INDEFINITE) { return waitForComm(false, timeout_ms); }

	template<typename T>
	void setOption(int level, int optname, T optval);

	template<typename T>
	void getOption(int level, int optname, T& optval) const;

	static void setVerbosity(int level) { verbosity = level; }
	static int  getVerbosity() { return verbosity; }

protected:
	Socket(int bsd_socketType, int a_bsd_socket);
	void cache_refresh();
	void bind(const SocketAddress& inboundAddr);

#if defined(__APPLE__)
	static double mach2ms();
#endif

protected:
	int bsd_socketType_;
	SOCKET bsd_socket_;
	size_t cache_sendBufSize;
	size_t cache_rcveBufSize;

	static int verbosity;

private:
	Socket(const Socket& s);
	Socket& operator=(const Socket& s);
};

// =============================================================================
// Nested classes declarations
// -----------------------------------------------------------------------------

class Socket::Exception : public ::utl::Exception
{
public:
	Exception();
	virtual ~Exception() throw();
	virtual void writeMessage(std::ostream& stream) const;
	virtual void writeContext(std::ostream& stream) const = 0;

	const int excep_errno;
};

class Socket::CreateException : public Socket::Exception
{
public:
	CreateException(int bsd_socketType) : bsd_socketType(bsd_socketType) {}
	virtual ~CreateException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const int bsd_socketType;
};

class Socket::BindException : public Socket::Exception
{
public:
	BindException(const int & socket, const SocketAddress& sockAddr) : socket(socket), sockAddr(sockAddr) {}
	virtual ~BindException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const int socket;
	const SocketAddress sockAddr;
};

class Socket::OptionException : public Socket::Exception
{
public:
	OptionException(bool whileSetting, int optName) : whileSetting(whileSetting), optName(optName) {}
	virtual ~OptionException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const bool whileSetting;
	const int optName;
};

class Socket::TimeoutException : public Socket::Exception
{
public:
	TimeoutException(const Milliseconds timeout_ms) : timeout_ms(timeout_ms) {}
	virtual ~TimeoutException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const Milliseconds timeout_ms;
};

class Socket::InvalidException : public Socket::Exception
{
public:
	InvalidException(int bsd_socket) : bsd_socket(bsd_socket) {}
	virtual ~InvalidException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const int bsd_socket;
};

class Socket::ErrorException : public Socket::Exception
{
public:
	ErrorException(int bsd_socket) : bsd_socket(bsd_socket) {}
	virtual ~ErrorException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const int bsd_socket;
};

class Socket::PeerShutdownException : public Socket::Exception
{
public:
	PeerShutdownException(const SocketAddress& sockAddr) : sockAddr(sockAddr) {}
	virtual ~PeerShutdownException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const SocketAddress sockAddr;
};

// =============================================================================
// Inline methods definitions
// -----------------------------------------------------------------------------

template < typename T >
void Socket::setOption(int level, int optname, T optval)
{
#if defined(VXWORKS) || defined(WIN32)
	int st = ::setsockopt(bsd_socket_, level, optname, reinterpret_cast<char*>(&optval), sizeof(T));
#else
	int st = ::setsockopt(bsd_socket_, level, optname, &optval, sizeof(T));
#endif

	if (st < 0)
	{
		throw OptionException(true, optname);
	}
}

template < typename T >
void Socket::getOption(int level, int optname, T& optval) const
{
#if defined(VXWORKS) || defined(WIN32)
	int optSize = sizeof(T);
	int st = ::getsockopt(bsd_socket_, level, optname, reinterpret_cast<char*>(&optval), &optSize);
#else
	socklen_t optSize = sizeof(T);
	int st = ::getsockopt(bsd_socket_, level, optname, &optval, &optSize);
	if(optSize < sizeof(optval))
	{
		if (optSize == sizeof(unsigned int))
			optval = (T) *(reinterpret_cast<unsigned int*>(&optval));
		else if (optSize == sizeof(unsigned long))
			optval = (T) *(reinterpret_cast<unsigned long*>(&optval));
		else if (optSize == sizeof(unsigned short))
			optval = (T) *(reinterpret_cast<unsigned short*>(&optval));
		else if (optSize == sizeof(unsigned char))
			optval = (T) *(reinterpret_cast<unsigned char*>(&optval));
	}
#endif

	if (st < 0)
	{
		throw OptionException(false, optname);
	}
}


} // namespace net
} // namespace utl
} // namespace mdaq

#endif /* mdaq_utl_net_Socket_h_INCLUDED */
