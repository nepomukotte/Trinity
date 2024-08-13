/**
 * @file TcpConnection.h
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

#ifndef mdaq_utl_net_TcpConnection_h_INCLUDED
#define mdaq_utl_net_TcpConnection_h_INCLUDED

#include "Socket.h"
#include "Buffer.hpp"
#include <utl/BinIO.hpp>
#include <vector>
#include <memory>

namespace mdaq
{
namespace utl
{
namespace net
{

class TcpConnection : public Socket
{
public:
	class ConnectException;
	class SendException;
	class RcveException;

public:
	TcpConnection(const ::utl::net::SocketAddress& sockAddr,
			const ::utl::net::SocketAddress & bindAddr = ::utl::net::SocketAddress(INADDR_ANY, 0));
	virtual ~TcpConnection();

	///@{ buffer level methods
	size_t trySend(const ConstBuffer& buf);
	size_t trySend(const ConstBuffer& buf, Milliseconds& timeout_ms);
	size_t tryRcve(const MutableBuffer& buf);
	size_t tryRcve(const MutableBuffer& buf, Milliseconds& timeout_ms);
	void doSend(const ConstBuffer& buf);
	Milliseconds doSend(const ConstBuffer& buf, Milliseconds timeout_ms);
	Milliseconds doRcve(const MutableBuffer& buf, Milliseconds timeout_ms = INDEFINITE);
	///@}

	///@{ raw-data level methods
	virtual size_t trySend(const void* data, size_t nBytes  );
	        size_t trySend(const void* data, size_t nBytes  , Milliseconds& timeout_ms);
	virtual size_t tryRcve(      void* data, size_t maxBytes);
	        size_t tryRcve(      void* data, size_t maxBytes, Milliseconds& timeout_ms);
	void doSend(const void* data, size_t nBytes);
	Milliseconds doSend(const void* data, size_t nBytes, Milliseconds timeout_ms);
	Milliseconds doRcve(void* data, size_t nBytes, Milliseconds timeout_ms = INDEFINITE);
	///@}

	template < typename T >
	Milliseconds doSendValue(const T& value, Milliseconds timeout_ms = INDEFINITE);
	template < typename T >
	Milliseconds doRcveValue(T& value, Milliseconds timeout_ms = INDEFINITE);

	static std::auto_ptr<TcpConnection> from_bsd_socket(int a_bsd_socket);

protected:
	TcpConnection(int a_bsd_socket);

private:
	TcpConnection(const TcpConnection&);
	TcpConnection& operator=(const TcpConnection&);
};

// =============================================================================
// Nested classes
// -----------------------------------------------------------------------------

class TcpConnection::ConnectException : public Socket::Exception
{
public:
	ConnectException(const SocketAddress& sockAddr) : sockAddr(sockAddr) {}
	virtual ~ConnectException() throw();
	virtual void writeContext(std::ostream& stream) const;
	const SocketAddress sockAddr;
};

class TcpConnection::SendException : public Socket::Exception
{
public:
	SendException(const SocketAddress& sockAddr, size_t nBytes) : sockAddr(sockAddr), nBytes(nBytes) {}
	virtual ~SendException() throw();
	virtual void writeContext(std::ostream& stream) const;
	const SocketAddress sockAddr;
	const size_t nBytes;
};

class TcpConnection::RcveException : public Socket::Exception
{
public:
	RcveException(const SocketAddress& sockAddr) : sockAddr(sockAddr) {}
	virtual ~RcveException() throw();
	virtual void writeContext(std::ostream& stream) const;
	const SocketAddress sockAddr;
};

// =============================================================================
// Inline methods definitions
// -----------------------------------------------------------------------------

inline size_t TcpConnection::trySend(const ConstBuffer& buf)
{
	return trySend(buf.data, buf.maxBytes);
}

inline size_t TcpConnection::trySend(const ConstBuffer& buf, Milliseconds& timeout_ms)
{
	return trySend(buf.data, buf.maxBytes, timeout_ms);
}

inline size_t TcpConnection::tryRcve(const MutableBuffer& buf)
{
	return tryRcve(buf.data, buf.maxBytes);
}

inline size_t TcpConnection::tryRcve(const MutableBuffer& buf, Milliseconds& timeout_ms)
{
	return tryRcve(buf.data, buf.maxBytes, timeout_ms);
}

inline void TcpConnection::doSend(const ConstBuffer& buf)
{
	doSend(buf.data, buf.maxBytes);
}

inline Milliseconds TcpConnection::doSend(const ConstBuffer& buf, Milliseconds timeout_ms)
{
	return doSend(buf.data, buf.maxBytes, timeout_ms);
}

inline Milliseconds TcpConnection::doRcve(const MutableBuffer& buf, Milliseconds timeout_ms)
{
	return doRcve(buf.data, buf.maxBytes, timeout_ms);
}

template < typename T >
Milliseconds TcpConnection::doSendValue(const T& scalar, Milliseconds timeout_ms)
{
	T n(scalar);
	::utl::BinIO::netswap(n); // host to network order byte swap
	timeout_ms = doSend(&n, sizeof(n), timeout_ms);
	return timeout_ms;
}

template < typename T >
Milliseconds TcpConnection::doRcveValue(T& scalar, Milliseconds timeout_ms)
{
	timeout_ms = doRcve(&scalar, sizeof(scalar), timeout_ms);
	::utl::BinIO::netswap(scalar); // network to host order byte swap
	return timeout_ms;
}

template <>
inline Milliseconds TcpConnection::doSendValue<std::string>(const std::string& s, Milliseconds timeout_ms)
{
	uint64_t strSize = s.size();
	timeout_ms = doSendValue(strSize, timeout_ms);
	timeout_ms = doSend(s.data(), s.size(), timeout_ms);
	return timeout_ms;
}

template <>
inline Milliseconds TcpConnection::doRcveValue<std::string>(std::string& s, Milliseconds timeout_ms)
{
	uint64_t strSize;
	timeout_ms = doRcveValue(strSize, timeout_ms);
	s.resize(static_cast<size_t>(strSize));
	for (size_t i=0u; i < strSize; ++i)
	{
		timeout_ms = doRcveValue(s[i], timeout_ms);
	}
	return timeout_ms;
}

} // namespace net
} // namespace utl
} // namespace mdaq

#endif /* mdaq_utl_net_TcpConnection_h_INCLUDED */
