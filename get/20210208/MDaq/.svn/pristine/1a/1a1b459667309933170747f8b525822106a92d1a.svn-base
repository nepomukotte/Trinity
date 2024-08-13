/* =============================================================================
 * UdpTransceiver.h
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

#ifndef mdaq_utl_net_UdpTransceiver_h_INCLUDED
#define mdaq_utl_net_UdpTransceiver_h_INCLUDED

#include "Socket.h"
#include "Buffer.hpp"

namespace mdaq
{
namespace utl
{
namespace net
{

class UdpTransceiver : public Socket
{
public:
	class SendException;
	class RcveException;

	UdpTransceiver();
	UdpTransceiver(const SocketAddress& inbound);
	virtual ~UdpTransceiver();

	///@{ buffer level methods
	size_t trySend(const SocketAddress& dest, const ConstBuffer& buf);
	size_t trySend(const SocketAddress& dest, const ConstBuffer& buf, Milliseconds& timeout_ms);
	size_t tryRcve(SocketAddress& source, const MutableBuffer& buf);
	size_t tryRcve(SocketAddress& source, const MutableBuffer& buf, Milliseconds& timeout_ms);
	void doSend(const SocketAddress& dest, const ConstBuffer& buf);
	Milliseconds doSend(const SocketAddress& dest, const ConstBuffer& buf, Milliseconds timeout_ms);
	Milliseconds doRcve(SocketAddress& source, const MutableBuffer& buf, Milliseconds timeout_ms = INDEFINITE);
	///@}

	///@{ raw-data level methods
	size_t trySend(const SocketAddress& dest, const void* data, size_t nBytes);
	size_t trySend(const SocketAddress& dest, const void* data, size_t nBytes, Milliseconds& timeout_ms);
	size_t tryRcve(SocketAddress& source, void* data, size_t maxBytes);
	size_t tryRcve(SocketAddress& source, void* data, size_t maxBytes, Milliseconds& timeout_ms);
	void doSend(const SocketAddress& dest, const void* data, size_t nBytes);
	Milliseconds doSend(const SocketAddress& dest, const void* data, size_t nBytes, Milliseconds timeout_ms);
	Milliseconds doRcve(SocketAddress& source, void* data, size_t nBytes, Milliseconds timeout_ms = INDEFINITE);
	///@}

private:
	UdpTransceiver(const UdpTransceiver&);
	UdpTransceiver& operator=(const UdpTransceiver&);
};

// =============================================================================
// Nested classes
// -----------------------------------------------------------------------------

class UdpTransceiver::SendException : public Socket::Exception
{
public:
	SendException(const SocketAddress& sockAddr, size_t nBytes) : sockAddr(sockAddr), nBytes(nBytes) {}
	virtual ~SendException() throw();
	virtual void writeContext(std::ostream& stream) const;

	const SocketAddress sockAddr;
	const size_t nBytes;
};

class UdpTransceiver::RcveException : public Socket::Exception
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

inline size_t UdpTransceiver::trySend(const SocketAddress& dest, const ConstBuffer& buf)
{
	return trySend(dest, buf.data, buf.maxBytes);
}

inline size_t UdpTransceiver::trySend(const SocketAddress& dest, const ConstBuffer& buf, Milliseconds& timeout_ms)
{
	return trySend(dest, buf.data, buf.maxBytes, timeout_ms);
}

inline size_t UdpTransceiver::tryRcve(SocketAddress& source, const MutableBuffer& buf)
{
	return tryRcve(source, buf.data, buf.maxBytes);
}

inline size_t UdpTransceiver::tryRcve(SocketAddress& source, const MutableBuffer& buf, Milliseconds& timeout_ms)
{
	return tryRcve(source, buf.data, buf.maxBytes, timeout_ms);
}

inline void UdpTransceiver::doSend(const SocketAddress& dest, const ConstBuffer& buf)
{
	doSend(dest, buf.data, buf.maxBytes);
}

inline Milliseconds UdpTransceiver::doSend(const SocketAddress& dest, const ConstBuffer& buf, Milliseconds timeout_ms)
{
	return doSend(dest, buf.data, buf.maxBytes, timeout_ms);
}

inline Milliseconds UdpTransceiver::doRcve(SocketAddress& source, const MutableBuffer& buf, Milliseconds timeout_ms)
{
	return doRcve(source, buf.data, buf.maxBytes, timeout_ms);
}


} /* namespace net */
} /* namespace utl */
} /* namespace mdaq */

#endif /* mdaq_utl_net_UdpTransceiver_h_INCLUDED */
