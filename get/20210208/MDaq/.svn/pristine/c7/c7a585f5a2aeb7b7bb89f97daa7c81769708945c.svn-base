/**
 * @file TcpAcceptor.cpp
 * @created Aug 9, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::net::TcpAcceptor
 * Created on: Aug 9, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "TcpAcceptor.h"

using ::utl::net::SocketAddress;

namespace mdaq
{
namespace utl
{
namespace net
{
// =============================================================================
// mdaq::utl::net::TcpAcceptor class methods body
// -----------------------------------------------------------------------------

TcpAcceptor::TcpAcceptor(const SocketAddress& sockAddr) : Socket(SOCK_STREAM)
{
	bind(sockAddr);

    if (::listen(bsd_socket_, TCP_BACKLOG) == -1)
    {
        throw ListenException(sockAddr);
    }
}

TcpAcceptor::~TcpAcceptor()
{
}

#ifdef VXWORKS_KERNEL
TcpAcceptor::TcpAcceptor(const TcpAcceptor&) {}
TcpAcceptor& TcpAcceptor::operator=(const TcpAcceptor&) {}
#endif

std::auto_ptr<TcpConnection>  TcpAcceptor::acceptConnection(Milliseconds timeout_ms)
{
	waitForRcve(timeout_ms);

	sockaddr_in sockaddr_client;
	SocklenType sockaddr_client_len = sizeof(sockaddr_client);
	int bsd_socket_client = ::accept(bsd_socket_,
			reinterpret_cast<sockaddr *>(&sockaddr_client),
			&sockaddr_client_len);
    if (bsd_socket_client < 0)
    {
    	throw AcceptException(getSocketAddress());
    }

    return TcpConnection::from_bsd_socket(bsd_socket_client);
}


// =============================================================================
// mdaq::utl::net::TcpAcceptor::ListenException class methods body
// -----------------------------------------------------------------------------

TcpAcceptor::ListenException::~ListenException() throw()
{
}

void TcpAcceptor::ListenException::writeContext(std::ostream& stream) const
{
	stream << "Listening to " << sockAddr.toString() << " failed";
}

// =============================================================================
// mdaq::utl::net::TcpAcceptor::AcceptException class methods body
// -----------------------------------------------------------------------------

TcpAcceptor::AcceptException::~AcceptException() throw()
{
}

void TcpAcceptor::AcceptException::writeContext(std::ostream& stream) const
{
	stream << "Failed to accept client on " << sockAddr.toString();
}

} // namespace net
} // namespace utl
} // namespace mdaq
