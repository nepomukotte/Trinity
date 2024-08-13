/**
 * @file ZBufConnection.cpp
 * @created Sep 10, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::utl::net::ZBufConnection
 * Created on: Sep 10, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "ZBufConnection.h"
#include "mdaq/mt/SysClock.h"
#include <utl/Logging.h>
#include <errno.h>
#include <vxWorks.h>
#include <netinet/tcp.h>
#include <zbufLib.h>
#include <zbufSockLib.h>

namespace mdaq
{
namespace utl
{
namespace net
{
ZBufConnection::ZBufConnection(const ::utl::net::SocketAddress& sockAddr, const ::utl::net::SocketAddress & bindAddr)
	: TcpConnection(sockAddr, bindAddr), dataSent(semCCreate(SEM_Q_PRIORITY, 0))
{
	int one(1);
	setOption(IPPROTO_TCP, TCP_NODELAY, one);
}

size_t ZBufConnection::trySend(const void* data, size_t nBytes)
{
	if (nBytes > sendBufSize()) nBytes = sendBufSize(); // Limiting sent data to socket buffer size

	ssize_t bytesSent = ::zbufSockBufSend (
		bsd_socket_,
		const_cast<char*>(static_cast<const char*>(data)),
		nBytes,
		(VOIDFUNCPTR) &ZBufConnection::sendCallback,
		reinterpret_cast<int>(this),
		0
	);

	if (ERROR == bytesSent)
	{
		throw SendException(getPeerAddress(), nBytes);
	}

	// Wait for call back to make the semaphore available, checking every 1s the state of the connection.
	const int timeout = mt::SysClock::getResolution_Hz();
	while (ERROR == semTake(dataSent, timeout))
	{
		// If timeout occurred while pending on semaphore, errno=S_objLib_OBJ_TIMEOUT
		// Check connection has not been lost
		const SocketAddress address = getSocketAddress();
		if (ECONNRESET == errno)
			throw PeerShutdownException(address);
	}

	return static_cast<size_t>(bytesSent);
}

size_t ZBufConnection::tryRcve(void* data, size_t maxBytes)
{
	if (maxBytes > sendBufSize()) maxBytes = sendBufSize(); // Limiting received data to socket buffer size
	ZBUF_ID zbRead = ::zbufSockRecv(bsd_socket_, 0, reinterpret_cast<int*>(&maxBytes));
	if (zbRead == NULL)
	{
		throw RcveException(getSocketAddress());
	}
	size_t nBytes = zbufLength(zbRead);
	zbufExtractCopy(zbRead, NULL, 0, static_cast<char*>(data), nBytes);
	zbufDelete(zbRead);
	return nBytes;
}

void ZBufConnection::sendCallback(caddr_t, ZBufConnection* zBufConnect)
{
	if (semGive(zBufConnect->dataSent) == ERROR)
		LOG_ERROR() << "Failed to give dataSent semaphore";
}

} // namespace net
} // namespace utl
} // namespace mdaq
