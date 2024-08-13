/**
 * @file TcpDataReceiver.cpp
 * @created Sep 13, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::daq::TcpDataReceiver
 * Created on: Sep 12, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "TcpDataReceiver.h"
#include "DataProcessorCore.h"
#include "mdaq/utl/net/TcpAcceptor.h"
#include "utl/Logging.h"
using mdaq::utl::net::TcpAcceptor;

namespace mdaq
{
namespace daq
{

TcpDataReceiver::TcpDataReceiver(const ::utl::net::SocketAddress& sockAddr, size_t rcveBufSize) : sockAddr(sockAddr)
{
	rcveBufSize_ = rcveBufSize;
}

TcpDataReceiver::~TcpDataReceiver()
{
	TcpDataReceiver::stop();
}

void TcpDataReceiver::start()
{
	TcpDataReceiver::stop();
	rcveThread = new RcveThread(*this);
	rcveThread->start();
}

TcpDataReceiver::RcveThread::RcveThread(TcpDataReceiver& receiver) : receiver(receiver)
{
	for(int i=1; i <= 10; ++i)
	{
		blips.addArg(std::string(i, '.') + std::string(10-i, ' '));
	}
	connectionWait_ms = DEFAULT_CONNECTION_WAIT;
	receptionWait_ms = DEFAULT_RECEPTION_WAIT;
}

void TcpDataReceiver::askForStop()
{
	if (rcveThread != 0)
	{
		rcveThread->shouldStop = true;
	}
}

void TcpDataReceiver::waitForStop()
{
	if (rcveThread != 0)
	{
		rcveThread->getThreadControl().join();
		rcveThread = 0;
	}
}

void TcpDataReceiver::RcveThread::run()
{
	LOG_DEBUG() << "Entering RcveThread::run()";
	try
	{
		shouldStop = false;
		TcpConnectionPtr connection;
		int iBlip = 0;
		std::auto_ptr< TcpAcceptor > tcpAcceptorPtr(0);
		DataProcessorCore& dataProc(receiver.dataProcessorCore());
		while (not shouldStop)
		{
			try
			{
				// Start accepting connections
				if (connection.get() == 0 and tcpAcceptorPtr.get() == 0)
					tcpAcceptorPtr.reset(new TcpAcceptor(receiver.sockAddr));

				// Wait until you get a connection
				if (connection.get() == 0)
				{
					// No data connection: deactivate data processor
					dataProc.setActive(false);

					// Clear any leftover
					dataProc.resetData();

					// Waiting for client connection
					std::cout << " Waiting for client connection " << blips[iBlip] << '\r' << std::flush;
					iBlip = (iBlip+1) % blips.size();
					connection = tcpAcceptorPtr->acceptConnection(connectionWait_ms);
					LOG_INFO() << "Connection established between " << connection->getPeerAddress() << " and " << connection->getSocketAddress();

					// Setting buffer sizes
					if (receiver.rcveBufSize_ != -1)
					{
						LOG_DEBUG() << "Setting socket receive buffer size to " << receiver.rcveBufSize_;
						connection->set_rcveBufSize(receiver.rcveBufSize_);
					}
					LOG_DEBUG() << "Setting data buffer size to " << std::dec << 4*connection->rcveBufSize() << " bytes";
					data.resize(4*connection->rcveBufSize());

					// now that connection is established, activating data processor
					dataProc.setActive(true);

					// Stop accepting new connections
					tcpAcceptorPtr.reset(0);
				}

				// Loop on incoming data
				const size_t dataSize(data.size());
				uint8_t* const dataPtr(&data[0]);
				while (true)
				{
					Milliseconds timeout_ms(receptionWait_ms); // Reset timeout every time
					size_t nBytes = connection->tryRcve(dataPtr, dataSize, timeout_ms); // Ignore modified timeout value to avoid blocking when zero
					dataProc.processData(dataPtr, nBytes);
				}
			}
			catch(const TimeoutException& e)
			{
				/* Nothing: just see if loop should stop */
			}
			catch(const PeerShutdownException & e)
			{
				LOG_INFO() << e;
				connection.reset();
			}
		}
	}
	catch (const std::exception& e)
	{
		LOG_ERROR() << e.what();
	}
	catch(const std::string& msg)
	{
		LOG_ERROR() << msg;
	}
	catch(const char* msg)
	{
		LOG_ERROR() << msg;
	}
	catch(...)
	{
		LOG_ERROR() << "Unknown exception";
	}
	LOG_DEBUG() << "Leaving RcveThread::run()";
}

} // namespace daq
} // namespace mdaq
