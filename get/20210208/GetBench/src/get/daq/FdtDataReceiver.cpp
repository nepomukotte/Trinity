/**
 * @file FdtDataReceiver.cpp
 * @date March 25th, 2013
 * @author sizun
 *
 * @note
 * SVN tag: $Id$
 * Contributors:
 *   - Shebli Anvar (shebli.anvar@cea.fr)
 *   - Frederic Saillant (saillant@ganil.fr)
 *   - Patrick Sizun (patrick.sizun@cea.fr)
 * @note
 * @note This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * @copyright © Grand Accelerateur National d'Ions Lourds (Ganil)
 *
 * @par FREE SOFTWARE LICENCING
 *   This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 *   software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 *   CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 *   and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 *   and the software's author, the holder of the economic rights, and the successive licensors have only limited
 *   liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 *   and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 *   mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 *   experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 *   software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 *   to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 *   you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 * @par COMMERCIAL SOFTWARE LICENCING
 *   You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 *   negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "FdtDataReceiver.h"
#include "mdaq/daq/FdtDataSender.hpp"
#include "mdaq/daq/DataProcessorCore.h"
#include "mdaq/utl/net/TcpAcceptor.h"
#include "mfm/PrimaryHeader.h"
#include "utl/Logging.h"
using mdaq::utl::net::TcpAcceptor;
using mdaq::daq::FdtDataSender;

#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/stream.hpp>
namespace io =  boost::iostreams;
#include <algorithm>

namespace get
{
namespace daq
{
//__________________________________________________________________________________________________
FdtDataReceiver::FdtDataReceiver(const ::utl::net::SocketAddress& sockAddr, size_t rcveBufSize) : sockAddr(sockAddr)
{
	rcveBufSize_ = rcveBufSize;
}
//__________________________________________________________________________________________________
FdtDataReceiver::~FdtDataReceiver()
{
	FdtDataReceiver::stop();
}
//__________________________________________________________________________________________________
void FdtDataReceiver::start()
{
	FdtDataReceiver::stop();
	rcveThread = new RcveThread(*this);
	rcveThread->start();
}
//__________________________________________________________________________________________________
FdtDataReceiver::RcveThread::RcveThread(FdtDataReceiver& receiver) : receiver(receiver),
		data(0), word(0u),
		frameStarted(false), remainingDataSize_B(0u),
		heartBeatReceptionWait_ms(5000u), frameReceptionWait_ms(1000u)
{
	for(int i=1; i <= 10; ++i)
	{
		blips.addArg(std::string(i, '.') + std::string(10-i, ' '));
	}
	connectionWait_ms = DEFAULT_CONNECTION_WAIT;
}
//__________________________________________________________________________________________________
void FdtDataReceiver::askForStop()
{
	if (rcveThread != 0)
	{
		rcveThread->shouldStop = true;
	}
}
//__________________________________________________________________________________________________
void FdtDataReceiver::waitForStop()
{
	if (rcveThread != 0)
	{
		rcveThread->getThreadControl().join();
		rcveThread = 0;
	}
}
//__________________________________________________________________________________________________
void FdtDataReceiver::RcveThread::run()
{
	LOG_DEBUG() << "Entering RcveThread::run()";
	try
	{
		shouldStop = false;
		TcpConnectionPtr connection;
		int iBlip = 0;
		std::auto_ptr< TcpAcceptor > tcpAcceptorPtr(0);
		mdaq::daq::DataProcessorCore& dataProc(receiver.dataProcessorCore());
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
					remainingDataSize_B = 0u;
					frameStarted = false;
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
					dataProc.resetData();

					// Stop accepting new connections
					tcpAcceptorPtr.reset(0);
				}

				// Loop on incoming data
				if (frameStarted)
				{
					while (remainingDataSize_B > 0)
					{
						Milliseconds timeout_ms(frameReceptionWait_ms);
						uint8_t* const dataPtr(&data[0]);
						const size_t dataSize = std::min(remainingDataSize_B, data.size());
						size_t nBytes = connection->tryRcve(dataPtr, dataSize, timeout_ms);
						dataProc.processData(dataPtr, nBytes);
						remainingDataSize_B -= nBytes;
					}
					frameStarted = false;
				}
				else
				{
					// Wait for 4 Byte word
					Milliseconds timeout_ms(heartBeatReceptionWait_ms);
					connection->tryRcve(&word, 4u, timeout_ms);

					// Filter out heart beat
					if (word == FdtDataSender<TcpConnection>::MSG_IMALIVE)
					{
						LOG_INFO() << "Received heart beat...";
					}
					// Handle goodbye
					else if (word == FdtDataSender<TcpConnection>::MSG_GOODBYE)
					{
						LOG_INFO() << "Closing connection (goodbye)...";
						connection.reset();
					}
					else
					{
						frameStarted = true;
						const uint8_t* dataPtr = reinterpret_cast< const uint8_t* >(&word);
						const size_t dataSize(4u);
						const size_t frameSize_B = decodeFrameSize_B(word);
						remainingDataSize_B = frameSize_B - dataSize;
						dataProc.processData(dataPtr, dataSize);
						LOG_INFO() << "Receiving frame of " << frameSize_B << " B";
					}

				}
			}
			catch (const TimeoutException& e)
			{
				if (connection.get() != 0)
				{
					LOG_INFO() << "Closing connection (timeout)...";
					connection.reset();
				}
			}
			catch (const PeerShutdownException & e)
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
//__________________________________________________________________________________________________
size_t FdtDataReceiver::RcveThread::decodeFrameSize_B(const uint32_t & chunk)
{
	const mfm::Byte* buffer = reinterpret_cast< const mfm::Byte* >(&chunk);

	// Construct stream from char array (using boost)
	io::array_source source(buffer, buffer + 4);
	io::stream_buffer< io::array_source > streambuf(source);
	std::istream stream(&streambuf);

	// Decode beginning of primary header
	mfm::PrimaryHeader header;
	try
	{
		header.decode(stream);
	}
	catch (const std::istream::failure & e)
	{
		; // Expected since we only provide four of the height Bytes in a 'primary header'
	}
	return header.frameSize_B();
}
//__________________________________________________________________________________________________
} // namespace daq
} // namespace get
