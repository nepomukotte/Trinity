/* =====================================================================================================================
 * IceDataSender.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::daq::IceDataSender
 * Created on: Dec 5, 2011 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
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

#include "IceDataSender.h"
#include "mdaq/utl/CmdException.h"
#include "utl/Logging.h"
#include <utility>
#include <iostream>
#include <algorithm>

namespace mdaq
{
namespace daq
{

IceDataSender::IceDataSender(const ::utl::net::SocketAddress& sockAddr)
{
	ic = Ice::initialize();
	init(sockAddr);
}

IceDataSender::IceDataSender(const ::utl::net::SocketAddress& sockAddr, Ice::CommunicatorPtr ic_) : ic(ic_)
{
	init(sockAddr);
}

IceDataSender::~IceDataSender()
{
}

void IceDataSender::init(const ::utl::net::SocketAddress& sockAddr)
{
	LOG_DEBUG() << "Allocating IceDataSender over proxy: " << sockAddr.iceProxyString("DataProcessor");
	Ice::ObjectPrx objPrx = ic->stringToProxy(sockAddr.iceProxyString("DataProcessor"))->ice_timeout(10000);
	dataProcessorPrx = DataProcessorPrx::checkedCast(objPrx); // checking object existence
	dataProcessorPrx = DataProcessorPrx::uncheckedCast(dataProcessorPrx->ice_oneway()); // keeping one-way proxy

	// Get max Ice message size and subtract 512 B for size of message header.
	maxDataChunkSize_B_ = 1024*ic->getProperties()->getPropertyAsIntWithDefault("Ice.MessageSizeMax", 1024) - 512;
	LOG_DEBUG() << "Max. data chunk size: " << maxDataChunkSize_B_ << " B";

	LOG_DEBUG() << "IceDataSender constructed";
}

void IceDataSender::sendData(mdaq::daq::CircularBuffer & circularBuffer)
{
	for (uint64_t chunkSize_B = std::min(circularBuffer.nextChunkSize_B(), maxDataChunkSize_B_);
			chunkSize_B > 0;
			chunkSize_B = std::min(circularBuffer.nextChunkSize_B(), maxDataChunkSize_B_))
	{
		const Ice::Byte* buffer = reinterpret_cast<const Ice::Byte*>(circularBuffer.readPtr.value());
		dataProcessorPrx->processData(std::make_pair(buffer, buffer + chunkSize_B));
		circularBuffer.readPtr.incr(chunkSize_B);
	}
}

void IceDataSender::flushData(mdaq::daq::CircularBuffer & circularBuffer)
{
	sendData(circularBuffer);
}

void IceDataSender::flushChunk(const void* data, uint64_t nBytes)
{
	const Ice::Byte* buffer = reinterpret_cast<const Ice::Byte*>(data);
	while (nBytes > 0u)
	{
		size_t bytesSent = std::min(nBytes, maxDataChunkSize_B_);
		dataProcessorPrx->processData(std::make_pair(buffer, buffer + bytesSent));
		nBytes -= bytesSent;
		buffer += bytesSent;
	}
}

} // namespace daq
} // namespace mdaq
