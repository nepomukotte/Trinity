/* =====================================================================================================================
 * MemRead.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::daq::MemRead
 * Created on: Apr 4, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * (c) Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GetBench package implementing the GET-specific software modules for the MDaq framework.
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */
#include "MemRead.h"
#include "utl/Logging.h"
#include "utl/BitFieldHelper.hpp"
#include "mdaq/utl/CmdException.h"
#include "get/daq/FakeDataSender.h"
#include "get/cobo/CtrlNodeI.h"
#include "utl/BitFieldHelper.hpp"
#include ICE_HEADER(LocalException.h)

#include <boost/cstdint.hpp>
#include <algorithm> // std::fill
#include <bitset>
#include <cstdio>
#include <cstring> // strerror
#include <utility>
#include <iostream>
#include <iomanip>

#ifdef VXWORKS
#include <vxWorks.h>
#include <logLib.h>
#endif

namespace get
{
namespace daq
{
//_________________________________________________________________________________________________
MemRead::MemRead(cobo::CtrlNodeI & ctrlNode) : ctrlNode_(ctrlNode)
{
	resetDataSender();
	createCircularBuffers();
	// By default, only use buffer for AsAd no. 0 (R-CoBo)
	setCircularBuffersEnabled(0x1);
	resetCircularBuffers();
}
//_________________________________________________________________________________________________
MemRead::~MemRead()
{
	deleteCircularBuffers();
}
//_________________________________________________________________________________________________
/**
 * Sets the data sender for actually sending data.
 * @param sender Data sender auto pointer. The memory reader takes ownership of the sender.
 */
void MemRead::resetDataSender(MemRead::DataSenderPtr sender)
{
	// Set or create new data sender
	if (0 != sender.get())
	{
		dataSenderPtr = sender;
		LOG_INFO() << "Firmware will be issuing an interrupt every " << bufferParam.interruptRate << " event(s)";
	}
	else
	{
		dataSenderPtr.reset(new FakeDataSender());
	}
}
//_________________________________________________________________________________________________
/**
  * Sets data sender flush policy.
  * @param enable Whether to sent all data immediately or to optimize data chunk size.
  * @throws mdaq::utl::CmdException Raised no data sender currently exists.
  **/
void MemRead::setAlwaysFlushData(bool enable)
{
	if  (0 != dataSenderPtr.get())
	{
		dataSenderPtr->setAlwaysFlushData(enable);
		LOG_DEBUG() << "Setting circular buffer 'always flush' policy to " << enable;
	}
	else
	{
		throw ::mdaq::utl::CmdException("Could not set data sender flush policy. Data sender does not exist yet.");
	}
}
//_________________________________________________________________________________________________
/**
 * Creates a circular buffer for all AsAd boards.
 */
void MemRead::createCircularBuffers()
{
	for (size_t asadId=0; asadId < 4u; ++asadId)
	{
		createCircularBuffer(asadId);
	}
}
//_________________________________________________________________________________________________
/**
 * Creates a circular buffer for a specific AsAd board.
 * @param asadId Identifier of the AsAd board (between 0 and 3).
 */
void MemRead::createCircularBuffer(const size_t asadId)
{
	LOG_DEBUG() << "Creating circular buffer for AsAd board no. " << asadId;
	circularBuffers_[asadId] = new SyncedCircularBuffer(
			bufferParam.asad[asadId].memStartPtrAddr,
			bufferParam.asad[asadId].memEndPtrAddr,
			bufferParam.asad[asadId].readPtrAddr,
			bufferParam.asad[asadId].writePtrAddr);
}
//_________________________________________________________________________________________________
/**
 * Deletes all circular buffers.
 */
void MemRead::deleteCircularBuffers()
{
	CircularBufferPtrMap::iterator buffIt;
	for (buffIt = circularBuffers_.begin(); buffIt != circularBuffers_.end(); buffIt++)
	{
		delete buffIt->second;
	}
	circularBuffers_.clear();
}
//_________________________________________________________________________________________________
bool MemRead::isBufferEnabled(const size_t asadId) const
{
	return ::utl::BitFieldHelper< uint16_t >::getBit(asadMask_, asadId);
}
//_________________________________________________________________________________________________
/**
 * Enables readout of circular buffers for a subset of the 4 AsAd boards.
 * @param asadMask Mask of active AsAd circular buffers.
 */
void MemRead::setCircularBuffersEnabled(const uint16_t asadMask)
{
	LOG_INFO() << "Enabling circular buffers with mask " << std::showbase << std::hex << asadMask << std::dec;
	asadMask_ = asadMask;
	asadCount_ = std::bitset<16>(asadMask).count();
}
//_________________________________________________________________________________________________
void MemRead::resetCircularBuffers()
{
	for (size_t asadId=0; asadId < 4u; ++asadId)
	{
		if (isBufferEnabled(asadId))
		{
			resetCircularBuffer(circularBuffer(asadId));
		}
	}
	statistics_.reset();
}
//_________________________________________________________________________________________________
void MemRead::resetCircularBuffer(SyncedCircularBuffer & buffer)
{
	buffer.reset();
	buffer.synchronize(true, true);
	buffer.print();
}
//_________________________________________________________________________________________________
SyncedCircularBuffer & MemRead::circularBuffer(const size_t asadId)
{
	return *circularBuffers_[asadId];
}
//_________________________________________________________________________________________________
const SyncedCircularBuffer & MemRead::circularBuffer(const size_t asadId) const
{
	return *(circularBuffers_.find(asadId)->second);
}
//_________________________________________________________________________________________________
void MemRead::sendData(SyncedCircularBuffer & buffer)
{
	//buffer.synchronize();

	try
	{
		dataSenderPtr->sendData(buffer);
	}
	catch (const mdaq::utl::CmdException & e)
	{
		LOG_ERROR() << "Error sending data: " << e.reason;
	}
	catch (const Ice::SocketException & e)
	{
		LOG_ERROR() << "Error sending data: Ice::SocketException: " <<  std::strerror(e.error);
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Error sending data: " << e.what();
	}
	//buffer.synchronize();
}
//_________________________________________________________________________________________________
/**
 * Sends contents of all enabled circular buffers, following sender flush policy.
 */
void MemRead::sendData()
{
	//LOG_DEBUG() << "MemRead::sendData()";
	flushData();
}
//_________________________________________________________________________________________________
/**
 * Sends all contents of a single circular buffer.
 * @param buffer Circular buffer whose contents to send.
 */
void MemRead::flushData(SyncedCircularBuffer & buffer)
{
	dataSenderPtr->flushData(buffer);
}
//_________________________________________________________________________________________________
/**
 * Sends all contents of all enabled circular buffers.
 */
void MemRead::flushData()
{
	// Get value of write pointers
	synchronizeCircularBuffers(true, false);

	// Update statistics
	updateStatistics();

	try
	{
		// Send data
		CircularBufferPtrMap::iterator buffIt;
		for (buffIt = circularBuffers_.begin(); buffIt != circularBuffers_.end(); buffIt++)
		{
			if (not isBufferEnabled(buffIt->first)) continue;
			flushData(*buffIt->second);
		}
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Error sending data: " << e.what();
	}

	// Set value of read pointers
	synchronizeCircularBuffers(false, true);
}
//_________________________________________________________________________________________________
/**
 * Synchronizes all circular buffers (checks value of write pointer, sets value of read pointer).
 */
void MemRead::synchronizeCircularBuffers(bool write, bool read)
{
	CircularBufferPtrMap::iterator buffIt;
	for (buffIt = circularBuffers_.begin(); buffIt != circularBuffers_.end(); buffIt++)
	{
		if (not isBufferEnabled(buffIt->first)) continue;
		buffIt->second->synchronize(write, read);
	}
}
//_________________________________________________________________________________________________
/**
 * Updates data sender statistics.
 */
void MemRead::updateStatistics()
{
	statistics_.byteCount =  0;
	statistics_.memSize_B = 0;
	statistics_.memEmpty_B = 0;

	for (size_t asadId=0; asadId < 4u; ++asadId)
	{
		if (not isBufferEnabled(asadId)) continue;

		SyncedCircularBuffer & buffer = circularBuffer(asadId);
		statistics_.byteCount +=  buffer.readPtr.dataSum_B();
		statistics_.memSize_B += buffer.size_B();
		statistics_.memEmpty_B += buffer.fpgaEmptySize_B();
	}
}
//_________________________________________________________________________________________________
/**
 * Invokes the data sender's action for timeouts.
 */
void MemRead::sendHeartBeat()
{
	try
	{
		dataSenderPtr->sendHeartBeat();
	}
	catch (const mdaq::utl::CmdException & e)
	{
		LOG_ERROR() << "Data sender failed to send heart beat: " << e.reason;
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Data sender failed to send heart beat: " << e.what();
	}

	// FIXME: just a test to see why small/empty frames are not sent immediately
	try
	{
		flushData();
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Data sender failed to flush data: " << e.what();
	}
}
//_________________________________________________________________________________________________
/**
 * Processes interrupt.
 *  - sends data
 *  - resets AsAd voltage generator's external mode.
 */
void MemRead::processInterrupt()
{
	sendData();
	ctrlNode_.onCoBoInterrupt();
}
//_________________________________________________________________________________________________
/**
 * Handler called after an interrupt, outside any task's context. Must NOT invoke routines that might cause the ISR to block!
 */
void MemRead::onInterrupt()
{
	try
	{
		/*
		if (UTL_LOG_LEVEL <= 0)
			logMsg("Interrupt count: %lu\n", (int) statistics_.intrCount, 0, 0, 0, 0, 0);
		*/
		statistics_.intrCount++;
	}
	catch (const std::exception & e)
	{
#ifdef VXWORKS
		logMsg((char*) e.what(), 0, 0, 0, 0, 0, 0);
#endif
	}
}
//_________________________________________________________________________________________________
/**
 * Sends a special frame to FDT receiver with topology information (active AsAd boards, 2p mode).
 */
void MemRead::sendTopology()
{
	if ("FDT" == dataSender().senderType())
	{
		// Build a frame (revision 0, version 2)
		const size_t frameSize_B = 12;
		uint8_t frameData[frameSize_B];
		// Fill basic header
		std::fill(frameData, frameData+frameSize_B, 0);
		frameData[0] = 0x40; // meta type: big-endian blob of 1 B blocks
		frameData[3] = 0xC; // frame size: 12 blocks of 1 B
		frameData[4] = ctrlNode_.dataSource();
		frameData[6] = 0x7;  // frame type: 7
		// Fill data section
		frameData[8] = ctrlNode_.coboIdx();
		frameData[9] = (uint8_t) asadMask_;
		frameData[10] = ctrlNode_.is2pModeEnabled();
		// Send frame
		LOG_INFO() << "Sending special topology frame of " << frameSize_B << " B";
		dataSender().flushChunk(frameData, frameSize_B);
	}
}
//_________________________________________________________________________________________________
} // namespace daq
} // namespace get
