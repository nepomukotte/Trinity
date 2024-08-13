/**
 * @file SyncedCircularBuffer.h
 * @date 2013-01-22
 * @author sizun
 * @author anvar
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Shebli Anvar, Patrick Sizun
 *
 * This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "SyncedCircularBuffer.h"
#include "utl/Logging.h"

namespace get {
namespace daq {

SyncedCircularBuffer::SyncedCircularBuffer(const uint32_t memStartPtrAddr,
			const uint32_t memEndPtrAddr,
			const uint32_t memReadPtrAddr,
			const uint32_t memWritePtrAddr)
	:
		::mdaq::daq::CircularBuffer(memStartPtrAddr, memEndPtrAddr),
		memReadPtr_(* ((uint8_t**) memReadPtrAddr)),
		memWritePtr_(* ((uint8_t**) memWritePtrAddr))
{
}

SyncedCircularBuffer::~SyncedCircularBuffer()
{
}

void SyncedCircularBuffer::print() const
{
	LOG_DEBUG() << "SyncedCircularBuffer: "
				<< std::showbase << std::hex
				<< "memStartPtr(" << (void*) &memStartPtr << ")=" << (void*) memStartPtr << "\t"
				<< "memEndPtr(" << (void*) &memEndPtr << ")=" << (void*) memEndPtr << "\t"
				<< "memReadPtr(" << (void*) &memReadPtr_ << ")=" << (void*) memReadPtr_ << "\t"
				<< "memWritePtr(" << (void*) &memWritePtr_ << ")=" << (void*) memWritePtr_ << std::dec;
				//<< "memStartPtr=" << (void*) memStartPtr << "\tmemEndPtr=" << (void*) memEndPtr
				//<< "\tmemReadPtr_=" << (void*) memReadPtr_ << "\tmemWritePtr_=" << (void*) memWritePtr_ << std::dec;
}

void SyncedCircularBuffer::synchronizeReadPtr()
{
	// Set the read pointer
	memReadPtr_ = readPtr.value();
}

void SyncedCircularBuffer::synchronizeWritePtr()
{
	// Check the write pointer points somewhere within the circular buffer
	if (memWritePtr_ > memEndPtr or memWritePtr_ < memStartPtr)
	{
		LOG_WARN() << "Out of bounds write pointer! memStart=" << (void*) memStartPtr
				<< " writePtr=" << (void*) memWritePtr_ << " memEnd=" << (void*) memEndPtr;
	}
	else
	{
		// Only copy the write pointer if it is valid!
		writePtr.setValue(reinterpret_cast< uint8_t* >(memWritePtr_));
	}
}

void SyncedCircularBuffer::synchronize(bool write, bool read)
{
	if (write) synchronizeWritePtr();
	if (read) synchronizeReadPtr();
}

/*
 * Returns size of fraction of actual FPGA circular buffer which is not filled, in Bytes.
 */
uint64_t SyncedCircularBuffer::fpgaEmptySize_B() const
{
	int64_t delta_B = memReadPtr_ - memWritePtr_;
	if (delta_B <= 0) delta_B += size_B();
	return delta_B;
}

} /* namespace daq */
} /* namespace get */
