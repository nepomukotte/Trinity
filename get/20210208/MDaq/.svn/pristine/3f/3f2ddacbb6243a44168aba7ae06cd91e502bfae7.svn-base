/**
 * @file CircularBuffer.h
 * @date 2013-01-22
 * @author sizun
 * @author anvar
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Shebli Anvar, Patrick Sizun
 *
 * This file is part of the MDaq software project.
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

#ifndef mdaq_daq_CircularBuffer_h_INCLUDED
#define mdaq_daq_CircularBuffer_h_INCLUDED

#include <boost/cstdint.hpp>

namespace mdaq {
namespace daq {

class CircularBuffer
{
public:
	struct Pointer
	{
		Pointer(CircularBuffer & circularBuffer) : circularBuffer(circularBuffer) { reset(); }
		void incr(const size_t nBytes);
		void reset() { ptrValue = circularBuffer.memStartPtr; dataSum_B_ = 0u; }
		const uint8_t* value() const { return ptrValue; }
		uint8_t* value() { return ptrValue; }
		uint64_t dataSum_B() { return dataSum_B_; }
		void setValue(uint8_t* value);
	private:
		uint8_t* ptrValue;
		uint64_t dataSum_B_;
		CircularBuffer & circularBuffer;
	};
public:
	CircularBuffer(const uint64_t memStartPtrAddr, const uint64_t memEndPtrAddr);
	virtual ~CircularBuffer();
	virtual void reset();
	uint64_t nextChunkSize_B(bool & bufferWrap) const;
	uint64_t nextChunkSize_B() const;
	uint64_t size_B() const { return memEndPtr - memStartPtr; }
	uint64_t occupiedSize_B();
protected:
	uint8_t* & memStartPtr; ///< Reference to actual FPGA circular buffer START pointer
	uint8_t* & memEndPtr; ///< Reference to actual FPGA circular buffer END pointer
public:
	Pointer readPtr; ///< Software copy of circular buffer READ pointer
	Pointer writePtr; ///< Software copy of circular buffer WRITE pointer
};

inline void CircularBuffer::Pointer::incr(const size_t nBytes)
{
	setValue(ptrValue + nBytes);
}

inline void CircularBuffer::Pointer::setValue(uint8_t* value)
{
	if (value >= ptrValue)
		dataSum_B_ += value - ptrValue;
	else
		dataSum_B_ += (circularBuffer.memEndPtr - ptrValue) + (value - circularBuffer.memStartPtr);

	ptrValue = value;
	while (ptrValue >= circularBuffer.memEndPtr)
	{
		ptrValue -= circularBuffer.size_B();
	}
}

inline uint64_t CircularBuffer::nextChunkSize_B(bool & bufferWrap) const
{
	int64_t occupiedSize_B = writePtr.value() - readPtr.value();
	if (occupiedSize_B < 0)
	{
		bufferWrap = true;
		return memEndPtr - readPtr.value();
	}
	else
	{
		bufferWrap = false;
		return (uint64_t) occupiedSize_B;
	}
}

inline uint64_t CircularBuffer::nextChunkSize_B() const
{
	bool bufferWrap;
	return nextChunkSize_B(bufferWrap);
}

} /* namespace daq */
} /* namespace mdaq */
#endif /* mdaq_daq_CircularBuffer_h_INCLUDED */
