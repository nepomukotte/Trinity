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

#ifndef get_daq_SyncedCircularBuffer_h_INCLUDED
#define get_daq_SyncedCircularBuffer_h_INCLUDED

#include "mdaq/daq/CircularBuffer.h"

namespace get {
namespace daq {

class SyncedCircularBuffer : public mdaq::daq::CircularBuffer
{
public:
	SyncedCircularBuffer(const uint32_t memStartPtrAddr = UINT32_C(0x20000010),
			const uint32_t memEndPtrAddr = UINT32_C(0x20000014),
			const uint32_t memReadPtrAddr = UINT32_C(0x20000020),
			const uint32_t memWritePtrAddr = UINT32_C(0x20000024));
	virtual ~SyncedCircularBuffer();
	void print() const;
	void synchronize(bool write, bool read);
	void synchronizeReadPtr();
	void synchronizeWritePtr();
	uint64_t fpgaEmptySize_B() const;
private:
	uint8_t* & memReadPtr_;  ///< reference to actual FPGA circular buffer READ pointer
	uint8_t* & memWritePtr_; ///< reference to actual FPGA circular buffer WRITE pointer
};

} /* namespace daq */
} /* namespace get */
#endif /* get_daq_SyncedCircularBuffer_h_INCLUDED */
