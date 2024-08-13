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

#ifndef get_daq_MemRead_h_INCLUDED
#define get_daq_MemRead_h_INCLUDED

#include "get/daq/SyncedCircularBuffer.h"
#include "get/CoBoCircularBufferRegisters.h"
#include "mdaq/daq/DataSender.h"
#include <boost/cstdint.hpp>
#include <string>
#include <map>
#include <memory>

namespace get
{

// Forward declaration
namespace cobo
{
	class CtrlNodeI;
}

namespace daq
{
//__________________________________________________________________________________________________
class MemRead
{
public:
	typedef std::auto_ptr< ::mdaq::daq::DataSender > DataSenderPtr;
	/// A class to store statistics about data flow
	struct Stats
	{
		Stats() { reset(); }
		void reset() { byteCount = 0; intrCount = 0; memSize_B = 0; memEmpty_B = 0; }
		uint64_t byteCount; ///< Total number of bytes sent
		uint64_t intrCount; ///< Total number of interrupts received
		uint64_t memEmpty_B; ///< Occupied memory in circular buffer.
		uint64_t memSize_B; ///< Size of memory in circular buffer.
	};
public:
	MemRead(cobo::CtrlNodeI & ctrlNode);
	virtual ~MemRead();
	void resetDataSender(DataSenderPtr sender = DataSenderPtr(0));
	::mdaq::daq::DataSender & dataSender() { return *dataSenderPtr.get(); };
	void setAlwaysFlushData(bool enable=false);
	void sendData();
	void sendData(SyncedCircularBuffer & buffer);
	void flushData();
	void flushData(SyncedCircularBuffer & buffer);
	void updateStatistics();
	void onInterrupt();
	void processInterrupt();
	void sendHeartBeat();
	void createCircularBuffers();
	void createCircularBuffer(const size_t asadId);
	void deleteCircularBuffers();
	void synchronizeCircularBuffers(bool write, bool read);
	bool isBufferEnabled(const size_t asadId) const;
	void setCircularBuffersEnabled(const uint16_t asadMask);
	void resetCircularBuffers();
	void resetCircularBuffer(SyncedCircularBuffer & buffer);
	SyncedCircularBuffer & circularBuffer(const size_t asadId);
	const SyncedCircularBuffer & circularBuffer(const size_t asadId) const;
	const Stats & statistics() const { return statistics_; }
	void sendTopology();
private:
	typedef std::map< size_t, ::get::daq::SyncedCircularBuffer* > CircularBufferPtrMap;
	CircularBufferPtrMap circularBuffers_; ///< Circular buffers indexed by AsAd ID
	uint16_t asadMask_; ///< Mask of AsAd boards whose buffer should be read
	size_t asadCount_;
	DataSenderPtr dataSenderPtr; ///< Data sender
	Stats statistics_; ///< Statistics about data flow
	CoBoCircularBufferRegisters bufferParam; ///< Parameters of the circular buffers
	cobo::CtrlNodeI & ctrlNode_; ///< Reference to CoBo control node
};
//__________________________________________________________________________________________________
} // namespace daq
} // namespace get

#endif //get_daq_MemRead_INCLUDED
