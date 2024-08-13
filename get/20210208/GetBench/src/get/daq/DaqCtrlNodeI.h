/* =====================================================================================================================
 * DaqCtrlNodeI.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::daq::DaqCtrlNodeI
 * Created on: Jun 21, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef get_daq_DaqCtrlNodeI_h_INCLUDED
#define get_daq_DaqCtrlNodeI_h_INCLUDED

#include "mdaq/mt/PollingMonitorAdapter.hpp"
#include "mdaq/mt/InterruptMonitorAdapter.hpp"
#include "mdaq/DaqControl.h"
#include "get/daq/SyncedCircularBuffer.h"
#include "MemRead.h"
#include ICE_HEADER(Time.h)
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
class DaqCtrlNodeI : public mdaq::DaqCtrlNode
{
public:
	typedef std::auto_ptr < ::get::daq::MemRead > MemReadPtr;
	typedef std::auto_ptr < ::mdaq::mt::InterruptMonitorAdapter < MemRead, MemRead > > InterruptMonitorAdapterPtr;
public:
	DaqCtrlNodeI(cobo::CtrlNodeI & ctrlNode, int timerPeriod_s=10);
	virtual ~DaqCtrlNodeI();
    virtual void connect(const std::string& dataRouterType, const std::string& dataRouterRef, const Ice::Current & c = Ice::Current());
    virtual void disconnect(const Ice::Current & c = Ice::Current());
    virtual void setAlwaysFlushData(bool enable=false, const Ice::Current & c = Ice::Current());
    virtual void setCircularBuffersEnabled(int asadMask, const Ice::Current & c = Ice::Current());
    virtual void daqStart(const Ice::Current & c = Ice::Current());
    virtual void daqStop(const Ice::Current & c = Ice::Current());
    virtual void shutdown(const Ice::Current & c = Ice::Current());
    virtual void flushData(const Ice::Current & c = Ice::Current());
    static void dispatch(int imParam, int timeout);
	virtual void printStatistics(int period_ticks);
    static void statistics_printer(int obj, int arg);
private:
    void createDataSender(const std::string & dataRouterType, const std::string & dataRouterRef, const Ice::Current & c = Ice::Current());
    void createStatisticsTask();
    void deleteStatisticsTask();
    MemRead memReader;
    InterruptMonitorAdapterPtr interruptMon;
	Ice::CommunicatorPtr ic;
	int timerPeriod_s;
	int statsTaskId;
	IceUtil::Time previousTime; ///< Date and time of previous print of statistics
	IceUtil::Time runStartTime; ///< Date and time of start of run
	MemRead::Stats previousStats;
	static const char* rateUnits[];
	bool isRunning_;
};

} // namespace daq
} // namespace get

#endif //get_daq_DaqCtrlNodeI_h_INCLUDED
