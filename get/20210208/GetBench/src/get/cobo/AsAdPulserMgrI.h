/**
 * @file AsAdPulserMgrI.h
 * @date Dec 19, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: AsAdPulserMgrI.h 1541 2016-03-14 12:49:18Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the GetBench software project.
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

#ifndef get_cobo_AsAdPulserMgrI_h_INCLUDED
#define get_cobo_AsAdPulserMgrI_h_INCLUDED

#include "get/cobo/CtrlNode.h"
#include "get/utl/PeriodicTask.h"
#include "mdaq/hw/Control.h"
#include <boost/cstdint.hpp>
#include <boost/system/config.hpp>
#include ICEUTIL_HEADER(Thread.h)

// Forward declarations
namespace mdaq
{
namespace hw
{
namespace server
{
	class BaseDynRegister;
}
}
}

namespace get
{
namespace cobo
{

class CtrlNodeI;
//__________________________________________________________________________________________________
/**
 * Interface for interacting with the AsAd voltage generator.
 */
class AsAdPulserMgrI
{
public:
	AsAdPulserMgrI(CtrlNodeI & node);
	virtual ~AsAdPulserMgrI() {};
	virtual void setDefaultVoltage(Ice::Long defaultVoltage_mV);
	virtual void triggerPulser(Ice::Long);
	virtual void triggerPulserAndRestore(Ice::Long voltage_mV);
	virtual void triggerPulserAndRestore_byDacCode(uint16_t dacCode);
	virtual void resetExternalPulser();
	virtual void configureExternalPulser(bool enabled, Ice::Long voltage_mV);
	virtual void configurePeriodicPulser(Ice::Long, Ice::Long);
	virtual void configureDoublePulserMode(Ice::Long mV=0);
	virtual void setRandomPulserEnabled(bool enabled=false);
	virtual void startPeriodicPulser();
	virtual void startPeriodicPulser(Ice::Long, Ice::Long);
	virtual void stopPeriodicPulser();
private:
	void disableSlowControlTriggerModes(::mdaq::hw::Device & ctrl, bool check=true);
	void restoreSlowControlTriggerModes(::mdaq::hw::Device & ctrl);
	uint16_t dacDecimalCode(const float & step_mV);
protected:
	class PeriodicPulser;
	typedef IceUtil::Handle<PeriodicPulser> PeriodicPulserPtr;
private:
	CtrlNodeI & node_; ///< Reference to CoBo control node for access to CoBo and AsAd registers.
	int16_t voltage_mV_; ///< Amplitude of voltage steps of periodic pulses [mV].
	int16_t voltage2_mV_; ///< Amplitude of second voltage steps of periodic pulses [mV].
	size_t period_ms_; ///< Period [ms] of periodic pulses.
	bool poissonMode_; ///< Whether the pulses follow a Poisson distribution or are periodic
	uint32_t triggerMask_; ///< Cache for value of triggerMask register
	bool slowControlMode_, slowControl2pMode_; ///< Flags to temporarily store trigger mode
	typedef std::map< float, uint16_t > DacCodeMap;
	DacCodeMap dacCodes_;
	uint16_t defaultDacCode_; ///< 14-bit DAC code for generator default output voltage.
	bool externalModeEnabled_;
	uint16_t externalDacCode_; ///< 14-bit DAC code for generator output voltage in external triggering mode.
	PeriodicPulserPtr periodicPulser_;
	static const float Vref_mV;

	struct RegisterCache
	{
		RegisterCache();
		void buildCache(CtrlNodeI &);
		bool isCacheReady_;
		mdaq::hw::server::BaseDynRegister *triggerMask_, *tcm0_, *tcm1_, *tcm2_;
	};

	RegisterCache registers_;
};
//__________________________________________________________________________________________________
class AsAdPulserMgrI::PeriodicPulser : public utl::PeriodicTask,  public IceUtil::Shared
{
public:
	PeriodicPulser(AsAdPulserMgrI & mgr, long voltage_mV = -1000, size_t period_ms = 1000, long voltage2_mV = 0, bool randomEnabled=false);
	~PeriodicPulser();
	virtual void doTask();
private:
	AsAdPulserMgrI & mgr_;
	long voltage_mV_; ///< Amplitude of voltage steps of periodic pulses [mV].
	long voltage2_mV_; ///< Amplitude of second voltage steps of periodic pulses [mV].
};
//__________________________________________________________________________________________________
} /* namespace cobo */
} /* namespace get */
#endif /* get_cobo_AsAdPulserMgrI_h_INCLUDED */
