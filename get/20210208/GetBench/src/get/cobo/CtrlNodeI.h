/**
 * @file CtrlNodeI.h
 * @date Dec 20, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: CtrlNodeI.h 1338 2014-09-02 12:08:58Z psizun $
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

#ifndef get_cobo_CtrlNodeI_h_INCLUDED
#define get_cobo_CtrlNodeI_h_INCLUDED

#include "get/cobo/CtrlNode.h"
#include "get/cobo/AsAdPulserMgrI.h"
#include "get/cobo/AsAdAlarmWatcher.h"
#include "get/cobo/LedManagerI.h"
#include "get/mt/AlarmService.h"
#include "mdaq/hw/server/HwNodeImpl.h"

namespace get
{
namespace cobo
{
//__________________________________________________________________________________________________
/**
 * Implementation of the get::cobo::CtrlNode ICE interface.
 */
class CtrlNodeI : public virtual CtrlNode, public virtual ::mdaq::hw::server::HwNodeImpl
{
public:
	CtrlNodeI();
	virtual ~CtrlNodeI();
	virtual void onCoBoInterrupt();
	bool is2pModeEnabled();
	uint8_t coboIdx();
	uint8_t dataSource();
	/// @name AsAd Voltage Generator
	/// @{
	virtual void setDefaultPulserVoltage(Ice::Long, const Ice::Current & c = ::Ice::Current());
	virtual void triggerPulser(::Ice::Long, const ::Ice::Current & c = ::Ice::Current());
	virtual void configureExternalPulser(bool enabled, Ice::Long voltage=0, const Ice::Current & c = ::Ice::Current());
	virtual void configurePeriodicPulser(::Ice::Long, ::Ice::Long, const ::Ice::Current & c = ::Ice::Current());
	virtual void setRandomPulserEnabled(bool enabled=false, const Ice::Current & c = ::Ice::Current());
	virtual void configureDoublePulserMode(Ice::Long voltage_mV=0, const Ice::Current & c = ::Ice::Current());
	virtual void startPeriodicPulser(const ::Ice::Current & c = ::Ice::Current());
	virtual void stopPeriodicPulser(const ::Ice::Current & c = ::Ice::Current());
	/// @}
	/// @name Monitoring of AsAd alarms
	/// @{
	virtual void setAsAdAlarmMonitoringEnabled(bool, const Ice::Current & c = ::Ice::Current());
	/// @}
	/// @name Management of CoBo LEDs
	/// @{
	LedManagerI & ledManager();
	virtual void setLED(int ledId, bool enabled=true, const Ice::Current & c = ::Ice::Current());
	virtual void setLEDs(bool enabled=true, const Ice::Current & c = ::Ice::Current());
	virtual void blinkLEDs(const Ice::Current & c = ::Ice::Current());
	virtual void flipLED(int ledId, const Ice::Current & c = ::Ice::Current());
	virtual void flipLEDs(const Ice::Current & c = ::Ice::Current());
	virtual void pulseLED(int ledId, int period_ms, const Ice::Current & c = ::Ice::Current());
	virtual void modifyLED(LedType, int, LedState, const ::Ice::Current & c = ::Ice::Current());
	/// @}
	mt::AlarmServicePrx getAlarmService(const Ice::Current & c = ::Ice::Current());
private:
	AsAdPulserMgrI asadPulserMgr;
	LedManagerI ledMgr_;
	AsAdAlarmWatcherPtr asadAlarmWatcher_;
};
//__________________________________________________________________________________________________
} /* namespace cobo */
} /* namespace get */
#endif /* get_cobo_CtrlNodeI_h_INCLUDED */
