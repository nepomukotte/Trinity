/* =====================================================================================================================
 * GetEccImpl.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::GetEccImpl
 * Created on: Jan 11, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
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
#ifndef get_GetEccImpl_INCLUDED
#define get_GetEccImpl_INCLUDED

// Ice Generated includes
#include "get/GetEcc.h"
#include "get/cobo/CtrlNode.h"

// Implementation includes
#include "get/rc/BackEnd.h"
#include "mdaq/DaqControl.h"
#include "mdaq/EccImpl.h"
#include "get/rc/BackEnd.h"

#include <Ice/Ice.h>

namespace get
{

class GetEccImpl : public virtual GetEcc, public virtual mdaq::EccImpl
{
protected:
	virtual rc::BackEnd & backEnd();
public:
	GetEccImpl();
	virtual ~GetEccImpl();
	virtual void setConfigManager(const std::string & repoUrl);
	virtual void createAlarmLogger(const ::utl::net::SocketAddress &);
	virtual void daqStart ( const Ice::Current & );
	virtual void daqStop  ( const Ice::Current & );

	/// @name CoBo specific
	/// @{
	virtual void systemReset(const Ice::Current & );
	virtual void coboSetWriteEnable(const Ice::Current & );
	virtual void coboSetReadEnable(const Ice::Current & );
	virtual int  coboReadMultiplicity(const Ice::Current & );
	virtual void coboSetWritingClockFrequency(float frequency, bool pllConfigEnabled, const cobo::PllRegisterMap & registers = cobo::PllRegisterMap(), const Ice::Current & c = Ice::Current());
	virtual void coboSetReadingClockFrequency(float frequency, const Ice::Current & c = Ice::Current());
	virtual void coboSwitchToSlowControlMode(bool slowMode=true, const Ice::Current & c = Ice::Current());
	virtual void coboSelectSingleAgetReadoutChannel(::Ice::Long chanId, const Ice::Current & c = Ice::Current());
	virtual void coboSelectSingleAgetTestChannel(Ice::Long chanId, const Ice::Current & c = Ice::Current());
	virtual void coboSetAgetChannelLsbThreshold(::Ice::Long chanId, ::Ice::Long threshold, const Ice::Current & c = Ice::Current());
	/// @}
	/// @name AsAd voltage generator
	/// @{
	virtual void asadInjectCharge(Ice::Long, const Ice::Current & c = Ice::Current());
	virtual void setAsAdPulserDefaultVoltage(Ice::Long, const Ice::Current & c = Ice::Current());
	virtual void triggerAsAdPulser(Ice::Long, const Ice::Current & c = Ice::Current());
	virtual void configureAsAdExternalPulser(bool enabled, Ice::Long voltage_mV, const Ice::Current & c = Ice::Current());
	virtual void configureAsAdPeriodicPulser(Ice::Long, Ice::Long, const Ice::Current & c = Ice::Current());
	virtual void setRandomAsAdPulserEnabled(bool poissonMode=false, const Ice::Current & c = Ice::Current());
	virtual void configureAsAdDoublePulserMode(Ice::Long mV, const Ice::Current & c = Ice::Current());
	virtual void startAsAdPeriodicPulser(const Ice::Current & c = Ice::Current());
	virtual void stopAsAdPeriodicPulser(const Ice::Current & c = Ice::Current());
	/// @}
	/// @name Alarms
	/// @{
	virtual void setAsAdAlarmMonitoringEnabled(bool enabled, const Ice::Current & c = Ice::Current());
	virtual get::mt::AlarmServicePrx alarmService(const Ice::Current & c = Ice::Current());
	virtual void subscribeAlarms(bool enabled, const Ice::Current & c = Ice::Current());
	/// @}
	/// @name AsAd specific
	/// @{
	virtual void asadSetISPM(const std::string& moduleName, const std::string& signalNameISP1, const std::string& signalNameISP2, const Ice::Current &);
	/// @}
	/// @name AGET specific
	/// @{
	virtual void agetReadHitRegister(std::vector<long>& hitRegisterValue, const Ice::Current & );
	virtual	void agetWriteHitRegister(const ::mdaq::LongList& hitRegisterValue, const Ice::Current & );
	virtual void agetCalibrateSlowControlDelay(const Ice::Current & );
	virtual void agetCalibrateFastControlDelay(const Ice::Current & );
	/// @}

	///@name Transitions of the state machine
	///@{
	virtual void describe(const std::string &, const Ice::Current &);
	virtual void prepare(const std::string &, const Ice::Current &);
	virtual void reset(const Ice::Current &);
	virtual void configure(const std::string &, const std::string &, const Ice::Current &);
	virtual void start(const Ice::Current &);
	virtual void stop(const Ice::Current &);
	virtual void pause(const Ice::Current &);
	virtual void resume(const Ice::Current &);
	virtual void breakup(const Ice::Current &);
	virtual void getStatus(rc::Status &, const Ice::Current &);
	virtual std::string getDataLinks(const Ice::Current &);
	///@}

	///@name Operations on nodes
	///@{
    virtual void describeNode(const std::string & nodeType, const std::string & instanceId, const std::string & descrId, const Ice::Current & c = Ice::Current());
    virtual void prepareNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId, const Ice::Current & c = Ice::Current());
    virtual void configureNode(const std::string & nodeType, const std::string & instanceId, const std::string & configId, const Ice::Current & c = Ice::Current());
	///@}

    /// @name Management of configurations
	///@{ Management of configurations
	virtual std::string addSubConfig(const std::string & subConfigType, const std::string & subConfig, const std::string & subConfigId = "", const Ice::Current & c = Ice::Current());
	virtual std::string getSubConfig(const std::string & subConfigType, const std::string & subConfigId, const Ice::Current & c = Ice::Current());
	virtual std::string getConfigSet(const std::string & xmlConfigId, const Ice::Current & c = Ice::Current());
	virtual std::string getConfigIds(const Ice::Current & c = Ice::Current());
	virtual std::string getSubConfigIds(const std::string & configType, const Ice::Current & c = Ice::Current());
	///@}

	/// @name CoBo PLL device
	///@{
	virtual void loadPLLConfig(const Ice::Float ckwFreq_MHz, const Ice::Current & c = Ice::Current());
	virtual void writePLLRegister(const Ice::Long value, const Ice::Current & c = Ice::Current());
	virtual void resetPLL(const Ice::Current & c = Ice::Current());
	///@}

	/// @name CoBo LED device
	///@{
	virtual get::cobo::LedManagerPrx ledManager(const Ice::Current & c = Ice::Current());
	///@}

	/// @name CoBo-Mutant link
	/// @{
	virtual void setupCoBoMutantLink(get::cobo::MutantLinkMode, const Ice::Current & c = Ice::Current());
	virtual mutant::AlignmentServerPrx mutantAlignmentServer(const Ice::Current & c = Ice::Current());
	/// @}
};
typedef ::IceUtil::Handle < GetEccImpl > GetEccImplPtr;

} // namespace get

#endif /* get_GetEccImpl_INCLUDED */
