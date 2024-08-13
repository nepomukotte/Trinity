/**
 * $Id: ControllerOptions.cpp 1741 2017-11-27 09:18:08Z psizun $
 * @file ControllerOptions.cpp
 * @date 6 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors:
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of
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

#include "get/ControllerOptions.h"
#include "mdaq/DefaultPortNums.h"
#include "get/DefaultPortNums.h"
#include "utl/Logging.h"
#include <QSettings>

using ::utl::net::IpAddress;
using ::utl::net::SocketAddress;

namespace get {
//__________________________________________________________________________________________________
ControllerOptions::ControllerOptions()
{
	loadSettings();
}
//__________________________________________________________________________________________________
ControllerOptions::~ControllerOptions()
{
	saveSettings();
}
//__________________________________________________________________________________________________
ControllerOptions & ControllerOptions::instance()
{
	static ControllerOptions options;
	return options;
}
//__________________________________________________________________________________________________
void ControllerOptions::loadSettings()
{
	LOG_DEBUG() << "Loading GetController preferences";
	QSettings settings;
	dataDisplayEnabled_ = settings.value("data/monitoring", true).toBool();
	headerDisplayEnabled_ = settings.value("data/showAsciiHeader", true).toBool();
	saveRawDataEnabled_ = settings.value("data/saveRawFrame", true).toBool();
	dataFlowType_ = settings.value("dataRouter/flowType", "ZBUF").toString().toStdString();
	dataFlowEndpoint_.ipAddress() = utl::net::IpAddress(settings.value("dataRouter/flowAddress", "127.0.0.1").toString().toStdString());
	dataFlowEndpoint_.port() = settings.value("dataRouter/flowPort", ::mdaq::Default::dataRouterFlowPortNum).toInt();
	dataCtrlEndpoint_ = SocketAddress(settings.value("dataRouter/ctrlAddress", "127.0.0.1").toString().toStdString(),
				settings.value("dataRouter/ctrlPort", ::mdaq::Default::dataRouterCtrlPortNum).toInt());
	const SocketAddress defaultAddress(INADDR_ANY, get::Default::alarmLoggerPortNum+1); // Use default address different from ECC server
	alarmLoggerEndpoint_.fromString(settings.value("alarmLoggerAddress", defaultAddress.toString().c_str()).toString().toStdString());
	isDataRouterExternal_ = settings.value("dataRouter/external", false).toBool();
}
//__________________________________________________________________________________________________
void ControllerOptions::saveSettings()
{
	LOG_DEBUG() << "Saving GetController preferences";
	QSettings settings;
	settings.setValue("data/monitoring", QVariant::fromValue(dataDisplayEnabled_));
	settings.setValue("data/showAsciiHeader", QVariant::fromValue(headerDisplayEnabled_));
	settings.setValue("data/saveRawFrame", QVariant::fromValue(saveRawDataEnabled_));
	settings.setValue("dataRouter/flowType", QString::fromStdString(dataFlowType_));
	settings.setValue("dataRouter/flowAddress", QString::fromStdString(dataFlowEndpoint_.ipAddress().toString()));
	settings.setValue("dataRouter/flowPort", dataFlowEndpoint_.port());
	settings.setValue("dataRouter/ctrlAddress", QString::fromStdString(dataCtrlEndpoint_.ipAddress().toString()));
	settings.setValue("dataRouter/ctrlPort", dataCtrlEndpoint_.port());
	settings.setValue("alarmLoggerAddress", QString::fromStdString(alarmLoggerEndpoint_.toString()));
	settings.setValue("dataRouter/external", isDataRouterExternal_);
}
//__________________________________________________________________________________________________
bool ControllerOptions::isDataDisplayEnabled() const
{
	return dataDisplayEnabled_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setDataDisplayEnabled(bool enabled)
{
	dataDisplayEnabled_ = enabled;
}
//__________________________________________________________________________________________________
bool ControllerOptions::isHeaderDisplayEnabled() const
{
	return headerDisplayEnabled_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setHeaderDisplayEnabled(bool enabled)
{
	headerDisplayEnabled_ = enabled;
}
//__________________________________________________________________________________________________
bool ControllerOptions::isSaveRawFrameEnabled() const
{
	return saveRawDataEnabled_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setSaveRawFrameEnabled(bool enabled)
{
	saveRawDataEnabled_ = enabled;
}
//__________________________________________________________________________________________________
::utl::net::SocketAddress ControllerOptions::dataFlowEndpoint() const
{
	return dataFlowEndpoint_;
}
//__________________________________________________________________________________________________
::utl::net::IpAddress ControllerOptions::dataFlowAddress() const
{
	return dataFlowEndpoint_.ipAddress();
}
//__________________________________________________________________________________________________
void ControllerOptions::setDataFlowAddress(const ::utl::net::IpAddress & addr)
{
	dataFlowEndpoint_.ipAddress() = addr;
}
//__________________________________________________________________________________________________
uint16_t ControllerOptions::dataFlowPort() const
{
	return dataFlowEndpoint().port();
}
//__________________________________________________________________________________________________
void ControllerOptions::setDataFlowPort(const uint16_t & port)
{
	dataFlowEndpoint_.port() = port;
}
//__________________________________________________________________________________________________
std::string ControllerOptions::dataFlowType() const
{
	return dataFlowType_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setDataFlowType(const std::string & flowType)
{
	dataFlowType_ = flowType;
}
//__________________________________________________________________________________________________
::utl::net::SocketAddress ControllerOptions::dataCtrlEndpoint() const
{
	return dataCtrlEndpoint_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setDataCtrlEndpoint(const ::utl::net::SocketAddress & endpoint)
{
	dataCtrlEndpoint_ = endpoint;
}
//__________________________________________________________________________________________________
::utl::net::SocketAddress ControllerOptions::alarmLoggerEndpoint() const
{
	return alarmLoggerEndpoint_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setAlarmLoggerEndpoint(const ::utl::net::SocketAddress & endpoint)
{
	alarmLoggerEndpoint_ = endpoint;
}
//__________________________________________________________________________________________________
bool ControllerOptions::isDataRouterExternal() const
{
	return isDataRouterExternal_;
}
//__________________________________________________________________________________________________
void ControllerOptions::setDataRouterExternal(bool enabled)
{
	isDataRouterExternal_ = enabled;
}
//__________________________________________________________________________________________________
} /* namespace get */
