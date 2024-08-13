/* =====================================================================================================================
 * EccImpl.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::ecc::server::EccImpl
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
#include "EccImpl.h"
#include "mdaq/cfg/NodeAtTop.h"
#include "mdaq/DefaultPortNums.h"
#include "mdaq/utl/CmdException.h"
#include "utl/Logging.h"

#include ICE_HEADER(ObjectAdapter.h)

#include <exception>

namespace mdaq
{
//_________________________________________________________________________________________________
EccImpl::EccImpl() : backEndPtr(0)
{
}
//_________________________________________________________________________________________________
EccImpl::~EccImpl()
{
}
//_________________________________________________________________________________________________
EccBackEnd & EccImpl::backEnd()
{
	if (backEndPtr.get() == 0)
	{
		backEndPtr.reset(new EccBackEnd());
	}
	return *backEndPtr.get();
}
//_________________________________________________________________________________________________
void EccImpl::testConnectionToEcc(const Ice::Current&)
{
	backEnd().testConnectionToEcc();
}
//_________________________________________________________________________________________________
void EccImpl::testConnectionToHardware(const Ice::Current&)
{
	backEnd().testConnectionToHardware();
}
//_________________________________________________________________________________________________
void EccImpl::shutdownEccServer(const Ice::Current& c)
{
	LOG_INFO() << "Shutting down ECC server...";
	c.adapter->getCommunicator()->shutdown();
}
//_________________________________________________________________________________________________
void EccImpl::shutdownHardwareServer(const Ice::Current&)
{
	backEnd().shutdownHardwareServer();
}
//_________________________________________________________________________________________________
void EccImpl::rebootHardware(const Ice::Current&)
{
	backEnd().rebootHardware();
}
//_________________________________________________________________________________________________
void EccImpl::loadHwDescription(const std::string& cfgFileName, const Ice::Current&)
{
	backEnd().loadHwDescription(cfgFileName);
}
//_________________________________________________________________________________________________
void EccImpl::writeField(const std::string& registerName, const std::string& fieldName, Ice::Long fieldValue, const Ice::Current&)
{
	backEnd().writeField(registerName, fieldName, fieldValue);
}
//_________________________________________________________________________________________________
Ice::Long EccImpl::readField(const std::string& registerName, const std::string& fieldName, const Ice::Current&)
{
	return backEnd().readField(registerName, fieldName);
}
//_________________________________________________________________________________________________
void EccImpl::addRegister(const ::mdaq::hw::RegisterConfig & cfg, const ::Ice::Current &)
{
	backEnd().addRegister(cfg);
}
//_________________________________________________________________________________________________
void EccImpl::addRegisters(const ::mdaq::hw::RegisterConfigList & listOfRegisters, const ::Ice::Current &)
{
	backEnd().addRegisters(listOfRegisters);
}
//_________________________________________________________________________________________________
void EccImpl::writeRegister(const std::string& registerName, Ice::Long regValue, const Ice::Current&)
{
	backEnd().writeRegister(registerName, regValue);
}
//_________________________________________________________________________________________________
Ice::Long EccImpl::readRegister(const std::string& registerName, const Ice::Current&)
{
	return backEnd().readRegister(registerName);
}
//_________________________________________________________________________________________________
/*
  * Sets data sender policy.
  * @param enable Whether to sent all data immediately or to optimize data chunk size.
  * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
  */
void EccImpl::setAlwaysFlushData(bool enable, const Ice::Current &)
{
	backEnd().setAlwaysFlushData(enable);
}
//_________________________________________________________________________________________________
/*
  * Selects which of the circular buffers should be read.
  * @param mask Mask of enabled buffers (1 bit per buffer).
  * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
  */
void EccImpl::setCircularBuffersEnabled(int mask, const Ice::Current &)
{
	backEnd().setCircularBuffersEnabled(mask);
}
//_________________________________________________________________________________________________
void EccImpl::daqStart(const Ice::Current&)
{
	backEnd().daqStart();
}
//_________________________________________________________________________________________________
void EccImpl::daqStop(const Ice::Current&)
{
	backEnd().daqStop();
}
//_________________________________________________________________________________________________
void EccImpl::connect(const std::string& objectType, const std::string& objectName, const Ice::Current&)
{
	backEnd().connect(objectType, objectName);
}
//_________________________________________________________________________________________________
void EccImpl::addNode(const std::string & nodeAddress, const Ice::Current&)
{
	backEnd().addHwNode(::utl::net::SocketAddress(nodeAddress));
}
//_________________________________________________________________________________________________
void EccImpl::renameNode(const std::string & nodeName, const Ice::Current&)
{
	backEnd().renameHwNode(nodeName);
}
//_________________________________________________________________________________________________
void EccImpl::removeNode(const std::string & nodeId, const Ice::Current&)
{
	backEnd().removeHwNode(nodeId);
}
//_________________________________________________________________________________________________
void EccImpl::removeAllNodes(const Ice::Current&)
{
	backEnd().removeAllHwNodes();
}
//_________________________________________________________________________________________________
void EccImpl::selectNodeById(const std::string & nodeId, const Ice::Current&)
{
	backEnd().selectHwNodeByEndpoint(nodeId);
}
//_________________________________________________________________________________________________
void EccImpl::selectNodeByName(const std::string & nodeName, const Ice::Current&)
{
	backEnd().selectHwNodeByName(nodeName);
}
//_________________________________________________________________________________________________
void EccImpl::addDevice(const mdaq::hw::DeviceDescription & descr, const Ice::Current &)
{
	backEnd().addDevice(descr);
}
//_________________________________________________________________________________________________
std::string EccImpl::checkDevice(const ::mdaq::hw::DeviceDescription & descr, Ice::Long regValue, Ice::Long firstRegOffset, Ice::Long regCount, const Ice::Current &)
{
	return backEnd().checkDevice(descr, regValue, firstRegOffset, regCount);
}
//_________________________________________________________________________________________________
void EccImpl::removeDevice(const std::string & deviceName, const Ice::Current&)
{
	backEnd().removeDevice(deviceName);
}
//_________________________________________________________________________________________________
void EccImpl::removeAllDevices(const Ice::Current&)
{
	backEnd().removeAllDevices();
}
//_________________________________________________________________________________________________
void EccImpl::selectCurrentDevice(const std::string & deviceName, const Ice::Current&)
{
	backEnd().selectCurrentDevice(deviceName);
}
//_________________________________________________________________________________________________
void EccImpl::getListOfNodes(hw::NodeList & nodeList, const Ice::Current&)
{
	backEnd().getListOfHwNodes(nodeList);
}
//_________________________________________________________________________________________________
void EccImpl::getListOfDevices(hw::DeviceList & devices, const Ice::Current&)
{
	backEnd().getListOfDevices(devices);
}
//_________________________________________________________________________________________________
void EccImpl::getListOfRegisters(hw::RegisterList & registers, const Ice::Current&)
{
	backEnd().getListOfRegisters(registers);
}
//_________________________________________________________________________________________________
void EccImpl::getListOfFields(const std::string & registerName, hw::FieldList & fields, const Ice::Current&)
{
	backEnd().getListOfFields(registerName, fields);
}
//_________________________________________________________________________________________________
void EccImpl::getNameOfCurrentNode(std::string & nodeName, const Ice::Current &)
{
	backEnd().getNameOfCurrentHwNode(nodeName);
}
//_________________________________________________________________________________________________
void EccImpl::getNameOfCurrentDevice(std::string & deviceName, const Ice::Current &)
{
	backEnd().getNameOfCurrentDevice(deviceName);
}
//_________________________________________________________________________________________________
void EccImpl::getEndpointOfCurrentNode(std::string & address, const Ice::Current &)
{
	backEnd().getEndpointOfCurrentHwNode(address);
}
//_________________________________________________________________________________________________
void EccImpl::getDeviceDescription(hw::DeviceDescription & deviceDescription, const Ice::Current&)
{
	backEnd().getDeviceDescription(deviceDescription);
}
//_________________________________________________________________________________________________
void EccImpl::getRegisterDescription(const std::string & registerName, hw::RegisterDescription & descr, const Ice::Current&)
{
	backEnd().getRegisterDescription(registerName, descr);
}
//_________________________________________________________________________________________________
void EccImpl::getFieldDescription(const std::string & registerName, const std::string & fieldName, hw::FieldDescription & descr, const Ice::Current&)
{
	backEnd().getFieldDescription(registerName, fieldName, descr);
}
//_________________________________________________________________________________________________
} //mdaq
