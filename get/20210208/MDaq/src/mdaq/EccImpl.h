/* =====================================================================================================================
 * EccImpl.h
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
#ifndef mdaq_EccImpl_INCLUDED
#define mdaq_EccImpl_INCLUDED

#include "mdaq/Ecc.h"
#include "mdaq/EccBackEnd.h"
#include <memory>

namespace mdaq
{

class EccImpl : public virtual Ecc
{
protected:
	virtual EccBackEnd & backEnd();
public:
	EccImpl();
	virtual ~EccImpl();

	/// @name Management
	/// @{
	virtual void testConnectionToEcc(const Ice::Current &);
	virtual void testConnectionToHardware(const Ice::Current &);
	virtual void shutdownEccServer(const Ice::Current &);
	virtual void shutdownHardwareServer(const Ice::Current &);
	virtual void rebootHardware(const Ice::Current &);
	///@}

	/// @name General
	///@{
	virtual void loadHwDescription(const std::string& cfgFileName, const Ice::Current & c = Ice::Current());
	virtual void daqStart(const Ice::Current & c = Ice::Current());
	virtual void daqStop(const Ice::Current & c = Ice::Current());
	virtual void connect(const std::string& objectType, const std::string& objectName, const Ice::Current & c = Ice::Current());
	virtual void setAlwaysFlushData(bool, const Ice::Current & c = Ice::Current());
	virtual void setCircularBuffersEnabled(int, const Ice::Current & c = Ice::Current());
	///@}

	/// @name Node
	/// @{
	virtual void addNode(const std::string& nodeAddress, const Ice::Current & c = Ice::Current());
	virtual void renameNode(const std::string & nodeName, const Ice::Current& c = Ice::Current());
	virtual void removeNode(const std::string& nodeName, const Ice::Current & c = Ice::Current());
	virtual void removeAllNodes(const Ice::Current & c = Ice::Current());
	virtual void selectNodeById(const std::string& nodeId, const Ice::Current & c = Ice::Current());
	virtual void selectNodeByName(const std::string& nodeName, const Ice::Current & c = Ice::Current());
	virtual void getListOfNodes(hw::NodeList& nodes, const Ice::Current & c = Ice::Current());
	virtual void getNameOfCurrentNode(std::string & nodeId, const Ice::Current & c = Ice::Current());
	virtual void getEndpointOfCurrentNode(std::string & address, const Ice::Current & c = Ice::Current());
	/// @}

	/// @name Device
	/// @{
	virtual void addDevice(const mdaq::hw::DeviceDescription&, const Ice::Current & c = Ice::Current());
	virtual void removeDevice(const std::string& deviceName, const Ice::Current & c = Ice::Current());
	virtual void removeAllDevices(const Ice::Current & c = Ice::Current());
	virtual void selectCurrentDevice(const std::string& deviceName, const Ice::Current & c = Ice::Current());
	virtual void getListOfDevices(hw::DeviceList& devices, const Ice::Current & c = Ice::Current());
	virtual void getNameOfCurrentDevice(std::string & deviceName, const Ice::Current & c = Ice::Current());
	virtual void getDeviceDescription(hw::DeviceDescription& deviceDescription, const Ice::Current & c = Ice::Current());
	virtual std::string checkDevice(const ::mdaq::hw::DeviceDescription & descr, Ice::Long, Ice::Long, Ice::Long, const Ice::Current & c = Ice::Current());
	/// @}

	/// @name Register
	/// @{
	virtual Ice::Long readRegister(const std::string& registerName, const Ice::Current & c = Ice::Current());
	virtual void writeRegister(const std::string& registerName, Ice::Long fieldValue, const Ice::Current & c = Ice::Current());
	virtual void addRegister(const ::mdaq::hw::RegisterConfig & cfg, const ::Ice::Current & c = Ice::Current());
	virtual void addRegisters(const ::mdaq::hw::RegisterConfigList & cfg, const ::Ice::Current & c = Ice::Current());
	virtual void getListOfRegisters(hw::RegisterList& registers, const Ice::Current & c = Ice::Current());
	virtual void getRegisterDescription(const std::string& registerName, hw::RegisterDescription & descr, const Ice::Current & c = Ice::Current());
	/// @}

	/// @name Field
	/// @{
	virtual Ice::Long readField(const std::string& registerName, const std::string& fieldName, const Ice::Current & c = Ice::Current());
	virtual void writeField(const std::string& registerName, const std::string& fieldName, Ice::Long fieldValue, const Ice::Current & c = Ice::Current());
	virtual void getListOfFields(const std::string& registerName, hw::FieldList& fields, const Ice::Current & c = Ice::Current());
	virtual void getFieldDescription(const std::string& registerName, const std::string& fieldName, hw::FieldDescription & descr, const Ice::Current & c = Ice::Current());
	/// @}
protected:
	std::auto_ptr< EccBackEnd > backEndPtr;
};

} // namespace mdaq

#endif /* mdaq_EccImpl_INCLUDED */
