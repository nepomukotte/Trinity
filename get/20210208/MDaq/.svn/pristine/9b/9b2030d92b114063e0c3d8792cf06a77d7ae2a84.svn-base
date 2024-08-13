/* =====================================================================================================================
 * HwNodeImpl.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::hw::server::HwNodeImpl
 * Created on: 19 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
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

#ifndef mdaq_hw_server_HwNodeImpl_h_INCLUDED
#define mdaq_hw_server_HwNodeImpl_h_INCLUDED

#include "mdaq/hw/Control.h"
#include ICE_HEADER(ICE_EDITION.h)
#include <vector>
#include <string>

namespace mdaq
{
namespace hw
{
namespace server
{
class DeviceImpl;
/**
 * This class implements the ICE interface allowing communication with a hardware
 * node containing devices.
 */
class HwNodeImpl : public virtual ::mdaq::hw::Node
{
public:
	HwNodeImpl();
	virtual ~HwNodeImpl();

	void testConnectionToHardware(const Ice::Current & c = Ice::Current());
	void shutdown(const Ice::Current & c = Ice::Current());
	void reboot(const Ice::Current & c = Ice::Current());

	std::string name(const Ice::Current & c = Ice::Current());
	void setName(const std::string & name, const ::Ice::Current & c = Ice::Current());

	void create(const NodeConfig& config, const Ice::Current & c = Ice::Current());
	void destroy(const Ice::Current & c = Ice::Current());
	void getNodeDescription(std::string & nodeID, const Ice::Current & c = Ice::Current());
	void getNodeConfig(NodeConfig& nodeConfig, const Ice::Current & c = Ice::Current());
	/// @name Devices
	/// @{
	virtual void addDevice(const mdaq::hw::DeviceDescription&, const Ice::Current & c = Ice::Current());
	void getMapOfDevices(DeviceMap & deviceMap, const ::Ice::Current & c = Ice::Current());
	void getListOfDevices(DeviceList& devices, const Ice::Current & c = Ice::Current());
	DevicePrx findDevice(const std::string & name, const Ice::Current & c = Ice::Current());
	DeviceImpl & findLocalDevice(const std::string & name);
	virtual void removeDevice(const std::string & deviceName, const Ice::Current & c = Ice::Current());
	void removeAllDevices(const Ice::Current & c = Ice::Current());
	/// @}
private:
	void createDevice(const DeviceConfig & config, const Ice::Current & c = Ice::Current());
	Ice::Identity makeDeviceIdentity(const std::string& deviceName);

	std::string nodeName; ///< Name of the node.
	DeviceMap deviceMap; ///< A map of Device proxies in the node with the device name as key.
	Ice::ObjectAdapterPtr adapter_; ///< Adapter for local access to device objects
};

}
}
}

#endif // mdaq_hw_server_HwNodeImpl_h_INCLUDED
