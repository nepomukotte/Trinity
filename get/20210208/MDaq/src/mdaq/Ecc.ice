/* =====================================================================================================================
 * Ecc.ice
 * ---------------------------------------------------------------------------------------------------------------------
 * mdaq::Ecc remote interface definition.
 * Created on: 4 May 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#ifndef mdaq_Ecc_ice_INCLUDED
#define mdaq_Ecc_ice_INCLUDED

#include "mdaq/utl/CmdException.ice"
#include "mdaq/hw/HardwareTypes.ice"

module mdaq
{	
	sequence<long> LongList;

	interface Ecc
	{
		//__________________________________Management______________________________________//	
		idempotent void testConnectionToEcc();
		idempotent void testConnectionToHardware() throws ::mdaq::utl::CmdException;
		void shutdownEccServer();
		void shutdownHardwareServer() throws ::mdaq::utl::CmdException;
		void rebootHardware() throws ::mdaq::utl::CmdException;
		
		/// General	
		idempotent void loadHwDescription(string cfgFileName) throws ::mdaq::utl::CmdException;
		/** This method is deprecated and remains only for backward compatibility. **/
		void connect(string objectType, string objectName) throws ::mdaq::utl::CmdException;

		/// Node
		idempotent void addNode(string nodeId) throws ::mdaq::utl::CmdException;
		idempotent void removeNode(string nodeId) throws ::mdaq::utl::CmdException;
		idempotent void removeAllNodes() throws ::mdaq::utl::CmdException;
		idempotent void renameNode(string nodeName) throws ::mdaq::utl::CmdException;
		idempotent void selectNodeById(string nodeId) throws ::mdaq::utl::CmdException;
		idempotent void selectNodeByName(string nodeName) throws ::mdaq::utl::CmdException;
		idempotent void getListOfNodes(out hw::NodeList nodes) throws ::mdaq::utl::CmdException;
		idempotent void getEndpointOfCurrentNode(out string endpoint) throws ::mdaq::utl::CmdException;
		idempotent void getNameOfCurrentNode(out string name) throws ::mdaq::utl::CmdException;

		/// Device
		idempotent void addDevice(hw::DeviceDescription deviceName) throws ::mdaq::utl::CmdException;
		idempotent void removeDevice(string deviceName) throws ::mdaq::utl::CmdException;
		idempotent void removeAllDevices() throws ::mdaq::utl::CmdException;
		idempotent void selectCurrentDevice(string deviceName) throws ::mdaq::utl::CmdException;
		idempotent void getListOfDevices(out hw::DeviceList devices) throws ::mdaq::utl::CmdException;
		idempotent void getNameOfCurrentDevice(out string name) throws ::mdaq::utl::CmdException;
		idempotent void getDeviceDescription(out hw::DeviceDescription deviceDescr);
		idempotent string checkDevice(hw::DeviceDescription descr, long regValue, long firstRegOffset, long regCount) throws ::mdaq::utl::CmdException;
		
		/// Register
		idempotent long readRegister(string registerName) throws ::mdaq::utl::CmdException;
		void writeRegister(string registerName, long registerValue) throws ::mdaq::utl::CmdException;
		void addRegister(hw::RegisterConfig cfg) throws ::mdaq::utl::CmdException;
		void addRegisters(hw::RegisterConfigList cfg) throws ::mdaq::utl::CmdException;
		idempotent void getListOfRegisters(out hw::RegisterList registers) throws ::mdaq::utl::CmdException;
		idempotent void getRegisterDescription(string registerName, out hw::RegisterDescription registerDescr);
		
		/// Field
		idempotent long readField(string registerName, string fieldName) throws ::mdaq::utl::CmdException;
		void writeField(string registerName, string fieldName, long fieldValue) throws ::mdaq::utl::CmdException;		
		idempotent void getListOfFields(string registerName, out hw::FieldList fields) throws ::mdaq::utl::CmdException;
		idempotent void getFieldDescription(string registerName, string fieldName, out hw::FieldDescription fieldDescr);

		/// Data Acquisition	
		void daqStart() throws ::mdaq::utl::CmdException;
		void daqStop() throws ::mdaq::utl::CmdException;
		/**
	      * Sets data sender policy.
	      * @param enable Whether to sent all data immediately or to optimize data chunk size.
	      * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
	      **/
	    void setAlwaysFlushData(bool enable) throws ::mdaq::utl::CmdException;
		/**
		 * Selects which of the circular buffers should be read.
		 * @param mask Mask of enabled buffers (1 bit per buffer).
		 * @throws ::mdaq::utl::CmdException Raised no data sender currently exists.
		 */
		void setCircularBuffersEnabled(int asadMask) throws ::mdaq::utl::CmdException;
	};
}; //mdaq

#endif /* Ecc_ice_INCLUDED */
