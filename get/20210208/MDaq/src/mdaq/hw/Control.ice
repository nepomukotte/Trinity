/* =====================================================================================================================
 * Control.ice
 * ---------------------------------------------------------------------------------------------------------------------
 * mdaq::hw::Device  remote interface definition.
 * mdaq::hw::Node    remote interface definition.
 * Created on: 4 May 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of
 * the Mordicus Real-Time Software Development Framework.
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

#ifndef mdaq_hw_Control_ice_INCLUDED
#define mdaq_hw_Control_ice_INCLUDED

#include "mdaq/utl/CmdException.ice"
#include "mdaq/hw/HardwareTypes.ice"

module mdaq
{
	module hw
	{	
		/** Forward declaration of the Device interface. **/
		interface Device;
		
		/** A key-value map of Device proxies indexed by their name. **/ 
		dictionary<string, Device*> DeviceMap;
		
		/** An interface for communication with an hardware node. **/
		interface Node
		{
			//__________________________________Management______________________________________//	
			idempotent void testConnectionToHardware();
			void shutdown();
			void reboot();	
			//__________________________________Hardware Node____________________________________//
			idempotent string name();
			idempotent void setName(string name);
		
			void create(NodeConfig config) throws ::mdaq::utl::CmdException;
			idempotent void destroy() throws ::mdaq::utl::CmdException;
			
			idempotent void getNodeConfig(out NodeConfig config);
			idempotent void getListOfDevices(out DeviceList devices);
			idempotent void getMapOfDevices(out DeviceMap devices);
			Device* findDevice(string name) throws ::mdaq::utl::CmdException;
			idempotent void removeDevice(string deviceName) throws ::mdaq::utl::CmdException;
			idempotent void removeAllDevices() throws ::mdaq::utl::CmdException;
			void addDevice(DeviceDescription config) throws ::mdaq::utl::CmdException;
		};

		// Register Command Types
		enum RegInstruction
		{
			Nop, // No operation is performed
			Sleep, // Sleep operation with parameter in milliseconds
			Read, // Register read operation
			Write // Register write operation
		};

		// Register Command Structure		
		struct RegCmd
		{
			RegInstruction instruction;
			string regName;
			string fieldName; // bit field name; if empty, the whole register is considered 
			long value; // value to be written in write commands, delay value in sleep command, return value for reads...
		};
		
		sequence<RegCmd> RegCmdSeq;
		
		interface Device
		{
			idempotent string name();
			idempotent long	baseAddress();
			idempotent string registerAccess();
			idempotent short registerWidth();
			
			void addRegister(RegisterConfig cfg) throws ::mdaq::utl::CmdException;
			void addRegisters(RegisterConfigList cfg) throws ::mdaq::utl::CmdException;
		
			idempotent long readRegister(string registerName) throws ::mdaq::utl::CmdException;
			idempotent long readField(string registerName, string fieldName) throws ::mdaq::utl::CmdException;
		
			void writeRegister(string registerName, long value) throws ::mdaq::utl::CmdException;
			void writeField(string registerName, string fieldName, long value) throws ::mdaq::utl::CmdException;
		
			void execBatch(RegCmdSeq regCmdBatch, out RegCmdSeq regCmdBatchStatus) throws ::mdaq::utl::CmdException;
		
			idempotent void getDeviceDescription(out DeviceDescription descr) throws ::mdaq::utl::CmdException;
			idempotent void getDeviceConfig(out DeviceConfig config);
		
			idempotent void getListOfRegisters(out RegisterList registers);
			idempotent void getListOfFields(string registerName, out FieldList fields) throws ::mdaq::utl::CmdException;
			
			idempotent void getRegisterDescription(string registerName, out RegisterDescription descr) throws ::mdaq::utl::CmdException;
			idempotent void getRegisterConfig(string registerName, out RegisterConfig config) throws ::mdaq::utl::CmdException;	
			
			idempotent void getFieldDescription(string registerName, string fieldName, out FieldDescription descr) throws ::mdaq::utl::CmdException;
			idempotent void setVerbose(bool verbosity);
			idempotent bool isVerbose();
		};
	};
};

#endif /* mdaq_hw_Control_ice_INCLUDED */
