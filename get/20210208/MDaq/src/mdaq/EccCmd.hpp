/* =====================================================================================================================
 * EccCmd.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::EccCmd
 * class mdaq::EccCmd::Engine
 * Created on: Jan 5, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq data acquisition framework @ CEA-Irfu.
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

#ifndef mdaq_EccCmd_hpp_INCLUDED
#define mdaq_EccCmd_hpp_INCLUDED

#include "mdaq/EccCmd.h"
#include "mdaq/utl/CmdException.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/Logging.h"

namespace mdaq
{

/**
 * This command loads the register configuration into the ECC.
 */
template <class Proxy = EccPrx>
struct Load : public EccCmd<Proxy>
{
	Load(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "load", "<cconfig file>", "Loads a register description configuration file into the current node") {}
	void exec(std::istream& cmdArgs)
	{
		std::string cfgFileName;
		cmdArgs >> cfgFileName;
		EccCmd<Proxy>::ecc()->loadHwDescription(cfgFileName);
	}
};

/**
 * This command prints the list of files in a given register of the current device of the current node.
 */
template <class Proxy = EccPrx>
struct ListFields : public EccCmd<Proxy>
{
	ListFields(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "field-list", "<register-name>", "Lists the bit fields in a given register of the current device") {}
	void exec(std::istream & cmdArgs)
	{
		std::string registerName;
		cmdArgs >> registerName;

		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		hw::FieldList fields;
		EccCmd<Proxy>::ecc()->getListOfFields(registerName, fields);
		std::cout << "Bit fields in register '" << registerName << "':";
		for (hw::RegisterList::const_iterator nameIt = fields.begin(); nameIt != fields.end(); ++nameIt)
		{
			std::cout << " '" << *nameIt << '\'';
		}
		std::cout << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command writes a value on a register bit field.
 */
template <class Proxy = EccPrx>
struct WriteField : public EccCmd<Proxy>
{
	WriteField(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "field-write", "<register name> <field name> <value>", "Writes to a bitfield of the current node") {}
	void exec(std::istream& cmdArgs)
	{
		std::string registerName, fieldName;
		uint64_t value;
		cmdArgs >> registerName;
		cmdArgs >> fieldName;
		cmdArgs.unsetf(std::ios::dec);
		cmdArgs >> value;
		EccCmd<Proxy>::ecc()->writeField(registerName, fieldName, value);
	}
};

/**
 * This command reads a value from a register bit field.
 */
template <class Proxy = EccPrx>
struct ReadField : public EccCmd<Proxy>
{
	ReadField(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "field-read", "<register name> <field name>", "Reads from a bitfield of the current node") {}
	void exec(std::istream& cmdArgs)
	{
		std::string registerName, fieldName;
		uint64_t value;
		cmdArgs >> registerName;
		cmdArgs >> fieldName;
		value = EccCmd<Proxy>::ecc()->readField(registerName, fieldName);
		std::cout << value << " : 0x" << std::hex << value << std::dec << std::endl;
	}
};

/**
 * This command adds a register to the current device.
 */
template <class Proxy = EccPrx>
struct AddRegister : public EccCmd<Proxy>
{
	AddRegister(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "reg-add", "<name> <offset>", "Adds a register to the current device") {}
	void exec(std::istream & cmdArgs)
	{
		mdaq::hw::RegisterConfig cfg;
		cfg.descr.offset=0;
		cfg.fields.clear();
		cmdArgs >> cfg.descr.name;
		cmdArgs.unsetf(std::ios::dec);
		cmdArgs >> cfg.descr.offset;

		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		LOG_INFO() << "Adding register '" << cfg.descr.name << "' with offset 0x" << std::hex << cfg.descr.offset;

		EccCmd<Proxy>::ecc()->addRegister(cfg);
	}
};

/**
 * This command writes a value on a whole register.
 */
template <class Proxy = EccPrx>
struct WriteReg : public EccCmd<Proxy>
{

	WriteReg(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "reg-write", "<register name> <value>", "Writes to a register of the current node") {}
	void exec(std::istream& cmdArgs)
	{
		std::string registerName;
		uint64_t value;
		cmdArgs >> registerName;
		cmdArgs.unsetf(std::ios::dec);
		cmdArgs >> value;
		if (cmdArgs.fail())
		{
			LOG_ERROR() << "Could not read integral value from command";
		}
		else
		{
			LOG_INFO() << "Writing into Register \'" << registerName << "\': value = 0x" << std::hex << value << " = " << std::dec << value;
			EccCmd<Proxy>::ecc()->writeRegister(registerName, value);
		}
	}
};

/**
 * This command read the value of a whole register.
 */
template <class Proxy = EccPrx>
struct ReadReg : public EccCmd<Proxy>
{
	ReadReg(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "reg-read", "<register name>", "Reads all fields of a register of the current node") {}
	void exec(std::istream& cmdArgs)
	{
		std::string registerName;
		uint64_t value;
		std::vector<std::string> fieldList;

		cmdArgs >> registerName;
		value = EccCmd<Proxy>::ecc()->readRegister(registerName);

		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		std::cout << std::left << std::setw(25) << std::setfill('-') <<registerName << "-> ";
		std::cout << "0x" << std::hex << std::right << std::setfill('0') << std::setw(8) << value << std::dec << " (" << value << ")" << std::endl;
		EccCmd<Proxy>::ecc()->getListOfFields(registerName, fieldList);
		std::vector<std::string>::iterator itF;
		for ( itF=fieldList.begin() ; itF < fieldList.end(); itF++ )
		{
			std::cout << std::right << std::setfill(' ') << std::setw(25) << *itF << " = " << EccCmd<Proxy>::ecc()->readField(registerName, *itF) << std::endl;
		}
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command reads name of current device of current node.
 */
template <class Proxy = EccPrx>
struct GetDeviceName : public EccCmd<Proxy>
{
	GetDeviceName(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "device-name", "", "Prints name of current device of current node") {}
	void exec(std::istream&)
	{
		std::string name;
		EccCmd<Proxy>::ecc()->getNameOfCurrentDevice(name);
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		std::cout << "Name of current device: '" << name << '\'' << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command adds a device to the current node.
 */
template <class Proxy = EccPrx>
struct AddDevice : public EccCmd<Proxy>
{
	AddDevice(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "device-add", "<name> <policy> <registerWidth_B> <baseAddress>", "Adds a device to the current node") {}
	void exec(std::istream & cmdArgs)
	{
		mdaq::hw::DeviceDescription descr;
		descr.registerWidth=0;
		cmdArgs >> descr.name >> descr.registerAccess >> descr.registerWidth;
		cmdArgs.unsetf(std::ios::dec);
		cmdArgs >> descr.baseAddress;

		if (cmdArgs.fail())
		{
			throw mdaq::utl::CmdException("Invalid command parameters");
		}
		LOG_INFO() << "Adding device '" << descr.name << "': policy='" << descr.registerAccess << "' regWidth=" << descr.registerWidth << " B address=0x" << std::hex << descr.baseAddress;

		EccCmd<Proxy>::ecc()->addDevice(descr);
	}
};


/**
 * This command removes a device from the current node.
 */
template <class Proxy = EccPrx>
struct RemoveDevice : public EccCmd<Proxy>
{
	RemoveDevice(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "device-remove", "<name>", "Removes a device from the current node") {}
	void exec(std::istream & cmdArgs)
	{
		std::string name;
		cmdArgs >> name;

		LOG_INFO() << "Removing device '" << name << '\'';
		EccCmd<Proxy>::ecc()->removeDevice(name);
	}
};

/**
 * This command reads out all the registers and fields within a device.
 */
template <class Proxy = EccPrx>
struct ReadDevice : public EccCmd<Proxy>
{
	ReadDevice(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "device-read", "<device name>", "Reads everything from a device of the current node") {}
	void exec(std::istream& cmdArgs)
	{
		std::string deviceName;
		uint64_t value;
		std::vector<std::string> registerList;
		std::vector<std::string> fieldList;

		cmdArgs >> deviceName;
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		EccCmd<Proxy>::ecc()->connect("device", deviceName);
		EccCmd<Proxy>::ecc()->getListOfRegisters(registerList);
		std::vector<std::string>::iterator it;
		std::cout << std::left << "== " << std::setw(10) << std::setfill(' ') << deviceName << " ========================================" << std::endl;
		for ( it=registerList.begin() ; it < registerList.end(); it++ )
		{
			try
			{
				value = EccCmd<Proxy>::ecc()->readRegister(*it);
				std::cout << std::left << std::setw(25) << std::setfill('-') <<*it << "-> ";
				std::cout << "0x" << std::hex << std::right << std::setfill('0') << std::setw(8) << value << std::dec << " (" << value << ")" << std::endl;
				EccCmd<Proxy>::ecc()->getListOfFields(*it, fieldList);
				std::vector<std::string>::iterator itF;
				for ( itF=fieldList.begin() ; itF < fieldList.end(); itF++ )
				{
					std::cout << std::right << std::setfill(' ') << std::setw(25) << *itF << " = " << EccCmd<Proxy>::ecc()->readField(*it, *itF) << std::endl;
				}
				std::cout << std::setfill(' ');
			}
			catch (const mdaq::utl::CmdException & e)
			{
			}
		}
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command reads out all the registers and fields of all the devices within current node.
 */
template <class Proxy = EccPrx>
struct ReadNode : public EccCmd<Proxy>
{
	ReadNode(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-read", "", "Reads everything from current node") {}
	void exec(std::istream&)
	{
		uint64_t value;
		std::vector<std::string> deviceList;
		std::vector<std::string> registerList;
		std::vector<std::string> fieldList;

		EccCmd<Proxy>::ecc()->getListOfDevices(deviceList);
		std::vector<std::string>::iterator itD;
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		for ( itD=deviceList.begin() ; itD < deviceList.end(); itD++ )
		{
			std::cout << std::endl << std::left << std::setfill(' ');
			std::cout << "== " << std::setw(10) << *itD << " ========================================" << std::endl << std::endl;
			EccCmd<Proxy>::ecc()->connect("device", *itD);
			EccCmd<Proxy>::ecc()->getListOfRegisters(registerList);
			std::vector<std::string>::iterator it;
			for ( it=registerList.begin() ; it < registerList.end(); it++ )
			{
				try
				{
					value = EccCmd<Proxy>::ecc()->readRegister(*it);
					std::cout << std::left << std::setw(25) << std::setfill('-') <<*it << "-> ";
					std::cout << "0x" << std::hex << std::right << std::setfill('0') << std::setw(8) << value << std::dec << " (" << value << ")" << std::endl;
					EccCmd<Proxy>::ecc()->getListOfFields(*it, fieldList);
					std::vector<std::string>::iterator itF;
					for ( itF=fieldList.begin() ; itF < fieldList.end(); itF++ )
					{
						std::cout << std::right << std::setfill(' ') << std::setw(25) << *itF << " = " << EccCmd<Proxy>::ecc()->readField(*it, *itF) << std::endl;
					}
				}
				catch (const mdaq::utl::CmdException & e)
				{
				}
			}
		}
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command reads name of current node.
 */
template <class Proxy = EccPrx>
struct GetNodeEndpoint : public EccCmd<Proxy>
{
	GetNodeEndpoint(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-endpoint", "", "Get endpoint (IP address and port) used for connecting to the current node") {}
	void exec(std::istream&)
	{
		std::string descr;
		EccCmd<Proxy>::ecc()->getEndpointOfCurrentNode(descr);
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		std::cout << "Endpoint of current node: '" << descr << '\'' << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command reads name of current node.
 */
template <class Proxy = EccPrx>
struct GetNodeName : public EccCmd<Proxy>
{
	GetNodeName(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-name", "", "Get name of current node") {}
	void exec(std::istream&)
	{
		std::string descr;
		EccCmd<Proxy>::ecc()->getNameOfCurrentNode(descr);
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		std::cout << "Name of current node: '" << descr << '\'' << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command starts the data acquisition.
 */
template <class Proxy = EccPrx>
struct DaqStart : public EccCmd<Proxy>
{
	DaqStart(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "daq-start", "", "Launches data acquisition") {}
	void exec(std::istream&)
	{
		EccCmd<Proxy>::ecc()->daqStart();
	}
};

/**
 * This command stops the data acquisition.
 */
template <class Proxy = EccPrx>
struct DaqStop : public EccCmd<Proxy>
{
	DaqStop(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "daq-stop", "", "Stops data acquisition") {}
	void exec(std::istream&)
	{
		EccCmd<Proxy>::ecc()->daqStop();
	}
};

/**
 * This command connects client to specified target.
 */
template <class Proxy = EccPrx>
struct Connect : public EccCmd<Proxy>
{
	Connect(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine,
		"connect",
		"<node|device|daq-ice|daq-tcp|daq-zbuf|daq-fdt> <target id> ",
		"Connects to a node, device or DAQ (of type 'ICE', 'TCP', 'ZBUF' or 'FDT') target"
		) {}
	void exec(std::istream& cmdArgs)
	{
		std::string objectType;
		std::string objectName;

		cmdArgs >> objectType;
		cmdArgs >> objectName;
		EccCmd<Proxy>::ecc()->connect(objectType, objectName);
	}
};

/**
 * This command selects the current device of the current node.
 */
template <class Proxy = EccPrx>
struct SelectDevice : public EccCmd<Proxy>
{
	SelectDevice(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine,
		"device-select",
		"<name>",
		"Selects a device of the current node"
		) {}
	void exec(std::istream & cmdArgs)
	{
		std::string deviceName;

		cmdArgs >> deviceName;
		EccCmd<Proxy>::ecc()->connect("device", deviceName);
	}
};

/**
 * This command reboots the specified target (connection is then lost).
 */
template <class Proxy = EccPrx>
struct Reboot : public EccCmd<Proxy>
{
	Reboot(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-reboot", "", "Reboots the current hardware node") {}
	void exec(std::istream &)
	{
		EccCmd<Proxy>::ecc()->rebootHardware();
	}
};

/**
 * This command prints the list of known nodes.
 */
template <class Proxy = EccPrx>
struct ListNodes : public EccCmd<Proxy>
{
	ListNodes(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-list", "", "Prints the names of known hardware nodes") {}
	void exec(std::istream &)
	{
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		hw::NodeList nodes;
		EccCmd<Proxy>::ecc()->getListOfNodes(nodes);
		std::cout << "Nodes:";
		for (hw::NodeList::const_iterator nameIt = nodes.begin(); nameIt != nodes.end(); ++nameIt)
		{
			std::cout << " '" << *nameIt << '\'';
		}
		std::cout << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command selects a node by its name and makes it the current node.
 */
template <class Proxy = EccPrx>
struct SelectNode : public EccCmd<Proxy>
{
	SelectNode(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-select", "<name>", "Selects a node by its name") {}
	void exec(std::istream & cmdArgs)
	{
		std::string nodeName;
		cmdArgs >> nodeName;

		EccCmd<Proxy>::ecc()->selectNodeByName(nodeName);
	}
};

/**
 * This command prints the list of devices in the current node.
 */
template <class Proxy = EccPrx>
struct ListDevices : public EccCmd<Proxy>
{
	ListDevices(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "device-list", "", "Lists the devices in the current node") {}
	void exec(std::istream &)
	{
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		hw::DeviceList devices;
		EccCmd<Proxy>::ecc()->getListOfDevices(devices);
		std::cout << "Devices:";
		for (hw::DeviceList::const_iterator nameIt = devices.begin(); nameIt != devices.end(); ++nameIt)
		{
			std::cout << " '" << *nameIt << '\'';
		}
		std::cout << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command prints the list of registers in the current device of the current node.
 */
template <class Proxy = EccPrx>
struct ListRegisters : public EccCmd<Proxy>
{
	ListRegisters(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "reg-list", "", "Lists the registers in the current device of the current node") {}
	void exec(std::istream &)
	{
		std::ios::fmtflags flags(std::cout.flags()); // saving cout flags
		hw::RegisterList registers;
		EccCmd<Proxy>::ecc()->getListOfRegisters(registers);
		std::cout << "Registers:";
		for (hw::RegisterList::const_iterator nameIt = registers.begin(); nameIt != registers.end(); ++nameIt)
		{
			std::cout << " '" << *nameIt << '\'';
		}
		std::cout << std::endl;
		std::cout.flags(flags); // restoring cout flags
	}
};

/**
 * This command renames the current node.
 */
template <class Proxy = EccPrx>
struct RenameNode : public EccCmd<Proxy>
{
	RenameNode(typename EccCmd<Proxy>::Engine & engine) : EccCmd<Proxy>(engine, "node-rename", "<new-name>", "Renames the current node") {}
	void exec(std::istream & cmdArgs)
	{
		std::string newName;
		cmdArgs >> newName;
		EccCmd<Proxy>::ecc()->renameNode(newName);
	}
};

template <class Proxy>
EccCmd<Proxy>::Engine::Engine(const std::string & name, const std::string & servantName) : cli::ClientCmd<Proxy>::Engine(name, servantName)
{
	cli::Cmd::Engine::subscribe< Load<Proxy> >();

	cli::Cmd::Engine::subscribe< ListFields<Proxy> >();
	cli::Cmd::Engine::subscribe< WriteField<Proxy> >();
	cli::Cmd::Engine::subscribe< ReadField<Proxy> >();

	cli::Cmd::Engine::subscribe< AddRegister<Proxy> >();
	cli::Cmd::Engine::subscribe< ListRegisters<Proxy> >();
	cli::Cmd::Engine::subscribe< WriteReg<Proxy> >();
	cli::Cmd::Engine::subscribe< ReadReg<Proxy> >();

	cli::Cmd::Engine::subscribe< AddDevice<Proxy> >();
	cli::Cmd::Engine::subscribe< RemoveDevice<Proxy> >();
	cli::Cmd::Engine::subscribe< GetDeviceName<Proxy> >();
	cli::Cmd::Engine::subscribe< ListDevices<Proxy> >();
	cli::Cmd::Engine::subscribe< ReadDevice<Proxy> >();

	cli::Cmd::Engine::subscribe< GetNodeName<Proxy> >();
	cli::Cmd::Engine::subscribe< GetNodeEndpoint<Proxy> >();
	cli::Cmd::Engine::subscribe< ListNodes<Proxy> >();
	cli::Cmd::Engine::subscribe< SelectNode<Proxy> >();
	cli::Cmd::Engine::subscribe< RenameNode<Proxy> >();
	cli::Cmd::Engine::subscribe< ReadNode<Proxy> >();

	cli::Cmd::Engine::subscribe< DaqStart<Proxy> >();
	cli::Cmd::Engine::subscribe< DaqStop<Proxy> >();

	cli::Cmd::Engine::subscribe< Connect<Proxy> >();
	cli::Cmd::Engine::subscribe< SelectDevice<Proxy> >();

	cli::Cmd::Engine::subscribe< Reboot<Proxy> >();

	// Aliases (for backward compatibility)
	cli::Cmd::Engine::addAlias("reboot", "node-reboot");
	cli::Cmd::Engine::addAlias("readreg", "reg-read");
	cli::Cmd::Engine::addAlias("writereg", "reg-write");
	cli::Cmd::Engine::addAlias("readdevice", "device-read");
	cli::Cmd::Engine::addAlias("readnode", "node-read");
	cli::Cmd::Engine::addAlias("daqstart", "daq-start");
	cli::Cmd::Engine::addAlias("daqstop", "daq-stop");
	cli::Cmd::Engine::addAlias("write", "field-write");
	cli::Cmd::Engine::addAlias("read", "field-read");
	cli::Cmd::Engine::addAlias("quit", "q");
	cli::Cmd::Engine::addAlias("exit", "q");

}

template <class Proxy>
typename cli::ClientCmd<Proxy>::IceProxy & EccCmd< Proxy >::Engine::servant()
{
	if (this->ip.port() == 0)
	{
		this->ip.port() = ::mdaq::Default::eccServerCtrlPortNum;
	}
	return cli::ClientCmd<Proxy>::Engine::servant();
}

} // namespace mdaq

#endif // mdaq_EccCmd_hpp_INCLUDED

