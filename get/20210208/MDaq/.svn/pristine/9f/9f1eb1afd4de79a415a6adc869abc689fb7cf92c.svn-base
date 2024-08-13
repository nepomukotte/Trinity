/* =============================================================================
 * HwNodeTest.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::hw::server::HwNodeFixture
 * Created on: 25 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the test classes the Mordicus Real-Time Software
 * Development Framework.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#include "mdaq/hw/Control.h"
#include "mdaq/hw/RemoteRegisterProcessor.h"
#include <utl/Logging.h>
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Initialize.h)
#include <boost/test/auto_unit_test.hpp>
#include <boost/lexical_cast.hpp>

#ifndef VXWORKS
#	include <boost/cstdint.hpp>
#endif

namespace mdaq
{
namespace hw
{
namespace server
{
//____________________________________________________________________________//

class HwNodeFixture
{
public:
	HwNodeFixture()
	{
		LOG_INFO() << std::string(120, '.');
		scReg[0] = 0;
		scReg[1] = 0;

		BOOST_REQUIRE_NO_THROW(ic = Ice::initialize(
				boost::unit_test::framework::master_test_suite().argc,
				boost::unit_test::framework::master_test_suite().argv));

		Ice::ObjectPrx base = ic->stringToProxy("HwNode:default -p 10000");
		BOOST_REQUIRE(base != 0);
		BOOST_REQUIRE_NO_THROW(node = hw::NodePrx::checkedCast(base));
		BOOST_REQUIRE(node != 0);

		if (instanceCount++ <= 0)
		{
			LOG_INFO() << std::string(120, '=');
			node->setName("HwServer-Test");
		}

		initNode(nodeConfig, "ScottTestBench", 1);
		DeviceConfig& device = initDevice(nodeConfig.devices[0], "scott", 0x1000, "MemBus", 4, 2);
		RegisterConfig& reg1 = initRegister(device.registers[0], "scReg1", 0, false, 5);
		initField(reg1.fields[0], "synchro_inv",	0, 1);
		initField(reg1.fields[1], "out_resync",		1, 1);
		initField(reg1.fields[2], "reset_dll",		2, 1);
		initField(reg1.fields[3], "reset_mem",		3, 1);
		initField(reg1.fields[4], "en1to8",			4, 1);

		RegisterConfig& reg2 = initRegister(device.registers[1], "scReg2", 4, false, 6);
		initField(reg2.fields[0], "elatch_L0",		0, 1);
		initField(reg2.fields[1], "test_mem_out",	1, 3);
		initField(reg2.fields[2], "ext_L0",			4, 1);
		initField(reg2.fields[3], "L0_resync",		5, 2);
		initField(reg2.fields[4], "L0_latency",		7, 4);
		initField(reg2.fields[5], "en_dll_out",		11, 1);

		BOOST_REQUIRE_NO_THROW(node->create(nodeConfig));
		BOOST_REQUIRE_NO_THROW(scott = node->findDevice("scott"));
	}
	~HwNodeFixture()
	{
		node->removeAllDevices();
		if(ic)
		{
			ic->destroy();
		}
		//LOG_DEBUG() << std::string(120, '.');
	}

	NodeConfig& initNode(NodeConfig& node, const char* name, size_t deviceCount)
	{
		node.id = name;
		node.devices.resize(deviceCount);
		return node;
	}
	DeviceConfig& initDevice(DeviceConfig& dev, const char* name, Ice::Long baseAddress, const char* access, unsigned regWidth, size_t regCount)
	{
		dev.descr.name = name;
		dev.descr.baseAddress = baseAddress;
		dev.descr.registerAccess = access;
		dev.descr.registerWidth = regWidth;
		dev.registers.resize(regCount);
		return dev;
	}
	RegisterConfig& initRegister(RegisterConfig& reg, const char* name, int64_t offset, bool /*isReadOnly*/, size_t fieldCount)
	{
		reg.descr.name = name;
		reg.descr.offset = offset;
		reg.fields.resize(fieldCount);
		return reg;
	}
	FieldDescription& initField(FieldDescription& field, const char* name, short offset, short width)
	{
		field.name = name;
		field.offset = offset;
		field.width = width;
		return field;
	}

public:
	Ice::CommunicatorPtr ic;
	NodeConfig nodeConfig;
	hw::NodePrx node;
	DevicePrx scott;
	uint32_t scReg[2];

private:
	static int instanceCount;
};

int HwNodeFixture::instanceCount = 0;
//____________________________________________________________________________//

BOOST_FIXTURE_TEST_SUITE(HwNodeTestSuite, HwNodeFixture)

BOOST_AUTO_TEST_CASE(name_test)
{
	BOOST_CHECK_EQUAL(node->name(), "ScottTestBench");
}

BOOST_AUTO_TEST_CASE(enumerateDevices_test)
{
	DeviceList devices;
	BOOST_CHECK_NO_THROW(node->getListOfDevices(devices));
	BOOST_CHECK_EQUAL(devices.size(), 1u);
	BOOST_CHECK_EQUAL(devices[0], "scott");
}

BOOST_AUTO_TEST_CASE(findDevice_test)
{
	DevicePrx device;
	BOOST_CHECK_NO_THROW(device = node->findDevice("scott"));
	BOOST_CHECK(device != 0);
}

BOOST_AUTO_TEST_CASE(findDevice_testHardwareNotFound)
{
	BOOST_CHECK_THROW(node->findDevice("foo"), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(deviceName_test)
{
	BOOST_CHECK_EQUAL(scott->name(), "scott");
}

BOOST_AUTO_TEST_CASE(deviceType_test)
{
	BOOST_CHECK_EQUAL(scott->registerAccess(), "MemBus");
}

BOOST_AUTO_TEST_CASE(deviceRegisterSize_test)
{
	BOOST_CHECK_EQUAL(scott->registerWidth(), 4);
}

BOOST_AUTO_TEST_CASE(writeRegister_test1)
{
	BOOST_CHECK_NO_THROW(scott->writeRegister("scReg1", 0x1234));
	BOOST_CHECK_EQUAL(scott->readRegister("scReg1"), 0x1234u);
}

BOOST_AUTO_TEST_CASE(writeRegister_testNotFound)
{
	BOOST_CHECK_THROW(scott->writeRegister("scReg42", 0x1234), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(readRegister_testNotFound)
{
	BOOST_CHECK_THROW(scott->readRegister("scReg42"), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(writeRegister_testOutOfRange)
{
	BOOST_CHECK_THROW(scott->writeRegister("scReg1", 0x1234567890LL), ::mdaq::utl::CmdException);
	BOOST_CHECK_THROW(scott->writeRegister("scReg1", -1), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(readField_test1)
{
	BOOST_CHECK_NO_THROW(scott->writeRegister("scReg1", 1u << 4 | 1 << 2 | 1));
	BOOST_CHECK_EQUAL(scott->readField("scReg1", "en1to8"), 1u);
	BOOST_CHECK_EQUAL(scott->readField("scReg1", "reset_mem"), 0u);
	BOOST_CHECK_EQUAL(scott->readField("scReg1", "reset_dll"), 1u);
	BOOST_CHECK_EQUAL(scott->readField("scReg1", "out_resync"), 0u);
	BOOST_CHECK_EQUAL(scott->readField("scReg1", "synchro_inv"), 1u);
}

BOOST_AUTO_TEST_CASE(readField_test2)
{
	BOOST_CHECK_NO_THROW(scott->writeRegister("scReg2", 1u << 11 | 9 << 7 | 2 << 5 | 1 << 4 | 5 << 1));
	BOOST_CHECK_EQUAL(scott->readField("scReg2", "en_dll_out"), 1u);
	BOOST_CHECK_EQUAL(scott->readField("scReg2", "L0_latency"), 9u);
	BOOST_CHECK_EQUAL(scott->readField("scReg2", "L0_resync"), 2u);
	BOOST_CHECK_EQUAL(scott->readField("scReg2", "ext_L0"), 1u);
	BOOST_CHECK_EQUAL(scott->readField("scReg2", "test_mem_out"), 5u);
	BOOST_CHECK_EQUAL(scott->readField("scReg2", "elatch_L0"), 0u);
}

BOOST_AUTO_TEST_CASE(readField_testNotFound)
{
	BOOST_CHECK_THROW(scott->readField("scReg2", "L0latency"), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(writeField_test1)
{
	BOOST_CHECK_NO_THROW(scott->writeField("scReg1", "en1to8", 1));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg1", "reset_mem", 0));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg1", "reset_dll", 1));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg1", "out_resync", 0));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg1", "synchro_inv", 1));
	BOOST_CHECK_EQUAL(scott->readRegister("scReg1"), 1u << 4 | 1 << 2 | 1);
}

BOOST_AUTO_TEST_CASE(writeField_test2)
{
	BOOST_CHECK_NO_THROW(scott->writeField("scReg2", "en_dll_out", 1));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg2", "L0_latency", 9));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg2", "L0_resync", 2));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg2", "ext_L0", 1));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg2", "test_mem_out", 5));
	BOOST_CHECK_NO_THROW(scott->writeField("scReg2", "elatch_L0", 0));
	BOOST_CHECK_EQUAL(scott->readRegister("scReg2"), 1u << 11 | 9 << 7 | 2 << 5 | 1 << 4 | 5 << 1);
}

BOOST_AUTO_TEST_CASE(writeField_testNotFound)
{
	BOOST_CHECK_THROW(scott->writeField("scReg2", "reset_DLL", 1), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(writeField_testOutOfRange)
{
	BOOST_CHECK_THROW(scott->writeField("scReg2", "L0_latency", 16), ::mdaq::utl::CmdException);
	BOOST_CHECK_THROW(scott->writeField("scReg2", "L0_latency", std::numeric_limits<int64_t>::max()), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(enumerateRegisters_test)
{
	RegisterList registers;
	BOOST_CHECK_NO_THROW(scott->getListOfRegisters(registers));
	BOOST_CHECK_EQUAL(registers.size(), 2u);
	BOOST_CHECK_EQUAL(registers[0], "scReg1");
	BOOST_CHECK_EQUAL(registers[1], "scReg2");
}

BOOST_AUTO_TEST_CASE(enumerateFields_test)
{
	FieldList fields;
	BOOST_CHECK_NO_THROW(scott->getListOfFields("scReg1", fields));
	BOOST_CHECK_EQUAL(fields.size(), 5u);
	BOOST_CHECK_EQUAL(fields[0], "en1to8");
	BOOST_CHECK_EQUAL(fields[1], "out_resync");
	BOOST_CHECK_EQUAL(fields[2], "reset_dll");
	BOOST_CHECK_EQUAL(fields[3], "reset_mem");
	BOOST_CHECK_EQUAL(fields[4], "synchro_inv");
}

BOOST_AUTO_TEST_CASE(enumerateFields_testNotFound)
{
	FieldList fields;
	BOOST_CHECK_THROW(scott->getListOfFields("scReg3", fields), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(batchExec_test)
{
	RegCmdSeq cmdBatch, cmdBatchStatus;

	cmdBatch.push_back(makeRegCmd(Write, "scReg1", "en1to8", 1));
	cmdBatch.push_back(makeRegCmd(Write, "scReg1", "reset_mem", 0));
	cmdBatch.push_back(makeRegCmd(Write, "scReg1", "reset_dll", 1));
	cmdBatch.push_back(makeRegCmd(Write, "scReg1", "out_resync", 0));
	cmdBatch.push_back(makeRegCmd(Write, "scReg1", "synchro_inv", 1));
	cmdBatch.push_back(makeRegCmd(Read, "scReg1", "", -1));

	BOOST_CHECK_NO_THROW(scott->execBatch(cmdBatch, cmdBatchStatus));
	BOOST_CHECK_EQUAL(cmdBatchStatus.size(), cmdBatch.size());
	BOOST_CHECK_EQUAL(cmdBatchStatus.back().instruction, Read);
	BOOST_CHECK_EQUAL(cmdBatchStatus.back().regName, "scReg1");
	BOOST_CHECK_EQUAL(cmdBatchStatus.back().fieldName, "");
	BOOST_CHECK_EQUAL(cmdBatchStatus.back().value, 1u << 4 | 1 << 2 | 1);

	cmdBatch.push_back(makeRegCmd(Read, "scReg1", "toto_batchExec", -1));
	BOOST_CHECK_THROW(scott->execBatch(cmdBatch, cmdBatchStatus), ::mdaq::utl::CmdException);
}

BOOST_AUTO_TEST_CASE(remoteProcessorUnbatched_test)
{
	RemoteRegisterProcessor rrp(scott);
	BOOST_CHECK_EQUAL(rrp.isBatched(), false);
	Ice::Long readBack = -1;

	BOOST_CHECK_NO_THROW(
	rrp.write("scReg1", "en1to8", 1)
	   .write("scReg1", "reset_mem", 0)
	   .write("scReg1", "reset_dll", 1)
	   .write("scReg1", "out_resync", 0)
	   .write("scReg1", "synchro_inv", 1)
	   .read("scReg1", readBack)
	);

	BOOST_CHECK_EQUAL(readBack, 1u << 4 | 1 << 2 | 1);

	readBack = -1;
	Ice::Long noReadBack = -1;
	BOOST_CHECK_THROW(
	rrp.write("scReg1", "en1to8", 1)
	   .write("scReg1", "reset_mem", 0)
	   .write("scReg1", "reset_dll", 1)
	   .write("scReg1", "out_resync", 0)
	   .write("scReg1", "synchro_inv", 1)
	   .read("scReg1", readBack)
	   .write("scReg1", "toto_remoteProcessorUnbatched", 42)
	   .read("scReg1", noReadBack)
	, ::mdaq::utl::CmdException
	);

	BOOST_CHECK_EQUAL(readBack, 1u << 4 | 1 << 2 | 1);
	BOOST_CHECK_EQUAL(noReadBack, -1);
}

BOOST_AUTO_TEST_CASE(remoteProcessorBatched_test)
{
	RemoteRegisterProcessor rrp(scott, true);
	BOOST_CHECK_EQUAL(rrp.isBatched(), true);
	Ice::Long readBack = -1;

	BOOST_CHECK_NO_THROW(
	rrp.write("scReg1", "en1to8", 1)
	   .write("scReg1", "reset_mem", 0)
	   .write("scReg1", "reset_dll", 1)
	   .write("scReg1", "out_resync", 0)
	   .write("scReg1", "synchro_inv", 1)
	   .read ("scReg1", readBack)
	);
	BOOST_CHECK_EQUAL(dynamic_cast<RemoteRegisterProcessor::BatchCore&>(rrp.core()).regCmdBatchStatus.size(), 6);
	BOOST_CHECK_EQUAL(readBack, 1u << 4 | 1 << 2 | 1);

	readBack = -1;
	Ice::Long noReadBack = -1;

	BOOST_CHECK_THROW(
	rrp.write("scReg1", "en1to8", 1)
	   .write("scReg1", "reset_mem", 0)
	   .write("scReg1", "reset_dll", 1)
	   .write("scReg1", "out_resync", 0)
	   .write("scReg1", "synchro_inv", 1)
	   .read ("scReg1", readBack)
	   .write("scReg1", "toto_remoteProcessorBatched", 42)
	, ::mdaq::utl::CmdException
	);
	BOOST_CHECK_EQUAL(dynamic_cast<RemoteRegisterProcessor::BatchCore&>(rrp.core()).regCmdBatchStatus.size(), 0);
	BOOST_CHECK_EQUAL(readBack, -1);
	BOOST_CHECK_EQUAL(noReadBack, -1);
}

BOOST_AUTO_TEST_CASE(remoteProcessorSteps_test)
{
	// Unbatched test first
	Ice::Long readBack = -1;
	RemoteRegisterProcessor rrp(scott, false);
	rrp << &RegisterProcessor::init
		<< Cmd<Write>("scReg2", "en_dll_out", 1)
		<< Cmd<Write>("scReg2", "L0_latency", 9)
		<< Cmd<Write>("scReg2", "L0_resync", 2)
		<< Cmd<Write>("scReg2", "ext_L0", 1)
		<< Cmd<Write>("scReg2", "test_mem_out", 5)
		<< Cmd<Write>("scReg2", "elatch_L0", 0)
		<< Cmd<Read >("scReg2", "", &readBack)
		<< &RegisterProcessor::flush;
	BOOST_CHECK_EQUAL(readBack, 1u << 11 | 9 << 7 | 2 << 5 | 1 << 4 | 5 << 1);

	BOOST_CHECK_THROW(
		rrp << &RegisterProcessor::init
	        << Cmd<Write>("scReg1", "en1to8", 1)
	        << Cmd<Write>("scReg1", "reset_mem", 0)
	        << Cmd<Write>("scReg1", "reset_dll", 1)
	        << Cmd<Write>("scReg1", "out_resync", 0)
	        << Cmd<Write>("scReg1", "synchro_inv", 1)
	        << Cmd<Read >("scReg1", "", &readBack)
	        << Cmd<Write>("scReg1", "toto_remoteProcessorBatched", 42)
			<< &RegisterProcessor::flush;
	, ::mdaq::utl::CmdException
	);
	BOOST_CHECK_EQUAL(readBack, 1u << 4 | 1 << 2 | 1);
}

BOOST_AUTO_TEST_CASE(remoteProcessorStepsBacthed_test)
{
	// Unbatched test first
	Ice::Long readBack = -1;
	RemoteRegisterProcessor rrp(scott, true);
	rrp << &RegisterProcessor::init
		<< Cmd<Write>("scReg2", "en_dll_out", 1)
		<< Cmd<Write>("scReg2", "L0_latency", 9)
		<< Cmd<Write>("scReg2", "L0_resync", 2)
		<< Cmd<Write>("scReg2", "ext_L0", 1)
		<< Cmd<Write>("scReg2", "test_mem_out", 5)
		<< Cmd<Write>("scReg2", "elatch_L0", 0)
		<< Cmd<Read >("scReg2", "", &readBack)
		<< &RegisterProcessor::flush;

	BOOST_CHECK_EQUAL(readBack, 1u << 11 | 9 << 7 | 2 << 5 | 1 << 4 | 5 << 1);

	readBack = -1;
	BOOST_CHECK_THROW(
		rrp << &RegisterProcessor::init
	        << Cmd<Write>("scReg1", "en1to8", 1)
	        << Cmd<Write>("scReg1", "reset_mem", 0)
	        << Cmd<Write>("scReg1", "reset_dll", 1)
	        << Cmd<Write>("scReg1", "out_resync", 0)
	        << Cmd<Write>("scReg1", "synchro_inv", 1)
	        << Cmd<Read >("scReg1", "", &readBack)
	        << Cmd<Write>("scReg1", "toto_remoteProcessorBatched", 42)
			<< &RegisterProcessor::flush;
	, ::mdaq::utl::CmdException
	);
	BOOST_CHECK_EQUAL(readBack, -1);
}

BOOST_AUTO_TEST_SUITE_END()

//____________________________________________________________________________//
}
}
}
