/* =============================================================================
 * DeviceFactoryTest.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::hw::server::TestDynDevice
 * Created on: 24 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "mdaq/hw/server/DeviceFactory.h"
#include <boost/test/auto_unit_test.hpp>
#include <boost/version.hpp>

namespace mdaq
{
namespace hw
{
namespace server
{
//____________________________________________________________________________//
/**
 * This Device class is only used for testing purpose to check if the correct
 * factory method is called.
 */
class TestDynDevice : public BaseDynDevice
{
public:
	TestDynDevice(const std::string& name) : _name(name) {}
	const std::string& name() const { return _name; }

	void	registerList(std::vector<std::string>& /*list*/) const {}
	size_t	registerCount() const { return 0; }
	void	registerClear()	{}
	BaseDynRegister& registerCreate	(const std::string& /*registerName*/, uint64_t /*addressOffset*/) { throw std::runtime_error("TestDynDevice does not support creating registers."); }
	BaseDynRegister* registerSearch	(const std::string& /*registerName*/) { return 0; };
	BaseDynRegister& registerFind	(const std::string& registerName) { throw RegisterNotFound(registerName); }

private:
	std::string _name;
};

//____________________________________________________________________________//

/**
 * This is the fixture class associated with the test suite for DeviceFactory.
 * It initializes the factory with several methods.
 */
class DeviceFactoryFixture
{
public:
	DeviceFactoryFixture()
		: factory(DeviceFactory::instance())
	{
		factory.registerFactoryMethod(DeviceCreationParams(4, "Mem"), DeviceFactoryFixture::createDeviceTest1);
		factory.registerFactoryMethod(DeviceCreationParams(8, "Mem"), DeviceFactoryFixture::createDeviceTest2);
		factory.registerFactoryMethod(DeviceCreationParams(4, "CanBus"), DeviceFactoryFixture::createDeviceTest3);
	}

	~DeviceFactoryFixture()
	{
		factory.reset();
	}

	static BaseDynDevicePtr createDeviceTest1(const uint64_t&)
	{
		return BaseDynDevicePtr(new TestDynDevice("Device:4:Mem"));
	}

	static BaseDynDevicePtr createDeviceTest1_aux(const uint64_t&)
	{
		return BaseDynDevicePtr(new TestDynDevice("Device:4:Mem:aux"));
	}

	static BaseDynDevicePtr createDeviceTest2(const uint64_t&)
	{
		return BaseDynDevicePtr(new TestDynDevice("Device:8:Mem"));
	}

	static BaseDynDevicePtr createDeviceTest3(const uint64_t&)
	{
		return BaseDynDevicePtr(new TestDynDevice("Device:4:CanBus"));
	}

	DeviceFactory& factory;
};

// This macro registers a factory method for a device with the specified parameters.
//REGISTER_DEVICE(MemBus4, 4, "Mem", DeviceFactoryFixture::createDeviceTest1)

//____________________________________________________________________________//

BOOST_FIXTURE_TEST_SUITE(DeviceFactoryTestSuite, DeviceFactoryFixture)

BOOST_AUTO_TEST_CASE(createDevice_success1)
{
	BaseDynDevicePtr device = factory.createDevice(DeviceCreationParams(4, "Mem"),0);
	BOOST_CHECK_EQUAL(dynamic_cast<TestDynDevice&>(*device).name(), "Device:4:Mem");
}

BOOST_AUTO_TEST_CASE(createDevice_success2)
{
	BaseDynDevicePtr device = factory.createDevice(DeviceCreationParams(8, "Mem"),0);
	BOOST_CHECK_EQUAL(dynamic_cast<TestDynDevice&>(*device).name(), "Device:8:Mem");
}

BOOST_AUTO_TEST_CASE(createDevice_success3)
{
	BaseDynDevicePtr device = factory.createDevice(DeviceCreationParams(4, "CanBus"),0);
	BOOST_CHECK_EQUAL(dynamic_cast<TestDynDevice&>(*device).name(), "Device:4:CanBus");
}

BOOST_AUTO_TEST_CASE(createDevice_replace)
{
	factory.registerFactoryMethod(DeviceCreationParams(4, "Mem"), DeviceFactoryFixture::createDeviceTest1_aux);
	BaseDynDevicePtr device = factory.createDevice(DeviceCreationParams(4, "Mem"),0);
	BOOST_CHECK_EQUAL(dynamic_cast<TestDynDevice&>(*device).name(), "Device:4:Mem:aux");
}

BOOST_AUTO_TEST_CASE(createDevice_failure1)
{
	BOOST_CHECK_THROW(factory.createDevice(DeviceCreationParams(8, "CanBus"),0), DeviceFactory::UnsupportedDevice);
}

BOOST_AUTO_TEST_CASE(createDevice_failure2)
{
	BOOST_CHECK_THROW(factory.createDevice(DeviceCreationParams(4, "Memm"),0), DeviceFactory::UnsupportedDevice);
}

BOOST_AUTO_TEST_SUITE_END()
#ifdef VXWORKS
;
#endif

//____________________________________________________________________________//

}
}
}
