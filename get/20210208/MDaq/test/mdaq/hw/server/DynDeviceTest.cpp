/* =============================================================================
 * DynDeviceTest.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::hw::server::DynDeviceFixture
 * Created on: 27 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "mdaq/hw/server/DynDevice.hpp"
#include <boost/test/auto_unit_test.hpp>
#include <boost/version.hpp>

namespace mdaq
{
namespace hw
{
namespace server
{
typedef DynDevice<uint32_t, NestedStoragePolicy> TestDevice;
typedef TestDevice::Register Register;

//____________________________________________________________________________//

class DynDeviceFixture
{
public:
	DynDeviceFixture() : device(0)
	{
		device.registerCreate("foo", 0x1000).setValue(0xABCD);
		device.registerCreate("bar", 0x1004).setValue(0xCDEF);
	}
	~DynDeviceFixture()
	{
	}

	TestDevice device;
};

//____________________________________________________________________________//

BOOST_FIXTURE_TEST_SUITE(DynDeviceTestSuite, DynDeviceFixture)

BOOST_AUTO_TEST_CASE(registerCount_test)
{
	BOOST_CHECK_EQUAL(device.registerCount(), 2u);
}

BOOST_AUTO_TEST_CASE(registerSearch_testSuccess)
{
	Register* reg = device.registerSearch("foo");
	BOOST_CHECK(reg != 0);
	BOOST_CHECK_EQUAL(reg->value(), 0xABCDu);
}

BOOST_AUTO_TEST_CASE(registerSearch_testFailure)
{
	Register* reg = device.registerSearch("fooo");
	BOOST_CHECK(reg == 0);
}

BOOST_AUTO_TEST_CASE(registerFind_testSuccess)
{
	BOOST_CHECK_EQUAL(device.registerFind("bar").value(), 0xCDEFu);
}

BOOST_AUTO_TEST_CASE(registerFind_testFailure)
{
	BOOST_CHECK_THROW(device.registerFind("foobar"), BaseDynDevice::RegisterNotFound);
}

BOOST_AUTO_TEST_CASE(operatorSqBrackets_test)
{
	BOOST_CHECK_EQUAL(device["bar"].value(), 0xCDEFu);
}

BOOST_AUTO_TEST_CASE(registerClear_test)
{
	device.registerClear();
	BOOST_CHECK_EQUAL(device.registerCount(), 0u);
}

BOOST_AUTO_TEST_CASE(registerList_test)
{
	std::vector<std::string> list;
	device.registerList(list);
	BOOST_CHECK_EQUAL(list.size(), 2u);
	BOOST_CHECK_EQUAL(list[0], "bar");
	BOOST_CHECK_EQUAL(list[1], "foo");
}

BOOST_AUTO_TEST_SUITE_END()
#ifdef VXWORKS
;
#endif
//____________________________________________________________________________//

}
}
}
