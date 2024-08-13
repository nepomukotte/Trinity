/* =============================================================================
 * DynRegisterTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 6 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "mdaq/hw/server/DynRegister.hpp"
#include "mdaq/hw/NestedStoragePolicy.hpp"
#include <boost/test/auto_unit_test.hpp>
#include <boost/version.hpp>

namespace mdaq
{
namespace hw
{
namespace server
{
typedef DynRegister<uint32_t, NestedStoragePolicy> Register;

//____________________________________________________________________________//

class DynRegisterFixture
{
public:
	DynRegisterFixture() : reg("scReg"), baseReg(reg)
	{
		reg.fieldAdd("isActive", 0, 1);
		reg.fieldAdd("threshold", 1, 12);
		reg.fieldAdd("boost", 13, 1);
		reg.fieldAdd("amplification", 14, 18);
	}

	~DynRegisterFixture()
	{
	}

	Register reg;
	BaseDynRegister& baseReg;
};

//____________________________________________________________________________//

BOOST_FIXTURE_TEST_SUITE(DynRegisterTestSuite, DynRegisterFixture)

BOOST_AUTO_TEST_CASE(fieldCount_test)
{
	BOOST_CHECK_EQUAL(reg.fieldCount(), 4u);
}

BOOST_AUTO_TEST_CASE(fieldRead_test)
{
	reg.setValue(1);
	BOOST_CHECK_EQUAL(reg.fieldRead("isActive"), 1u);

	reg.setValue(50 << 1);
	BOOST_CHECK_EQUAL(reg.fieldRead("isActive"), 0u);
	BOOST_CHECK_EQUAL(reg.fieldRead("threshold"), 50u);
	BOOST_CHECK_EQUAL(reg.fieldRead("boost"), 0u);
	BOOST_CHECK_EQUAL(reg.fieldRead("amplification"), 0u);

	reg.setValue(500u << 1 | 1u << 13 | 10u << 14);
	BOOST_CHECK_EQUAL(reg.fieldRead("isActive"), 0u);
	BOOST_CHECK_EQUAL(reg.fieldRead("threshold"), 500u);
	BOOST_CHECK_EQUAL(reg.fieldRead("boost"), 1u);
	BOOST_CHECK_EQUAL(reg.fieldRead("amplification"), 10u);
}

bool isNotFoundMessageOk_isActived(const BaseDynRegister::FieldNotFound& e)
{
	return e.getMessage().compare("No field is named \"isActived\"") == 0;
}

BOOST_AUTO_TEST_CASE(fieldRead_testNotFound)
{
	BOOST_CHECK_EXCEPTION(reg.fieldRead("isActived"), BaseDynRegister::FieldNotFound, isNotFoundMessageOk_isActived);
}

BOOST_AUTO_TEST_CASE(fieldWrite_test)
{
	reg.fieldWrite("isActive", 1);
	reg.fieldWrite("threshold", 4095);
	reg.fieldWrite("amplification", 20000);
	BOOST_CHECK_EQUAL(reg.value(), 1u | 4095u << 1 | 20000u << 14);
}

bool isOutOfRangeMessageOk_isActive(const Register::OutOfRange& e)
{
	return e.getMessage().compare("Error setting the field \"isActive\" value to 2 (valid range is 0 to 1)") == 0;
}

bool isOutOfRangeMessageOk_threshold(const Register::OutOfRange& e)
{
	return e.getMessage().compare("Error setting the field \"threshold\" value to 4096 (valid range is 0 to 4095)") == 0;
}

BOOST_AUTO_TEST_CASE(fieldWrite_testOutOfRange)
{

	BOOST_CHECK_EXCEPTION(reg.fieldWrite("isActive", 2), Register::OutOfRange, isOutOfRangeMessageOk_isActive);
	BOOST_CHECK_EXCEPTION(reg.fieldWrite("threshold", 4096), Register::OutOfRange, isOutOfRangeMessageOk_threshold);
}

BOOST_AUTO_TEST_CASE(fieldWrite_testNotFound)
{
	BOOST_CHECK_THROW(reg.fieldWrite("Amplification", 20), BaseDynRegister::FieldNotFound);
}

BOOST_AUTO_TEST_CASE(getValueInt64_test)
{
	reg.setValue(0x7ABCDEF0u);
	BOOST_CHECK_EQUAL(baseReg.getValueInt64(), 0x7ABCDEF0u);
}

BOOST_AUTO_TEST_CASE(setValueInt64_test)
{
	baseReg.setValueInt64(0x789ABCDEu);
	BOOST_CHECK_EQUAL(reg.value(), 0x789ABCDEu);
}

BOOST_AUTO_TEST_CASE(fieldWriteInt64_test)
{
	baseReg.fieldWriteInt64("isActive", 1);
	baseReg.fieldWriteInt64("threshold", 4095);
	baseReg.fieldWriteInt64("amplification", 20000);
	BOOST_CHECK_EQUAL(reg.value(), 1u | 4095u << 1 | 20000u << 14);
}

BOOST_AUTO_TEST_CASE(fieldList_test)
{
	std::vector<std::string> list;
	baseReg.fieldList(list);
	BOOST_CHECK_EQUAL(list.size(), 4u);
	BOOST_CHECK_EQUAL(list[0], "amplification");
	BOOST_CHECK_EQUAL(list[1], "boost");
	BOOST_CHECK_EQUAL(list[2], "isActive");
	BOOST_CHECK_EQUAL(list[3], "threshold");
}

BOOST_AUTO_TEST_SUITE_END()
#ifdef VXWORKS
;
#endif

//____________________________________________________________________________//

}
}
}
