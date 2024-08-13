/* =============================================================================
 * SmartRegisterTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include "mdaq/hw/StdRegister.hpp"
#include <boost/test/auto_unit_test.hpp>
#include <iostream>

namespace mdaq
{
namespace hw
{

BOOST_AUTO_TEST_CASE(address_test)
{
	uint16_t mem;
	RegisterUINT16 reg(&mem);
	BOOST_CHECK_EQUAL(&mem, reg.address());
}

BOOST_AUTO_TEST_CASE(value_test)
{
	uint16_t mem;
	RegisterUINT16 reg(&mem);
	reg.setValue(0);
	BOOST_CHECK_EQUAL(mem, reg.value());
	BOOST_CHECK_EQUAL(mem, 0);
	reg.setValue(0xFFFF);
	BOOST_CHECK_EQUAL(mem, reg.value());
	BOOST_CHECK_EQUAL(mem, 0xFFFFu);
}

BOOST_AUTO_TEST_CASE(assignment_test)
{
	uint8_t mem;
	RegisterUINT8 reg(&mem);
	reg = 0xAB;
	BOOST_CHECK_EQUAL(reg.value(), 0xABu);
}

BOOST_AUTO_TEST_CASE(binaryand_test)
{
	uint32_t mem = 0xABCDEF;
	RegisterUINT32 reg(&mem);
	reg &= 0xFEDCBA;
	BOOST_CHECK_EQUAL(reg.value(), 0xABCDEFu & 0xFEDCBAu);
}

BOOST_AUTO_TEST_CASE(binaryor_test)
{
	uint32_t mem = 0xABCDEF;
	RegisterUINT32 reg(&mem);
	reg |= 0xFEDCBA;
	BOOST_CHECK_EQUAL(reg.value(), 0xABCDEFu | 0xFEDCBAu);
}

BOOST_AUTO_TEST_CASE(binaryxor_test)
{
	uint32_t mem = 0xABCDEF;
	RegisterUINT32 reg(&mem);
	reg ^= 0xFEDCBA;
	BOOST_CHECK_EQUAL(reg.value(), 0xABCDEFu ^ 0xFEDCBAu);
}

BOOST_AUTO_TEST_CASE(leftshift_test)
{
	uint32_t mem = 0x10;
	RegisterUINT32 reg(&mem);
	reg <<= 2;
	BOOST_CHECK_EQUAL(reg.value(), 0x10u << 2);
}

BOOST_AUTO_TEST_CASE(rightshift_test)
{
	uint32_t mem = 0x10;
	RegisterUINT32 reg(&mem);
	reg >>= 2;
	BOOST_CHECK_EQUAL(reg.value(), 0x10u >> 2);
}

BOOST_AUTO_TEST_CASE(unaryadd_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	reg += 0x7890;
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u + 0x7890);
}

BOOST_AUTO_TEST_CASE(unarysub_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	reg -= 0x7890;
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u - 0x7890);
}

BOOST_AUTO_TEST_CASE(unarydiv_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	reg /= 5;
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u / 5);
}

BOOST_AUTO_TEST_CASE(unarymod_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	reg %= 5;
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u % 5);
}

BOOST_AUTO_TEST_CASE(preincrement_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	uint32_t val = ++reg;
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u + 1);
	BOOST_CHECK_EQUAL(reg.value(), val);
}

BOOST_AUTO_TEST_CASE(predecrement_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	uint32_t val = --reg;
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u - 1);
	BOOST_CHECK_EQUAL(reg.value(), val);
}

BOOST_AUTO_TEST_CASE(postincrement_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	uint32_t val = reg++;
	BOOST_CHECK_EQUAL(val, 0x123456u);
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u + 1);
}

BOOST_AUTO_TEST_CASE(postdecrement_test)
{
	uint32_t mem = 0x123456;
	RegisterUINT32 reg(&mem);
	uint32_t val = reg--;
	BOOST_CHECK_EQUAL(val, 0x123456u);
	BOOST_CHECK_EQUAL(reg.value(), 0x123456u - 1);
}

}
}
