/* =============================================================================
 * BitFieldHelperTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 12 nov. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the Utilities_clean project.
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

#include <utl/BitFieldHelper.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/cstdint.hpp>

using namespace utl;

struct BitFieldFixture
{
	uint32_t pattern1;
	uint32_t pattern2;
	int32_t pattern3;
	uint16_t u16_1;

	BitFieldFixture()
		: pattern1(0xAF3C7832),
		  pattern2(0x7FFFFFEE),
		  pattern3(0xFF6CE738),
		  u16_1(0xF742)
	{
	}

	template<unsigned pos>
	void checkSetBit()
	{
		uint32_t var = 0;
		BitFieldHelper<uint32_t>::setBit<pos>(var, true);
		BOOST_CHECK_EQUAL(var, 1 << pos);

		uint32_t var2 = pattern2;
		BitFieldHelper<uint32_t>::setBit<pos>(var2, true);
		BOOST_CHECK_EQUAL(var2, pattern2 | (1 << pos));
	}
};

BOOST_FIXTURE_TEST_SUITE(utl_BitFieldHelper, BitFieldFixture)

BOOST_AUTO_TEST_CASE(testSetBitTpl_u32)
{
	checkSetBit< 0>();
	checkSetBit<31>();
	checkSetBit< 4>();
}

// =============================================================================
// tests cases for getBit methods
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testGetBitTpl_u32)
{
	uint32_t b1 = BitFieldHelper<uint32_t>::getBit<1>(pattern1);
	BOOST_CHECK_EQUAL(b1, 1);
	uint32_t b28 = BitFieldHelper<uint32_t>::getBit<28>(pattern1);
	BOOST_CHECK_EQUAL(b28, 0);
	uint32_t b31 = BitFieldHelper<uint32_t>::getBit<31>(pattern1);
	BOOST_CHECK_EQUAL(b31, 1);
}

BOOST_AUTO_TEST_CASE(testGetBit_u32)
{
	uint32_t b1 = BitFieldHelper<uint32_t>::getBit(pattern1, 1);
	BOOST_CHECK_EQUAL(b1, 1);
	uint32_t b28 = BitFieldHelper<uint32_t>::getBit(pattern1, 28);
	BOOST_CHECK_EQUAL(b28, 0);
	uint32_t b31 = BitFieldHelper<uint32_t>::getBit(pattern1, 31);
	BOOST_CHECK_EQUAL(b31, 1);
}

BOOST_AUTO_TEST_CASE(testGetBitTpl_u16)
{
	uint16_t b0 = BitFieldHelper<uint16_t>::getBit<0>(u16_1);
	BOOST_CHECK_EQUAL(b0, 0);
	uint16_t b15 = BitFieldHelper<uint16_t>::getBit<15>(u16_1);
	BOOST_CHECK_EQUAL(b15, 1);
	uint16_t b9 = BitFieldHelper<uint16_t>::getBit<9>(u16_1);
	BOOST_CHECK_EQUAL(b9, 1);
}
// =============================================================================
// tests cases for getField methods
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testGetFieldTpl_u32)
{
	uint32_t f1 = BitFieldHelper<uint32_t>::getField<4, 8>(pattern1);
	BOOST_CHECK_EQUAL(f1, 0x83);
	uint32_t f2 = BitFieldHelper<uint32_t>::getField<12, 12>(pattern1);
	BOOST_CHECK_EQUAL(f2, 0x3C7);
}

BOOST_AUTO_TEST_CASE(testGetFieldTpl_u16)
{
	uint16_t f1 = BitFieldHelper<uint16_t>::getField<4, 8>(u16_1);
	BOOST_CHECK_EQUAL(f1, 0x74);
	uint16_t f2 = BitFieldHelper<uint16_t>::getField<12, 4>(u16_1);
	BOOST_CHECK_EQUAL(f2, 0xF);
}

BOOST_AUTO_TEST_CASE(testGetFieldTpl_u32_i32)
{
	// checking getting a field as a signed integer
	int32_t f1;
	BitFieldHelper<uint32_t>::getField<16, 12>(pattern1, f1);
	BOOST_CHECK_EQUAL(f1, -196);
}

BOOST_AUTO_TEST_CASE(testGetFieldTpl_u32_i16)
{
	// checking getting a field as a smaller signed integer
	int16_t f1;
	BitFieldHelper<uint32_t>::getField<16, 12>(pattern1, f1);
	BOOST_CHECK_EQUAL(f1, -196);
}

BOOST_AUTO_TEST_CASE(testGetField_u32)
{
	uint32_t f1 = BitFieldHelper<uint32_t>::getField(pattern1, 4, 8);
	BOOST_CHECK_EQUAL(f1, 0x83);
	uint32_t f2 = BitFieldHelper<uint32_t>::getField(pattern1, 12, 12);
	BOOST_CHECK_EQUAL(f2, 0x3C7);
}

BOOST_AUTO_TEST_CASE(testGetField_u16)
{
	uint16_t f1 = BitFieldHelper<uint16_t>::getField(u16_1, 4, 8);
	BOOST_CHECK_EQUAL(f1, 0x74);
	uint16_t f2 = BitFieldHelper<uint16_t>::getField(u16_1, 12, 4);
	BOOST_CHECK_EQUAL(f2, 0xF);
}

BOOST_AUTO_TEST_CASE(testGetField_u32_i32)
{
	// checking getting a field as a signed integer
	int32_t f1;
	BitFieldHelper<uint32_t>::getField(pattern1, 16, 12, f1);
	BOOST_CHECK_EQUAL(f1, -196);
}

BOOST_AUTO_TEST_CASE(testGetField_u32_i16)
{
	// checking getting a field as a smaller signed integer
	int16_t f1;
	BitFieldHelper<uint32_t>::getField(pattern1, 16, 12, f1);
	BOOST_CHECK_EQUAL(f1, -196);
}

BOOST_AUTO_TEST_CASE(testGetFieldTpl_i32)
{
	int32_t f1 = BitFieldHelper<int32_t>::getField<28, 4>(pattern3);
	BOOST_CHECK_EQUAL(f1, -1);
	int32_t f2 = BitFieldHelper<int32_t>::getField<20, 8>(pattern3);
	BOOST_CHECK_EQUAL(f2, -10);
	int32_t f3 = BitFieldHelper<int32_t>::getField<11, 9>(pattern3);
	BOOST_CHECK_EQUAL(f3, -100);
	int32_t f4 = BitFieldHelper<int32_t>::getField<0, 11>(pattern3);
	BOOST_CHECK_EQUAL(f4, -200);
}

BOOST_AUTO_TEST_CASE(testGetFieldTpl_i32_u32)
{
	uint32_t f1;
	BitFieldHelper<int32_t>::getField<11, 9>(pattern3, f1);
	BOOST_CHECK_EQUAL(f1, 0x19C);
}

BOOST_AUTO_TEST_CASE(testGetFieldTpl_i32_u8)
{
	uint8_t f1 = BitFieldHelper<int32_t>::getField<20, 8, uint8_t>(pattern3);
	BOOST_CHECK_EQUAL(f1, 0xF6);
}

BOOST_AUTO_TEST_CASE(testGetField_i32)
{
	int32_t f1 = BitFieldHelper<int32_t>::getField(pattern3, 28, 4);
	BOOST_CHECK_EQUAL(f1, -1);
	int32_t f2 = BitFieldHelper<int32_t>::getField(pattern3, 20, 8);
	BOOST_CHECK_EQUAL(f2, -10);
	int32_t f3 = BitFieldHelper<int32_t>::getField(pattern3, 11, 9);
	BOOST_CHECK_EQUAL(f3, -100);
	int32_t f4 = BitFieldHelper<int32_t>::getField(pattern3, 0, 11);
	BOOST_CHECK_EQUAL(f4, -200);
}

BOOST_AUTO_TEST_CASE(testGetField_i32_u32)
{
	uint32_t f1;
	BitFieldHelper<int32_t>::getField(pattern3, 11, 9, f1);
	BOOST_CHECK_EQUAL(f1, 0x19C);
}

BOOST_AUTO_TEST_CASE(testGetField_i32_u8)
{
	uint8_t f1 = BitFieldHelper<int32_t>::getField<uint8_t>(pattern3, 20, 8);
	BOOST_CHECK_EQUAL(f1, 0xF6);
}

// =============================================================================
// tests cases for setField methods
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testSetFieldTpl_u32)
{
	uint32_t var = 0;
	BitFieldHelper<uint32_t>::setField<0, 5>(var, 0x15);
	BOOST_CHECK_EQUAL(var, 0x15);

	uint32_t var2 = pattern1;
	BitFieldHelper<uint32_t>::setField<0, 5>(var2, 0x15);
	BOOST_CHECK_EQUAL(var2, 0xAF3C7835);

	var = 0;
	BitFieldHelper<uint32_t>::setField<16, 4>(var, 0xD);
	BOOST_CHECK_EQUAL(var, 0x000D0000);

	var2 = pattern1;
	BitFieldHelper<uint32_t>::setField<8, 16>(var2, 0x1234);
	BOOST_CHECK_EQUAL(var2, 0xAF123432);

	// clearing a single bit
	var2 = pattern2;
	BitFieldHelper<uint32_t>::setField<29, 1>(var2, 0);
	BOOST_CHECK_EQUAL(var2, 0x5FFFFFEE);

	// setting a single bit
	var2 = pattern2;
	BitFieldHelper<uint32_t>::setField<31, 1>(var2, 1);
	BOOST_CHECK_EQUAL(var2, 0xFFFFFFEE);

	// checking if field overflow is clipped
	var2 = pattern1;
	BitFieldHelper<uint32_t>::setField<12, 11>(var2, 0x1F55);
	BOOST_CHECK_EQUAL(var2, 0xAF755832);
}

BOOST_AUTO_TEST_CASE(testSetFieldTpl_u32_i16)
{
	uint32_t var = pattern1;
	BitFieldHelper<uint32_t>::setField<16, 12>(var, int16_t(-1234));
	BOOST_CHECK_EQUAL(var, 0xAB2E7832);
}

BOOST_AUTO_TEST_CASE(testSetField_u32)
{
	uint32_t var = 0;
	BitFieldHelper<uint32_t>::setField(var, 0, 5, 0x15);
	BOOST_CHECK_EQUAL(var, 0x15);

	uint32_t var2 = pattern1;
	BitFieldHelper<uint32_t>::setField(var2, 0, 5, 0x15);
	BOOST_CHECK_EQUAL(var2, 0xAF3C7835);

	var = 0;
	BitFieldHelper<uint32_t>::setField(var, 16, 4, 0xD);
	BOOST_CHECK_EQUAL(var, 0x000D0000);

	var2 = pattern1;
	BitFieldHelper<uint32_t>::setField(var2, 8, 16, 0x1234);
	BOOST_CHECK_EQUAL(var2, 0xAF123432);

	// clearing a single bit
	var2 = pattern2;
	BitFieldHelper<uint32_t>::setField(var2, 29, 1, 0);
	BOOST_CHECK_EQUAL(var2, 0x5FFFFFEE);

	// setting a single bit
	var2 = pattern2;
	BitFieldHelper<uint32_t>::setField(var2, 31, 1, 1);
	BOOST_CHECK_EQUAL(var2, 0xFFFFFFEE);

	// checking if field overflow is clipped
	var2 = pattern1;
	BitFieldHelper<uint32_t>::setField(var2, 12, 11, 0x1F55);
	BOOST_CHECK_EQUAL(var2, 0xAF755832);
}

BOOST_AUTO_TEST_CASE(testSetField_u32_i16)
{
	uint32_t var = pattern1;
	BitFieldHelper<uint32_t>::setField(var, 16, 12, int16_t(-1234));
	BOOST_CHECK_EQUAL(var, 0xAB2E7832);
}

BOOST_AUTO_TEST_CASE(testSetFieldTpl_i32)
{
	int32_t var = 0xAAAAAAAA;
	BitFieldHelper<int32_t>::setField<28, 4>(var, -1);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFAAAAAAA));
	BitFieldHelper<int32_t>::setField<20, 8>(var, -10);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFF6AAAAA));
	BitFieldHelper<int32_t>::setField<11, 9>(var, -100);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFF6CE2AA));
	BitFieldHelper<int32_t>::setField<0, 11>(var, -200);
	BOOST_CHECK_EQUAL(var, pattern3);
}

BOOST_AUTO_TEST_CASE(testSetFieldTpl_i32_u32)
{
	int32_t var = 0xAAAAAAAA;
	BitFieldHelper<int32_t>::setField<0, 31>(var, 0xFEDCBA98);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFEDCBA98));
}

BOOST_AUTO_TEST_CASE(testSetField_i32)
{
	int32_t var = 0xAAAAAAAA;
	BitFieldHelper<int32_t>::setField(var, 28, 4, -1);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFAAAAAAA));
	BitFieldHelper<int32_t>::setField(var, 20, 8, -10);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFF6AAAAA));
	BitFieldHelper<int32_t>::setField(var, 11, 9, -100);
	BOOST_CHECK_EQUAL(var, static_cast<int32_t>(0xFF6CE2AA));
	BitFieldHelper<int32_t>::setField(var, 0, 11, -200);
	BOOST_CHECK_EQUAL(var, pattern3);
}


BOOST_AUTO_TEST_SUITE_END()
