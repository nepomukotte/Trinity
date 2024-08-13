/* =============================================================================
 * BitFieldHelperVectorTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 13 nov. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include <utl/BitFieldHelperVector.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/cstdint.hpp>

using namespace utl;

struct BitFieldHelperVectorFixture
{
	typedef std::vector<uint32_t> VecUINT32;
	typedef std::vector<int32_t> VecINT32;
	typedef std::vector<uint8_t> VecUINT8;

	VecUINT32 pattern1;
	VecUINT32 pattern2;
	VecINT32 pattern3;
	VecUINT8 pattern4;

	BitFieldHelperVectorFixture()
		: pattern1(1, 0),
		  pattern2(3, 0),
		  pattern3(2, 0),
		  pattern4(4, 0)
	{
		WordIterator<VecUINT32::value_type> p1 = WordIterator<VecUINT32::value_type>::getBegin(pattern1);
		*p1++ = 0xAF3C7832;

		WordIterator<VecUINT32::value_type> p2 = WordIterator<VecUINT32::value_type>::getBegin(pattern2);
		*p2++ = 0xA0E07CA1;
		*p2++ = 0xFC05A003;
		*p2++ = 0xAC018005;

		WordIterator<VecINT32::value_type> p3 = WordIterator<VecINT32::value_type>::getBegin(pattern3);
		*p3++ = 0x0B1E0000;
		*p3++ = 0xFF6CE58F;

		WordIterator<VecUINT8::value_type> p4 = WordIterator<VecUINT8::value_type>::getBegin(pattern4);
		*p4++ = 0x32;
		*p4++ = 0x78;
		*p4++ = 0x3C;
		*p4++ = 0xAF;
	}
};

BOOST_FIXTURE_TEST_SUITE(utl_BitFieldHelperVector, BitFieldHelperVectorFixture)

// =============================================================================
// tests cases for getField method
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testGetField1)
{
	// checking if retrieving the whole vector works
	uint32_t f1 = BitFieldHelper<VecUINT32>::getField<uint32_t>(pattern1, 0, 32);
	BOOST_CHECK_EQUAL(f1, 0xAF3C7832);
	// checking if retrieving an un-aligned 8-bit field works
	uint8_t f2 = BitFieldHelper<VecUINT32>::getField<uint8_t>(pattern1, 4, 8);
	BOOST_CHECK_EQUAL(f2, 0x83);
	// checking if retrieving an un-aligned 16-bit field works
	uint16_t f3 = BitFieldHelper<VecUINT32>::getField<uint16_t>(pattern1, 12, 12);
	BOOST_CHECK_EQUAL(f3, 0x3C7);
}

BOOST_AUTO_TEST_CASE(testGetField2)
{
	uint32_t f1 = BitFieldHelper<VecUINT32>::getField<uint32_t>(pattern2, 1, 24);
	BOOST_CHECK_EQUAL(f1, 0x703E50u);
	uint32_t f2 = BitFieldHelper<VecUINT32>::getField<uint32_t>(pattern2, 28, 20);
	BOOST_CHECK_EQUAL(f2, 0xA003Au);
	uint32_t f3 = BitFieldHelper<VecUINT32>::getField<uint32_t>(pattern2, 62, 5);
	BOOST_CHECK_EQUAL(f3, 0x17u);
	uint32_t f4 = BitFieldHelper<VecUINT32>::getField<uint32_t>(pattern2, 88, 8);
	BOOST_CHECK_EQUAL(f4, 0xACu);

#ifndef BOOST_NO_INT64_T
	uint64_t f5 = BitFieldHelper<VecUINT32>::getField<uint64_t>(pattern2, 16, 64);
	BOOST_CHECK_EQUAL(f5, static_cast<uint64_t>(0x8005FC05) << 32 | 0xA003A0E0);
#endif
}

BOOST_AUTO_TEST_CASE(testGetField3)
{
	int8_t f1 = BitFieldHelper<VecINT32>::getField<int8_t>(pattern3, 62, 4);
	BOOST_CHECK_EQUAL(f1, -1);
	int8_t f2 = BitFieldHelper<VecINT32>::getField<int8_t>(pattern3, 52, 8);
	BOOST_CHECK_EQUAL(f2, -10);
	int16_t f3 = BitFieldHelper<VecINT32>::getField<int16_t>(pattern3, 43, 9);
	BOOST_CHECK_EQUAL(f3, -100);
	int32_t f4 = BitFieldHelper<VecINT32>::getField<int32_t>(pattern3, 28, 15);
	BOOST_CHECK_EQUAL(f4, -10000);
	int16_t f5 = BitFieldHelper<VecINT32>::getField<int16_t>(pattern3, 16, 13);
	BOOST_CHECK_EQUAL(f5, 0xB1E);
}

BOOST_AUTO_TEST_CASE(testGetField4)
{
	uint32_t f1 = BitFieldHelper<VecUINT8>::getField<uint32_t>(pattern4, 0, 32);
	BOOST_CHECK_EQUAL(f1, 0xAF3C7832);
	uint8_t f2 = BitFieldHelper<VecUINT8>::getField<uint8_t>(pattern4, 4, 8);
	BOOST_CHECK_EQUAL(f2, 0x83);
	int8_t f3 = BitFieldHelper<VecUINT8>::getField<int8_t>(pattern4, 20, 8);
	BOOST_CHECK_EQUAL(f3, -13);
	uint16_t f4 = BitFieldHelper<VecUINT8>::getField<uint16_t>(pattern4, 12, 16);
	BOOST_CHECK_EQUAL(f4, 0xF3C7);
}

// =============================================================================
// tests cases for setField method
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testSetField1)
{
	VecUINT32 vec = pattern1;
	BitFieldHelper<VecUINT32>::setField(vec, 0, 32, 0x12345678);
	BOOST_CHECK_EQUAL(vec[0], 0x12345678u);

	vec = pattern1;
	BitFieldHelper<VecUINT32>::setField(vec, 8, 16, uint16_t(0xBAF0));
	BOOST_CHECK_EQUAL(vec[0], 0xAFBAF032u);
}

BOOST_AUTO_TEST_CASE(testSetField_overflow)
{
	// testing field width overflow
	VecUINT32 vec = pattern1;
	BitFieldHelper<VecUINT32>::setField(vec, 24, 16, uint16_t(0xF00D));
	BOOST_CHECK_EQUAL(vec[0], 0x0D3C7832u);

	// testing field pos overflow
	vec = pattern1;
	BitFieldHelper<VecUINT32>::setField(vec, 32, 16, uint16_t(0xF00D));
	BOOST_CHECK_EQUAL(vec[0], 0xAF3C7832u);
}

BOOST_AUTO_TEST_CASE(testSetField3)
{
	VecINT32 vec = pattern3;
	BitFieldHelper<VecINT32>::setField(vec, 16, 12, -150);
	BOOST_CHECK_EQUAL(vec[0], 0x0F6A0000);

	vec = pattern3;
	BitFieldHelper<VecINT32>::setField(vec, 24, 12, -450);
	BOOST_CHECK_EQUAL(vec[0], 0x3E1E0000);
	BOOST_CHECK_EQUAL(vec[1], static_cast<int32_t>(0xFF6CE58E));
}

BOOST_AUTO_TEST_CASE(testSetField4)
{
	VecUINT8 vec = pattern4;
	BitFieldHelper<VecUINT8>::setField(vec, 4, 8, 0xAB);
	BOOST_CHECK_EQUAL(vec[0], 0xB2);
	BOOST_CHECK_EQUAL(vec[1], 0x7A);
	BOOST_CHECK_EQUAL(vec[2], 0x3C);
	BOOST_CHECK_EQUAL(vec[3], 0xAF);

	vec = pattern4;
	BitFieldHelper<VecUINT8>::setField(vec, 11, 14, 0x1234);
	BOOST_CHECK_EQUAL(vec[0], 0x32);
	BOOST_CHECK_EQUAL(vec[1], 0xA0);
	BOOST_CHECK_EQUAL(vec[2], 0x91);
	BOOST_CHECK_EQUAL(vec[3], 0xAE);
}

// =============================================================================
// tests cases for getBit method
// -----------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(testGetBit)
{

}

BOOST_AUTO_TEST_SUITE_END()
