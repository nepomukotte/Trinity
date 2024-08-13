/**
 * @file BitSetTest.cpp
 * @created Oct 8, 2013
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the Utilities software project.
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

#include "utl/BitSet.hpp"
#include <boost/test/auto_unit_test.hpp>
#include <boost/cstdint.hpp>
#include <vector>
#include <utl/Logging.h>
#include <iostream>

using namespace utl;
using namespace std;

struct BitSetFixture
{
	typedef  int16_t Atom1;
	typedef  int32_t Atom2;
	typedef vector<Atom1> Buffer1;
	typedef vector<Atom2> Buffer2;

	typedef boost::make_signed  <Atom1>::type AtomS1;
	typedef boost::make_unsigned<Atom1>::type AtomU1;
	typedef boost::make_signed  <Atom2>::type AtomS2;
	typedef boost::make_unsigned<Atom2>::type AtomU2;

	static const unsigned BIT_CNT_1 =  17u;
	static const unsigned BIT_CNT_2 = 125u;
	static const unsigned BIT_MOD_1 = BIT_CNT_1 % (sizeof(Atom1)*CHAR_BIT);
	static const unsigned BIT_MOD_2 = BIT_CNT_2 % (sizeof(Atom2)*CHAR_BIT);
	static const unsigned BIT_REM_1 = ((sizeof(Atom1)*CHAR_BIT) - BIT_MOD_1) % (sizeof(Atom1)*CHAR_BIT);
	static const unsigned BIT_REM_2 = ((sizeof(Atom2)*CHAR_BIT) - BIT_MOD_2) % (sizeof(Atom2)*CHAR_BIT);

	typedef BitSet< BIT_CNT_1, AtomS1 > BitSetS1;
	typedef BitSet< BIT_CNT_1, AtomU1 > BitSetU1;
	typedef BitSet< BIT_CNT_2, AtomS2 > BitSetS2;
	typedef BitSet< BIT_CNT_2, AtomU2 > BitSetU2;

	BitSetFixture();
	~BitSetFixture() {}

	BitSetS1 bitSetS1;
	BitSetU1 bitSetU1;
	BitSetS2 bitSetS2;
	BitSetU2 bitSetU2;

	Buffer1 buf1;
	Buffer2 buf2;

	template <typename W>
	static void binPrint(ostream& stream, const vector<W>& v);

#ifdef WORDS_BIG_ENDIAN
	Atom1& buf1_lsa() { return buf1.back (); }
	Atom1& buf1_msa() { return buf1.front(); }
	Atom2& buf2_lsa() { return buf2.back (); }
	Atom2& buf2_msa() { return buf2.front(); }
#else
	Atom1& buf1_lsa() { return buf1.front(); }
	Atom1& buf1_msa() { return buf1.back (); }
	Atom2& buf2_lsa() { return buf2.front(); }
	Atom2& buf2_msa() { return buf2.back (); }
#endif

	Atom1* buf1_ptr() { return &buf1.front(); }
	Atom2* buf2_ptr() { return &buf2.front(); }
};

BitSetFixture::BitSetFixture()
{
	buf1.clear();
	buf1.resize(BitSetS1::atomCount());
	buf2.clear();
	buf2.resize(BitSetU2::atomCount());
}

template <typename W>
void BitSetFixture::binPrint(ostream& stream, const vector<W>& v)
{
	stream << "(v)" << hex << setfill('0');

	size_t iBit = v.size() * sizeof(W) * CHAR_BIT - 1;
#ifdef WORDS_BIG_ENDIAN
	for (typename vector<W>::const_iterator it = v.begin(); it != v.end(); ++it)
#else
	for (typename vector<W>::const_reverse_iterator it = v.rbegin(); it != v.rend(); ++it)
#endif
	{
		for (size_t j=sizeof(W) * CHAR_BIT; j > 0 ; --j)
		{
			stream << (((*it) >> (j-1)) & 1);
			if (iBit != 0)
			{
				if (iBit % (sizeof(W) * CHAR_BIT) == 0) stream << '|';
				else if(iBit % 8 == 0) stream << '\'';
				else if(iBit % 4 == 0) stream << ' ';
			}
			--iBit;
		}
	}
	stream << " @ " << (v.size() > 0 ? ((const void*)&v[0]) : ((void*) 0));
	stream << dec << endl;
}

BOOST_FIXTURE_TEST_SUITE(utl_BitSet, BitSetFixture)

BOOST_AUTO_TEST_CASE(reset_BitSet)
{
	bitSetS1.resetBuffer(buf1_ptr()).fill(static_cast<Atom1>(0xDEAD));
	bitSetS2.resetBuffer(buf2_ptr()).fill(static_cast<Atom2>(0xDEADFACE));
	BOOST_CHECK_EQUAL(buf1_msa(), static_cast<Atom1>(0xFFFF));
	BOOST_CHECK_EQUAL(buf1_lsa(), static_cast<Atom1>(0xDEAD));
	BOOST_CHECK_EQUAL(buf2_msa(), static_cast<Atom2>(0xFEADFACE));
	BOOST_CHECK_EQUAL(buf2_lsa(), static_cast<Atom2>(0xDEADFACE));

	bitSetU1.resetBuffer(buf1_ptr()).fill(static_cast<Atom1>(0xDEAD));
	bitSetU2.resetBuffer(buf2_ptr()).fill(static_cast<Atom2>(0xDEADFACE));
	BOOST_CHECK_EQUAL(buf1_msa(), static_cast<Atom1>(1));
	BOOST_CHECK_EQUAL(buf1_lsa(), static_cast<Atom1>(0xDEAD));
	BOOST_CHECK_EQUAL(buf2_msa(), static_cast<AtomU2>(0x1EADFACE));
	BOOST_CHECK_EQUAL(buf2_lsa(), static_cast<AtomU2>(0xDEADFACE));
}

BOOST_AUTO_TEST_CASE(fill_BitSet)
{
	bitSetS1.resetBuffer(buf1_ptr());
	bitSetS1.fill(-1);
	bitSetU2.resetBuffer(buf2_ptr());
	bitSetU2.fill(-1);
	BOOST_CHECK_EQUAL(buf1_lsa(), -1);
	BOOST_CHECK_EQUAL(buf1_msa(), -1);
	BOOST_CHECK_EQUAL(buf2_lsa(), -1);
	BOOST_CHECK_EQUAL(buf2_msa(), 0x1fffffff);
}

BOOST_AUTO_TEST_CASE(setValue_BitSet)
{
	bitSetS1.resetBuffer(buf1_ptr());
	bitSetS1.setValue(static_cast<AtomS1>(0xDEADFACE));
	BOOST_CHECK_EQUAL(buf1_msa(), -1);
	BOOST_CHECK_EQUAL(buf1_lsa(), static_cast<AtomS1>(0xDEADFACE));

	bitSetU1.resetBuffer(buf1_ptr());
	bitSetU1.setValue(static_cast<AtomU1>(0xDEADFACE));
	BOOST_CHECK_EQUAL(buf1_msa(), 0);
	BOOST_CHECK_EQUAL(buf1_lsa(), static_cast<AtomS1>(0xDEADFACE));

	bitSetS2.resetBuffer(buf2_ptr());
	bitSetS2 = 0xDEADFACE;
	BOOST_CHECK_EQUAL(buf2_msa(), 0);
	BOOST_CHECK_EQUAL(buf2_lsa(), static_cast<AtomS2>(0xDEADFACE));

	bitSetU2.resetBuffer(buf2_ptr());
	bitSetU2 = 0xDEADFACE;
	BOOST_CHECK_EQUAL(buf2_msa(), 0);
	BOOST_CHECK_EQUAL(buf2_lsa(), static_cast<AtomU2>(0xDEADFACE));

	bitSetS1 = bitSetS2;
	BOOST_CHECK_EQUAL(buf1_msa(), -1);
	BOOST_CHECK_EQUAL(buf1_lsa(), static_cast<AtomS1>(0xDEADFACE));

	bitSetU2 = bitSetS1;
	BOOST_CHECK_EQUAL(buf2_msa(), -1u >> BIT_REM_2);
	BOOST_CHECK_EQUAL(buf2_lsa(), static_cast<AtomU2>(static_cast<AtomS1>(0xDEADFACE)));
}

BOOST_AUTO_TEST_CASE(getValue_BitSet)
{
	bitSetS1.resetBuffer(buf1_ptr()).fill(static_cast<AtomS1>(0xDEAD));
	int8_t   value_sc = bitSetS1;
	uint8_t  value_uc = bitSetS1;
	int16_t  value_ss = bitSetS1;
	uint16_t value_us = bitSetS1;
	int32_t  value_si = bitSetS1;
	uint32_t value_ui = bitSetS1;
	int64_t  value_sl = bitSetS1;
	uint64_t value_ul = bitSetS1;
	BOOST_CHECK_EQUAL(value_sc, static_cast<int8_t>  (0xDEAD));
	BOOST_CHECK_EQUAL(value_uc, static_cast<uint8_t> (0xDEAD));
	BOOST_CHECK_EQUAL(value_ss, static_cast<int16_t> (0xDEAD));
	BOOST_CHECK_EQUAL(value_us, static_cast<uint16_t>(0xDEAD));
	BOOST_CHECK_EQUAL(value_si, static_cast<int32_t> (0xFFFFDEAD));
	BOOST_CHECK_EQUAL(value_ui, static_cast<uint32_t>(0xFFFFDEAD));
	BOOST_CHECK_EQUAL(value_sl, static_cast<int64_t> (0xFFFFFFFFFFFFDEADll));
	BOOST_CHECK_EQUAL(value_ul, static_cast<uint64_t>(0xFFFFFFFFFFFFDEADllu));

	bitSetU2.resetBuffer(buf2_ptr()).fill(static_cast<AtomU2>(0xDEADFACE));
	value_sc = bitSetU2;
	value_uc = bitSetU2;
	value_ss = bitSetU2;
	value_us = bitSetU2;
	value_si = bitSetU2;
	value_ui = bitSetU2;
	value_sl = bitSetU2;
	value_ul = bitSetU2;
	BOOST_CHECK_EQUAL(value_sc, static_cast<int8_t>  (0xDEADFACE));
	BOOST_CHECK_EQUAL(value_uc, static_cast<uint8_t> (0xDEADFACE));
	BOOST_CHECK_EQUAL(value_ss, static_cast<int16_t> (0xDEADFACE));
	BOOST_CHECK_EQUAL(value_us, static_cast<uint16_t>(0xDEADFACE));

	BOOST_CHECK_EQUAL(value_si, static_cast<int32_t> (0xDEADFACE));
	BOOST_CHECK_EQUAL(value_ui, static_cast<uint32_t>(0xDEADFACE));
	BOOST_CHECK_EQUAL(value_sl, static_cast<int64_t> (0xDEADFACEDEADFACEll));
	BOOST_CHECK_EQUAL(value_ul, static_cast<uint64_t>(0xDEADFACEDEADFACEllu));
}

BOOST_AUTO_TEST_CASE(bitwise_BitSet)
{
	bitSetS1.resetBuffer(buf1_ptr()).fill(static_cast<AtomS1>(0xDEAD));
	bitSetU2.resetBuffer(buf2_ptr()).fill(static_cast<AtomU2>(0xDEADFACE));

	uint64_t result = bitSetS1 & bitSetU2;
	BOOST_CHECK_EQUAL(result, 0xffffffffffffda8cllu);

	result = bitSetS1 | bitSetU2;
	BOOST_CHECK_EQUAL(result, 0xfffffffffffffeefllu);

	result = bitSetS1 ^ bitSetU2;
	BOOST_CHECK_EQUAL(result, 0x0000000000002463llu);

	result = bitSetS1 >> 3ul;
	BOOST_CHECK_EQUAL(result, 0xfffffffffffffbd5llu);

	result = bitSetS1 >> 16ul;
	BOOST_CHECK_EQUAL(result, -1);

	result = bitSetS1 << 10ul;
	BOOST_CHECK_EQUAL(result, 0x000000000000b400llu);

	result = bitSetS1 << 32ul;
	BOOST_CHECK_EQUAL(result, 0x0000000000000000llu);

	result = bitSetS1 >> 32ul;
	BOOST_CHECK_EQUAL(result, -1);
}

BOOST_AUTO_TEST_SUITE_END()
