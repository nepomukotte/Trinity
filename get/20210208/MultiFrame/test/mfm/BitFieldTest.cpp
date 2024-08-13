/*
 * $Id: BitFieldTest.cpp 740 2012-06-05 14:59:08Z psizun $
 * @file BitFieldTest.cpp
 * @created 28 mars 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 * -----------------------------------------------------------------------------
 */

#include "mfm/Frame.h"
#include "mfm/FrameDictionary.h"
#include "mfm/BitField.h"
using namespace mfm;
#include "tut.h"
#include <string>
#include <fstream>
using namespace std;

#ifndef VXWORKS
#define TEMP_DIR "/tmp/"
#else
#define TEMP_DIR "tmp/"
#endif

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////
struct bitfield_test_data
{
};

typedef test_group<bitfield_test_data> testgroup;
typedef testgroup::object testobject;
testgroup bitfield_testgroup("BitField");

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Checks setting value of item bit-fields.
 */
template<> template<>
void testobject::test<1>()
{
    set_test_name("Checks setting value of item bit-fields");

    // Decode frame
    string path(TEMP_DIR "cobo_frame3.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    for (size_t i=0; i < coboFrame->itemCount(); ++i)
    {
    	// Set CoBo 3.0 item bit field 'agetIdx'
    	coboFrame->itemAt(i).field(0, 4).bitField(30, 2).setValue(3);
    	// Set CoBo 3.0 item bit field 'chanIdx'
    	coboFrame->itemAt(i).field(0, 4).bitField(23, 7).setValue(67);
    	// Set CoBo 3.0 item bit field 'buckIdx'
    	coboFrame->itemAt(i).field(0, 4).bitField(14, 9).setValue(511);
    	// Set CoBo 3.0 item bit field 'UNUSED'
    	coboFrame->itemAt(i).field(0, 4).bitField(12, 2).setValue(0);
    	// Set CoBo 3.0 item bit field 'sample'
    	coboFrame->itemAt(i).field(0, 4).bitField(0, 12).setValue(4095);
    }

    // Re-encode frame
    string path2(TEMP_DIR "cobo_frame4.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    coboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks getting value of item bit-fields.
 */
template<> template<>
void testobject::test<2>()
{
    set_test_name("Checks getting value of item bit-fields");

    // Decode frame
    string path(TEMP_DIR "cobo_frame4.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    ensure_equals("sample bit-field", coboFrame->itemAt(1).field(0, 4).bitField(0, 12).value<uint32_t>(), 4095u);
    ensure_equals("UNUSED bit-field", coboFrame->itemAt(1).field(0, 4).bitField(12, 2).value<uint32_t>(), 0u);
    ensure_equals("buckIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(14, 9).value<uint32_t>(), 511u);
    ensure_equals("chanIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(23, 7).value<uint32_t>(), 67u);
    ensure_equals("agetIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(30, 2).value<uint32_t>(), 3u);
}

/**
 * Checks getting header bit field.
 */
template<> template<>
void testobject::test<3>()
{
    set_test_name("Checks getting header bit field");

    // Decode frame
    string path(TEMP_DIR "cobo_frame4.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    for (size_t blockIndex=8; blockIndex > 0; --blockIndex)
    {
    	Field block = coboFrame->headerField(31 + blockIndex, 1);
    	for (size_t bitPos = 0; bitPos < 8; ++bitPos)
    		block.setBit(bitPos, (bitPos%2 == 0)?false:true);
    }
    Field lastBlock = coboFrame->headerField(31,1);
    lastBlock.setValue<uint8_t>(0x0A);

    std::string hitPat_0;
    boost::to_string(coboFrame->headerField(31,9).value< DynamicBitset >(), hitPat_0);
    ensure_equals("hitPat_0", hitPat_0, std::string("000010101010101010101010101010101010101010101010101010101010101010101010"));

    // Re-encode frame
    string path2(TEMP_DIR "cobo_frame5.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    coboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks setting header bit field.
 */
template<> template<>
void testobject::test<4>()
{
    set_test_name("Checks setting header bit field");

    // Decode frame
    string path(TEMP_DIR "cobo_frame4.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Set 'hitPat_0' header bit pattern
    DynamicBitset pattern(std::string("000010101010101010101010101010101010101010101010101010101010101010101010"), 0, 72);
    coboFrame->headerField(31,9).setValue(pattern);

    ensure_equals("hitPat_0", coboFrame->headerField(31,9).value< DynamicBitset >(), pattern);

    // Re-encode frame
    string path2(TEMP_DIR "cobo_frame6.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    coboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks getting value item bit field as dynamic bitset.
 */
template<> template<>
void testobject::test<5>()
{
    set_test_name("Checks getting value item bit field as dynamic bitset");

    // Decode frame
    string path(TEMP_DIR "cobo_frame4.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Get item bitfields
    ensure_equals("sample bit-field", coboFrame->itemAt(1).field(0, 4).bitField(0, 12).value<DynamicBitset>().to_ulong(), 4095u);
    ensure_equals("UNUSED bit-field", coboFrame->itemAt(1).field(0, 4).bitField(12, 2).value<DynamicBitset>().to_ulong(), 0u);
    ensure_equals("buckIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(14, 9).value<DynamicBitset>().to_ulong(), 511u);
    ensure_equals("chanIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(23, 7).value<DynamicBitset>().to_ulong(), 67u);
    ensure_equals("agetIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(30, 2).value<DynamicBitset>().to_ulong(), 3u);
}

/**
 * Checks setting item bit field.
 */
template<> template<>
void testobject::test<6>()
{
    set_test_name("Checks setting item bit field");

    // Decode frame
    string path(TEMP_DIR "cobo_frame4.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Modify 'buckIdx'
    coboFrame->itemAt(1).field(0, 4).bitField(14, 9).setValue(DynamicBitset(std::string("101010101")));

    // Get item bitfields
    ensure_equals("sample bit-field", coboFrame->itemAt(1).field(0, 4).bitField(0, 12).value<uint32_t>(), 4095u);
    ensure_equals("UNUSED bit-field", coboFrame->itemAt(1).field(0, 4).bitField(12, 2).value<uint32_t>(), 0u);
    ensure_equals("buckIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(14, 9).value<uint32_t>(), 341u);
    ensure_equals("chanIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(23, 7).value<uint32_t>(), 67u);
    ensure_equals("agetIdx bit-field", coboFrame->itemAt(1).field(0, 4).bitField(30, 2).value<uint32_t>(), 3u);
}

/**
 * Checks setting value of item bit field by name.
 */
template<> template<>
void testobject::test<7>()
{
    set_test_name("Checks setting value of item bit field by name");

    // Decode frame
    string path(TEMP_DIR "cobo_frame5.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    FrameDictionary::instance().addFormats("mfm/sample_formats.xml");

	try
	{
		for (size_t i = 0; i < coboFrame->itemCount(); ++i)
		{
			// Set CoBo 3.0 item bit field 'agetIdx'
			DynamicBitset agetIdx(2, 3);
			coboFrame->itemAt(i).field("").bitField("agetIdx").setValue(agetIdx);
			// Set CoBo 3.0 item bit field 'chanIdx'
			DynamicBitset chanIdx(7, 0x7F);
			coboFrame->itemAt(i).field(0, 4).bitField("chanIdx").setValue(chanIdx);
			// Set CoBo 3.0 item bit field 'buckIdx'
			DynamicBitset buckIdx(9, 0xAA);
			coboFrame->itemAt(i).field("").bitField("buckIdx").setValue(buckIdx);
			// Set CoBo 3.0 item bit field 'sample'
			DynamicBitset sample(12, 0xFFF);
			coboFrame->itemAt(i).field(0, 4).bitField("sample").setValue(sample);
		}
    }
    catch (mfm::BitFieldNotFound & e)
    {
    	fail("Could not find bit field by name");
    }

    // Re-encode frame
    string path2(TEMP_DIR "cobo_frame7.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    coboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks getting value of item bit field by name.
 */
template<> template<>
void testobject::test<8>()
{
    set_test_name("Checks getting value of item bit field by name");

    // Decode frame
    string path(TEMP_DIR "cobo_frame7.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    FrameDictionary::instance().addFormats("mfm/sample_formats.xml");

    for (size_t i=0; i < coboFrame->itemCount(); ++i)
    {
        ensure_equals(coboFrame->itemAt(i).field("").value<uint32_t>(), 0xFFAA8FFF);
    	// Get CoBo 3.0 item bit field 'agetIdx'
    	ensure_equals(coboFrame->itemAt(i).field("").bitField("agetIdx").value<uint32_t>(), 3u);
    	// Get CoBo 3.0 item bit field 'chanIdx'
    	ensure_equals(coboFrame->itemAt(i).field(0, 4).bitField("chanIdx").value<uint32_t>(), 0x7Fu);
    	// Get CoBo 3.0 item bit field 'buckIdx'
    	ensure_equals(coboFrame->itemAt(i).field(0, 4).bitField("buckIdx").value<uint32_t>(), 0xAAu);
    	// Get CoBo 3.0 item bit field 'sample'
    	ensure_equals(coboFrame->itemAt(i).field(0, 4).bitField("sample").value<uint32_t>(), 0xFFFu);
    }
}

} // namespace tut


