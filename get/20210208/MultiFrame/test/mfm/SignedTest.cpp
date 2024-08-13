/*
 * $Id: SignedTest.cpp 705 2012-04-04 09:03:20Z psizun $
 * @file SignedTest.cpp
 * @created 22 mars 2012
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
#include "mfm/FrameFormat.h"
#include "mfm/Exception.h"
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
struct signed_test_data
{
};

typedef test_group<signed_test_data> testgroup;
typedef testgroup::object testobject;
testgroup signed_testgroup("Signed");

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Checks encoding a signed variable into a larger field.
 */
template<> template<>
void testobject::test<1>()
{
    set_test_name("Checks encoding a signed variable into a larger field");

	try
	{
		// Create frame
	    FrameFormat format("mfm/sample_formats.xml", "CoBo", 0);
		std::auto_ptr<Frame> frame = format.createFrame();

		// Set negative value into field larger than the value container
		int8_t signed_value = -5;
		frame->headerField("readOffset").setValue(signed_value);
		ensure_equals("Value of negative signed variable", frame->headerField("readOffset").value<int16_t>(), int16_t(signed_value));
	    std::string readOffset;
	    boost::to_string(frame->headerField("readOffset").value<DynamicBitset>(), readOffset);
		ensure_equals("Bits of negative signed variable", readOffset, "1111111111111011");

		// Set unsigned value into field larger than the value container
		uint8_t unsigned_value = 0xFB;
		frame->headerField("readOffset").setValue(unsigned_value);
		ensure_equals("Value of unsigned variable", frame->headerField("readOffset").value<uint16_t>(), uint16_t(unsigned_value));
	    boost::to_string(frame->headerField("readOffset").value<DynamicBitset>(), readOffset);
		ensure_equals("Bits of unsigned variable", readOffset, "0000000011111011");

	    // Encode frame
	    string path(TEMP_DIR "signed_frame1.bin");
	    ofstream ofs(path.c_str(), ios::out | ios::binary);
	    frame->write(ofs);
	    ofs.close();
	}
	catch (mfm::Exception & e)
	{
		fail("Test of signed variable failed");
	}
}

/**
 * Checks encoding a signed variable into a same-length field.
 */
template<> template<>
void testobject::test<2>()
{
    set_test_name("Checks encoding a signed variable into a same-length field");

	try
	{
		// Create frame
	    FrameFormat format("mfm/sample_formats.xml", "CoBo", 0);
		std::auto_ptr<Frame> frame = format.createFrame();

		// Set negative value
		int16_t signed_value = -5;
		frame->headerField("readOffset").setValue(signed_value);
		ensure_equals("Value of negative signed variable", frame->headerField("readOffset").value<int16_t>(), int16_t(signed_value));
	    std::string readOffset;
	    boost::to_string(frame->headerField("readOffset").value<DynamicBitset>(), readOffset);
		ensure_equals("Bits of negative signed variable", readOffset, "1111111111111011");

		// Set unsigned value
		uint16_t unsigned_value = 0xFB;
		frame->headerField("readOffset").setValue(unsigned_value);
		ensure_equals("Value of unsigned variable", frame->headerField("readOffset").value<uint16_t>(), uint16_t(unsigned_value));
	    boost::to_string(frame->headerField("readOffset").value<DynamicBitset>(), readOffset);
		ensure_equals("Bits of unsigned variable", readOffset, "0000000011111011");

	    // Encode frame
	    string path(TEMP_DIR "signed_frame2.bin");
	    ofstream ofs(path.c_str(), ios::out | ios::binary);
	    frame->write(ofs);
	    ofs.close();
	}
	catch (mfm::Exception & e)
	{
		fail("Test of signed variable failed");
	}
}

/**
 * Checks encoding a signed variable into a larger bit field.
 */
template<> template<>
void testobject::test<3>()
{
    set_test_name("Checks encoding a signed variable into a larger bit field");

	try
	{
		// Create frame
	    FrameFormat format("mfm/sample_formats.xml", "CoBo", 0);
		std::auto_ptr<Frame> frame = format.createFrame();
/*
		// Set negative value
		int8_t signed_value = -5;
		size_t pos = 3, width = 10;
		frame->headerField("eventIdx").setValue(0);
		frame->headerField("eventIdx").setBitField<int32_t>(pos, width, signed_value);
	    // Encode frame
	    string path(TEMP_DIR "signed_frame.bin");
	    ofstream ofs(path.c_str(), ios::out | ios::binary);
	    frame->write(ofs);
	    ofs.close();
	    std::string eventIdx;
	    boost::to_string(frame->headerField("eventIdx").bitField<DynamicBitset>(pos, width), eventIdx);
		ensure_equals("Bits of negative signed variable", eventIdx, "1111111011");
		*/

		// Set negative value
		int8_t signed_value = -5;
		size_t pos = 3, width = 10;
		frame->headerField("eventIdx").setValue(0);
		DynamicBitset bitset(width, signed_value);
		frame->headerField("eventIdx").bitField(pos, width).setValue(bitset);
	    // Encode frame
	    string path(TEMP_DIR "signed_frame3.bin");
	    ofstream ofs(path.c_str(), ios::out | ios::binary);
	    frame->write(ofs);
	    ofs.close();
	    std::string eventIdx;
	    boost::to_string(frame->headerField("eventIdx").bitField(pos, width).value<DynamicBitset>(), eventIdx);
		ensure_equals("Bits of negative signed variable", eventIdx, "1111111011");
/*
		// Set unsigned value
		uint8_t unsigned_value = 0xFB;
		frame->headerField("eventIdx").setValue(0);
		frame->headerField("eventIdx").setBitField(pos, width, unsigned_value);
		ensure_equals("Value of unsigned variable", frame->headerField("eventIdx").bitField<uint16_t>(pos, width), uint16_t(unsigned_value));
	    boost::to_string(frame->headerField("eventIdx").bitField<DynamicBitset>(pos, width), eventIdx);
		ensure_equals("Bits of unsigned variable", eventIdx, "0011111011");

	    // Encode frame
	    string path(TEMP_DIR "signed_frame.bin");
	    ofstream ofs(path.c_str(), ios::out | ios::binary);
	    frame->write(ofs);
	    ofs.close();*/
	}
	catch (mfm::Exception & e)
	{
		fail("Test of signed variable failed");
	}
}

} // namespace tut


