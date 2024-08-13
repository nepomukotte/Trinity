/*
 * $Id: FrameTest.cpp 763 2012-07-02 09:10:14Z psizun $
 * @file FrameTest.cpp
 * @created 13 févr. 2012
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
#include "mfm/StandardHeader.h"
#include "mfm/Serializer.h"
#include "mfm/Exception.h"
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
struct frame_test_data
{
    Frame frame;
};

typedef test_group<frame_test_data> testgroup;
typedef testgroup::object testobject;
testgroup frame_testgroup("Frame");

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Checks default constructor.
 */
template<> template<>
void testobject::test<1>()
{
    set_test_name("Checks default constructor");

    ensure_equals("isLayeredFrame", frame.header().isLayeredFrame(), true);
    ensure_equals("isBigEndian", frame.header().isBigEndian(), true);
}

/**
 * Checks serialization of CoBo 3.0 format frame (with no items).
 */
template<> template<>
void testobject::test<2>()
{
    set_test_name("Checks serialization of CoBo 3.0 format frame (with no items)");

    StandardHeader coboHeader(utl::BigEndian);
    coboHeader.setMetaType(BASIC);
    coboHeader.set_blkSize_B(64u);
    coboHeader.set_headerSize_B(128u);
    coboHeader.set_dataSize_B(0u); // frame size with no items
    coboHeader.set_dataSource(0u);
    coboHeader.set_frameType(1u);
    coboHeader.set_revision(1u);
    coboHeader.set_itemSize_B(4u);

    Frame coboFrame(coboHeader);

    string path(TEMP_DIR "cobo_empty_frame.bin");
    ofstream ofs(path.c_str(), ios::out | ios::binary);
    coboFrame.write(ofs);
    ofs.close();

    // Check size of binary file
    ifstream ifs(path.c_str(), ios::in | ios::binary);
    size_t begin = ifs.tellg();
    ifs.seekg (0, ios::end);
    size_t end = ifs.tellg();
    size_t size = end - begin;
    ifs.close();
    ensure_equals("Size of Cobo 3.0 header", size, 128u);
}

/**
 * Checks de-serialization of CoBo 3.0 format frame (with no items).
 */
template<> template<>
void testobject::test<3>()
{
    set_test_name("Checks de-serialization of CoBo 3.0 format frame (with no items)");

    string path(TEMP_DIR "cobo_empty_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > emptyCoboFrame = Frame::read(ifs);
    ifs.close();

    ensure_equals("isBigEndian", emptyCoboFrame->header().isBigEndian(), true);
    ensure_equals("endianness", emptyCoboFrame->header().endianness(), utl::BigEndian);
    ensure_equals("isBlobFrame", emptyCoboFrame->header().isBlobFrame(), false);
    ensure_equals("blkSize_B", emptyCoboFrame->header().blkSize_B(), 64u);
    ensure_equals("frameSize_blk", emptyCoboFrame->header().frameSize_blk(), 2u);
    ensure_equals("frameSize_B", emptyCoboFrame->header().frameSize_B(), 128u);
    ensure_equals("headerSize_B", emptyCoboFrame->header().headerSize_B(), 128u);
    ensure_equals("dataSize_B", emptyCoboFrame->header().dataSize_B(), 0u);
    ensure_equals("dataSource", emptyCoboFrame->header().dataSource(), 0u);
    ensure_equals("frameType", emptyCoboFrame->header().frameType(), 1u);
    ensure_equals("revision", emptyCoboFrame->header().revision(), 1u);
    ensure_equals("metaType", emptyCoboFrame->header().metaType(), BASIC);
    ensure_equals("itemSize_B", emptyCoboFrame->header().itemSize_B(), 4u);
    ensure_equals("itemCount", emptyCoboFrame->header().itemCount(), 0u);

    emptyCoboFrame.reset();
}

/**
 * Checks Frame::addItem.
 */
template<> template<>
void testobject::test<4>()
{
    set_test_name("Checks Frame::addItem");

    string path(TEMP_DIR "cobo_empty_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > emptyCoboFrame = Frame::read(ifs);
    ifs.close();

    emptyCoboFrame->addItem();
    emptyCoboFrame->addItem();
    emptyCoboFrame->addItem();

    ensure_equals("itemCount", emptyCoboFrame->header().itemCount(), 3u);
    ensure_equals("headerSize_B", emptyCoboFrame->header().headerSize_B(), 128u);
    ensure_equals("dataSize_B", emptyCoboFrame->header().dataSize_B(), 64u);
    ensure_equals("frameSize_B", emptyCoboFrame->header().frameSize_B(), 192u);

    string path2(TEMP_DIR "cobo_frame.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    emptyCoboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks de-serialization of CoBo 3.0 format frame (with items).
 */
template<> template<>
void testobject::test<5>()
{
    set_test_name("Checks de-serialization of CoBo 3.0 format frame (with items)");

    // Check size of binary file
    string path(TEMP_DIR "cobo_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    size_t begin = ifs.tellg();
    ifs.seekg (0, ios::end);
    size_t end = ifs.tellg();
    size_t size = end - begin;
    ensure_equals("Size of Cobo 3.0 frame", size, 192u);

    // Decode frame
    ifs.seekg (0, ios::beg);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    ensure_equals("isBigEndian", coboFrame->header().isBigEndian(), true);
    ensure_equals("endianness", coboFrame->header().endianness(), utl::BigEndian);
    ensure_equals("isBlobFrame", coboFrame->header().isBlobFrame(), false);
    ensure_equals("blkSize_B", coboFrame->header().blkSize_B(), 64u);
    ensure_equals("frameSize_blk", coboFrame->header().frameSize_blk(), 3u);
    ensure_equals("frameSize_B", coboFrame->header().frameSize_B(), 192u);
    ensure_equals("headerSize_B", coboFrame->header().headerSize_B(), 128u);
    ensure_equals("dataSize_B", coboFrame->header().dataSize_B(), 64u);
    ensure_equals("dataSource", coboFrame->header().dataSource(), 0u);
    ensure_equals("frameType", coboFrame->header().frameType(), 1u);
    ensure_equals("revision", coboFrame->header().revision(), 1u);
    ensure_equals("metaType", coboFrame->header().metaType(), BASIC);
    ensure_equals("itemSize_B", coboFrame->header().itemSize_B(), 4u);
    ensure_equals("itemCount", coboFrame->header().itemCount(), 3u);
}
/**
 * Checks access at non-existing item raises exception.
 */
template<> template<>
void testobject::test<6>()
{
    set_test_name("Checks access at non-existing item raises exception");

    // Decode frame
    string path(TEMP_DIR "cobo_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Get field values
    try
    {
    	coboFrame->itemAt(5);
    	fail("exception expected");
    }
    catch (std::exception & e)
    {
    	;
    }
}

/**
 * Checks reading an incomplete frame raises an exception.
 */
template<> template<>
void testobject::test<7>()
{
    set_test_name("Checks reading an incomplete frame raises an exception");

    try
    {
		// Decode frame
		string path(TEMP_DIR "cobo_standard_header.bin");
		ifstream ifs(path.c_str(), ios::binary);
		std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
		ifs.close();
    	fail("exception expected");
    }
    catch (std::exception & e)
    {
    	;
    }
}

/**
 * Checks searching for a field by name.
 */
template<> template<>
void testobject::test<8>()
{
    set_test_name("Checks searching for a field by name");

    try
    {
        // Decode frame
        string path(TEMP_DIR "cobo_frame.bin");
        ifstream ifs(path.c_str(), ios::binary);
        std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
        ifs.close();

        // Search for field ''
        size_t fieldPos, fieldSize;
        try
        {
        	coboFrame->findItemField("", 0, fieldPos, fieldSize);
        	ensure_equals("findItemField.size", fieldSize, 4u);
        	ensure_equals("findItemField.pos", fieldPos, 128u);
        }
        catch (mfm::FieldNotFound & e)
        {
        	fail("findItemField failed");
        }
    }
    catch (std::exception & e)
    {
    	;
    }
}

/**
 * Checks Frame constructor with Serializer argument.
 */
template<> template<>
void testobject::test<9>()
{
    set_test_name("Checks Frame constructor with Serializer argument");

        // Read frame data
        string path(TEMP_DIR "cobo_frame.bin");
        ifstream ifs(path.c_str(), ios::binary);
        // Get size
        size_t begin = ifs.tellg();
        ifs.seekg (0, ios::end);
        size_t end = ifs.tellg();
        size_t byteCount = end - begin;
        // Initialize buffer
        mfm::Serializer dataSerializer(byteCount);
        // Read data into buffer
        ifs.seekg (0, ios::beg);
        dataSerializer.read(ifs, byteCount);
        ifs.close();

        // Create frame
        Frame frame(dataSerializer);

        ensure_equals("itemCount", frame.header().itemCount(), 3u);
}

/**
 * Checks de-serialization from a char buffer.
 */
template<> template<>
void testobject::test<10>()
{
    set_test_name("Checks de-serialization from a char buffer");

    // Check size of binary file
    string path(TEMP_DIR "cobo_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    size_t begin = ifs.tellg();
    ifs.seekg (0, ios::end);
    size_t end = ifs.tellg();
    const size_t size = end - begin;
    ensure_equals("Size of Cobo 3.0 frame", size, 192u);
    // Read file contents into char buffer
    std::vector< char > buffer(size);
    ifs.seekg(0, ios::beg);
    ifs.read(&buffer[0], size);
    ifs.close();

    // Decode frame
    const char* beg = &buffer[0];
    std::auto_ptr< Frame > coboFrame = Frame::read(beg, beg+end);

    ensure_equals("isBigEndian", coboFrame->header().isBigEndian(), true);
    ensure_equals("endianness", coboFrame->header().endianness(), utl::BigEndian);
    ensure_equals("isBlobFrame", coboFrame->header().isBlobFrame(), false);
    ensure_equals("blkSize_B", coboFrame->header().blkSize_B(), 64u);
    ensure_equals("frameSize_blk", coboFrame->header().frameSize_blk(), 3u);
    ensure_equals("frameSize_B", coboFrame->header().frameSize_B(), 192u);
    ensure_equals("headerSize_B", coboFrame->header().headerSize_B(), 128u);
    ensure_equals("dataSize_B", coboFrame->header().dataSize_B(), 64u);
    ensure_equals("dataSource", coboFrame->header().dataSource(), 0u);
    ensure_equals("frameType", coboFrame->header().frameType(), 1u);
    ensure_equals("revision", coboFrame->header().revision(), 1u);
    ensure_equals("metaType", coboFrame->header().metaType(), BASIC);
    ensure_equals("itemSize_B", coboFrame->header().itemSize_B(), 4u);
    ensure_equals("itemCount", coboFrame->header().itemCount(), 3u);
}

} // namespace tut


