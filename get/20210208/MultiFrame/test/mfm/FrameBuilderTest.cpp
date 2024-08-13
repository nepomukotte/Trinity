/*
 * $Id: FrameBuilderTest.cpp 742 2012-06-12 13:36:02Z psizun $
 * @file FrameBuilderTest.cpp
 * @created 16 avr. 2012
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
#include "mfm/FrameBuilder.h"
#include "mfm/Frame.h"
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

struct TestFrameBuilder : public FrameBuilder
{
	TestFrameBuilder() : frameCount(0) {}
	void reset() { FrameBuilder::reset(); frameCount = 0; }
	void processFrame(mfm::Frame &)
	{
		frameCount++;
	}
	size_t frameCount;
};

namespace tut
{
///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////
struct frame_builder_test_data
{
	TestFrameBuilder builder;
};

typedef test_group<frame_builder_test_data> testgroup;
typedef testgroup::object testobject;
testgroup frame_builder_testgroup("FrameBuilder");

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Checks frame reconstruction.
 */
template<> template<>
void testobject::test<1>()
{
    set_test_name("Checks frame reconstruction.");

    // Dump 5 frames to binary file
    string path(TEMP_DIR "cobo_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();
    path = string(TEMP_DIR "cobo_frames.bin");
    ofstream ofs(path.c_str(), ios::binary | ios::trunc);
    for (size_t i=0; i< 5u; ++i)
    {
    	coboFrame->write(ofs);
    }
    ofs.close();

    // Read from files by chunks of 3 Bytes
    ifs.open(path.c_str(), ios::binary);
    size_t const CHUNK_SIZE = 3;
    char buffer[CHUNK_SIZE];
    while (ifs.good())
    {
        try
        {
        	ifs.read(buffer, CHUNK_SIZE);
        }
        catch (const ifstream::failure & e)
        {
        	if (not ifs.eof()) throw;
        }
        size_t const numRead = ifs.gcount();
        builder.addDataChunk(buffer, buffer + numRead);
    }

    ensure_equals("frameCount", builder.frameCount, 5u);
}

template<> template<>
void testobject::test<2>()
{
    set_test_name("Checks frame reconstruction with big chunks.");

    // Read from files by chunks of 3 Bytes
    string path = string(TEMP_DIR "cobo_frames.bin");
    ifstream ifs(path.c_str(), ios::binary);
    size_t const CHUNK_SIZE = 390;
    char buffer[CHUNK_SIZE];
    builder.reset();
    while (ifs.good())
    {
        try
        {
        	ifs.read(buffer, CHUNK_SIZE);
        }
        catch (const ifstream::failure & e)
        {
        	if (not ifs.eof()) throw;
        }
        size_t const numRead = ifs.gcount();
        builder.addDataChunk(buffer, buffer + numRead);
    }

    ensure_equals("frameCount", builder.frameCount, 5u);
}

} // namespace tut
