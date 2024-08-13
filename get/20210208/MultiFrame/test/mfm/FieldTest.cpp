/*
 * $Id: FieldTest.cpp 705 2012-04-04 09:03:20Z psizun $
 * @file FieldTest.cpp
 * @created 21 févr. 2012
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
struct field_test_data
{
};

typedef test_group<field_test_data> testgroup;
typedef testgroup::object testobject;
testgroup field_testgroup("Field");

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Checks setting value of item fields.
 */
template<> template<>
void testobject::test<1>()
{
    set_test_name("Checks field-related methods");

    // Decode frame
    string path(TEMP_DIR "cobo_frame.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Set field value
    Item secondItem = coboFrame->itemAt(1);
    secondItem.field(1, 2).setValue(13u);
    Item thirdItem = coboFrame->itemAt(2);
    thirdItem.field(3, 1).setValue(255u);
    thirdItem.field(0, 2).setValue(63000u);

    // Re-encode frame
    string path2(TEMP_DIR "cobo_frame2.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    coboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks reading back field values.
 */
template<> template<>
void testobject::test<2>()
{
    set_test_name("Checks reading back field values");

    // Decode frame
    string path(TEMP_DIR "cobo_frame2.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Get field values
    ensure_equals("First field value", coboFrame->itemAt(1).field(1, 2).value<size_t>(), 13u);
    ensure_equals("Second field value", coboFrame->itemAt(2).field(3, 1).value<size_t>(), 255u);
    ensure_equals("Third field value", coboFrame->itemAt(2).field(0, 2).value<size_t>(), 63000u);
}

/**
 * Checks setting value of header fields.
 */
template<> template<>
void testobject::test<3>()
{
    set_test_name("Checks setting value of header fields");

    // Decode frame
    string path(TEMP_DIR "cobo_frame2.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    // Set CoBo 3.0 header field 'eventIdx'
    coboFrame->headerField(22,4).setValue(4000000000u);
    // Set CoBo 3.0 header field 'readOffset'
    coboFrame->headerField(28,2).setValue(148);

    // Re-encode frame
    string path2(TEMP_DIR "cobo_frame3.bin");
    ofstream ofs(path2.c_str(), ios::out | ios::binary);
    coboFrame->write(ofs);
    ofs.close();
}

/**
 * Checks getting item field by name.
 */
template<> template<>
void testobject::test<4>()
{
    set_test_name("Checks getting item field by name");

    // Decode frame
    string path(TEMP_DIR "cobo_frame3.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    FrameDictionary::instance().addFormats("mfm/sample_formats.xml");

    // Get field
    Field field = coboFrame->itemAt(1).field("");
    ensure_equals("fieldSize", field.size_B(), 4u);
}

/**
 * Checks getting header field by name.
 */
template<> template<>
void testobject::test<5>()
{
    set_test_name("Checks getting header field by name");

    // Decode frame
    string path(TEMP_DIR "cobo_frame3.bin");
    ifstream ifs(path.c_str(), ios::binary);
    std::auto_ptr< Frame > coboFrame = Frame::read(ifs);
    ifs.close();

    FrameDictionary::instance().addFormats("mfm/sample_formats.xml");

    // Get field
    Field field = coboFrame->headerField("eventTime");
    ensure_equals("fieldSize", field.size_B(), 6u);
    field = coboFrame->headerField("hitPat_0");
    ensure_equals("fieldSize", field.size_B(), 9u);
}

} // namespace tut
