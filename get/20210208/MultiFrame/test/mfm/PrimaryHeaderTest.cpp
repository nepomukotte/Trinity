/* =============================================================================
 * $Id: PrimaryHeaderTest.cpp 656 2012-02-17 16:45:12Z psizun $
 * @file PrimaryHeaderTest.cpp
 * @author sizun
 * @date 15 déc. 2011
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
 * =============================================================================
 */

#include "mfm/PrimaryHeader.h"
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
struct primary_header_test_data
{
    PrimaryHeader header;
};

typedef test_group<primary_header_test_data> testgroup;
typedef testgroup::object testobject;
testgroup primary_header_testgroup("PrimaryHeader");

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

    ensure_equals("isBigEndian", header.isBigEndian(), true);
    ensure_equals("endianness", header.endianness(), utl::BigEndian);
    ensure_equals("isBlobFrame", header.isBlobFrame(), true);
    ensure_equals("blkSize_B", header.blkSize_B(), 1u);
    ensure_equals("frameSize_blk", header.frameSize_blk(), 8u);
    ensure_equals("frameSize_B", header.frameSize_B(), 8u);
    ensure_equals("dataSource", header.dataSource(), 0);
    ensure_equals("frameType", header.frameType(), 0);
    ensure_equals("revision", header.revision(), 0);
}

/**
 * Checks serialization.
 */
template<> template<>
void testobject::test<2>()
{
    set_test_name("Checks serialization");

    string path(TEMP_DIR "primary_header.bin");
    ofstream ofs(path.c_str(), ios::out | ios::binary);
    header.encode(ofs);
    ofs.close();

    // Check size of binary file
    ifstream ifs(path.c_str(), ios::in | ios::binary);
    size_t begin = ifs.tellg();
    ifs.seekg (0, ios::end);
    size_t end = ifs.tellg();
    size_t size = end - begin;
    ifs.close();
    ensure_equals("Size of primary header", size, 8u);
}

/**
 * Checks unserialization.
 */
template<> template<>
void testobject::test<3>()
{
    set_test_name("Checks unserialization");

    string path(TEMP_DIR "primary_header.bin");
    ifstream ifs(path.c_str(), ios::in | ios::binary);
    PrimaryHeader header2;
    header2.decode(ifs);
    ifs.close();

    ensure_equals("isBigEndian", header2.isBigEndian(), true);
    ensure_equals("endianness", header2.endianness(), utl::BigEndian);
    ensure_equals("isBlobFrame", header2.isBlobFrame(), true);
    ensure_equals("blkSize_B", header2.blkSize_B(), 1u);
    ensure_equals("frameSize_blk", header2.frameSize_blk(), 8u);
    ensure_equals("frameSize_B", header2.frameSize_B(), 8u);
    ensure_equals("dataSource", header2.dataSource(), 0);
    ensure_equals("frameType", header2.frameType(), 0);
    ensure_equals("revision", header2.revision(), 0);
}

/**
 * Checks equality operator.
 */
template<> template<>
void testobject::test<4>()
{
    set_test_name("Checks equality operator.");

    string path(TEMP_DIR "primary_header.bin");
    ifstream ifs(path.c_str(), ios::in | ios::binary);
    PrimaryHeader header2;
    header2.decode(ifs);
    ifs.close();

    ensure(header2 == header);
}

/**
 * Checks de-serialization of CoBo 3.0 frame primary header.
 */
template<> template<>
void testobject::test<5>()
{
    set_test_name("Checks de-serialization of CoBo 3.0 frame primary header");

    PrimaryHeader coboPrimaryHeader(utl::BigEndian);
    coboPrimaryHeader.set_blobNess(false);
    coboPrimaryHeader.set_blkSize_B(64u);
    coboPrimaryHeader.set_frameSize_B(128u);
    coboPrimaryHeader.set_dataSource(0);
    coboPrimaryHeader.set_frameType(1u);
    coboPrimaryHeader.set_revision(1u);

    string path(TEMP_DIR "cobo_primary_header.bin");
    ofstream ofs(path.c_str(), ios::out | ios::binary);
    coboPrimaryHeader.encode(ofs);
    ofs.close();

    // Check size of binary file
    ifstream ifs(path.c_str(), ios::in | ios::binary);
    size_t begin = ifs.tellg();
    ifs.seekg (0, ios::end);
    size_t end = ifs.tellg();
    size_t size = end - begin;
    ensure_equals("Size of primary header", size, 8u);

    ifs.seekg (0, ios::beg);
    std::auto_ptr< PrimaryHeader > decodedHeader = PrimaryHeader::decodePrimaryHeader(ifs);
    ifs.close();

    ensure_equals("isBigEndian", decodedHeader->isBigEndian(), true);
    ensure_equals("endianness", decodedHeader->endianness(), utl::BigEndian);
    ensure_equals("blkSize_B", decodedHeader->blkSize_B(), 64u);
    ensure_equals("isBlobFrame", decodedHeader->isBlobFrame(), false);
    ensure_equals("frameSize_blk", decodedHeader->frameSize_blk(), 2u);
    ensure_equals("frameSize_B", decodedHeader->frameSize_B(), 128u);
    ensure_equals("dataSource", decodedHeader->dataSource(), 0);
    ensure_equals("frameType", decodedHeader->frameType(), 1u);
    ensure_equals("revision", decodedHeader->revision(), 1u);
}

} // namespace tut


