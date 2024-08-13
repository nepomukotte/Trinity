/* =============================================================================
 * PrimaryHeader.cpp
 * -----------------------------------------------------------------------------
 * class mfm::PrimaryHeader
 * Created on: 16 March 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the Multiframe Metaformat (MFM2.x) Package of the
 * Mordicus Real-Time Software Development Framework, whose purpose is to encode
 * structured data into binary,recursive, configurable and evolving application
 * formats that allow both backward and forward compatibility to the software
 * that use this data.
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
 * commercial purposes. For thism you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#include "mfm/PrimaryHeader.h"
#include "mfm/FrameFormat.h"
#include "utl/BitFieldHelper.hpp"
#include "utl/BinIO.hpp"
#include <iostream>
#include <cmath>

using namespace utl;
using std::ios_base;

namespace mfm
{

/**
 * Default constructor builds an empty "Blob" frame primary header.
 */
PrimaryHeader::PrimaryHeader(utl::Endianness byteOrder)
	: bigEndian_(utl::BigEndian == byteOrder), // By default, serialized data is in network order.
	  blobNess_(true),               // Default Meta-type of the Frame is the simplest: BLOB
	  blkSize_B_(1u),                // By default, block size is 1 byte (= 2^0).
	  frameSize_blk_(SPEC_SIZE_B),   // By default, frame is reduced to one PrimaryHeader, its size aligned with block size.
	  dataSource_(0u),               // Defaults to 0. Meaning is application dependent.
	  frameType_(0u),                // Defaults to 0. Meaning is application dependent.
	  revision_(0u)                  // Defaults to 0. Meaning is application dependent.
{
		;
}

/**
 * Copy constructor.
 */
PrimaryHeader::PrimaryHeader(PrimaryHeader const & r)
	: bigEndian_(r.isBigEndian()),
	  blobNess_(r.isBlobFrame()),
	  blkSize_B_ (r.blkSize_B()),
	  frameSize_blk_(r.frameSize_blk()),
	  dataSource_(r.dataSource()),
	  frameType_(r.frameType()),
	  revision_(r.revision())
{
	;
}

/**
 * This constructor fills the PrimaryHeader components according to the serialized
 * data that comes from the specified input. Decoding is based on MFM2.x spec.
 * @param dataIn Input data stream.
 */
PrimaryHeader::PrimaryHeader(std::istream& dataIn)
	: bigEndian_(utl::BigEndian),    // By default, serialized data is in network order.
	  blobNess_(true),               // Default Meta-type of the Frame is the simplest: BLOB
	  blkSize_B_(1u),                // By default, block size is 1 byte (= 2^0).
	  frameSize_blk_(SPEC_SIZE_B),   // By default, frame is reduced to one PrimaryHeader, its size aligned with block size.
	  dataSource_(0u),               // Defaults to 0. Meaning is application dependent.
	  frameType_(0u),                // Defaults to 0. Meaning is application dependent.
	  revision_(0u)                  // Defaults to 0. Meaning is application dependent.
{
	decode(dataIn);
}

/** Destructor. */
PrimaryHeader::~PrimaryHeader()
{
	;
}

/**
 * Equality comparison operator.
 * @param rvalue Primary header to compare this header to.
 * @return Whether both headers are identical.
 */
bool PrimaryHeader::operator==(const PrimaryHeader& rvalue) const
{
    return bigEndian_ == rvalue.bigEndian_ and
            blobNess_ == rvalue.blobNess_ and
            blkSize_B_ == rvalue.blkSize_B_ and
            frameSize_blk_ == rvalue.frameSize_blk_ and
            dataSource_ == rvalue.dataSource_ and
            frameType_ == rvalue.frameType_ and
            revision_ == rvalue.revision_;
}

/**
 * Unequality comparison operator.
 */
bool PrimaryHeader::operator!=(const PrimaryHeader& rvalue) const
{
    return !operator==(rvalue);
}

/**
 * @fn bool PrimaryHeader::isBlobFrame() const
 * Returns true if this is a Blob Frame (See MFM2.x).
 */

/**
 * @fn uint64_t PrimaryHeader::frameSize_B() const
 * Returns size of frame in bytes (as opposed to size in blocks).
 */

/**
 * @fn void PrimaryHeader::set_frameSize_B(uint64_t _frameSize_B)
 * Sets the frame size. Since the actual frame size has to be aligned with the
 * frame's block size, this call will set the frame to the smallest size that is
 * aligned with the block size and greater or equal than the specified size.
 * @param _frameSize_B Frame size in bytes. Actually treated as a hint since the
 * final value can be greater than what is specified by this parameter.
 */

/**
 * Decodes the serialized data that comes from the specified input and fills
 * PrimaryHeader object accordingly.
 * @param dataIn Input data stream.
 * @return true if data was byte swapped during decoding, false otherwise.
 */
void PrimaryHeader::decode(std::istream& dataIn)
{
	// Make sure that input stream is in exception mode
	dataIn.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);

	// First read out metaType byte (see MFM2.x specification)
	uint8_t metaType;
	BinIO::decode<uint8_t>(metaType, dataIn);

	// Decode meta-information in first byte
	bigEndian_ =            not BitFieldHelper<uint8_t>::getBit  (metaType, 7);
	blkSize_B_   = pow2<uint16_t>(BitFieldHelper<uint8_t>::getField(metaType, 0, 4));
	blobNess_    =                BitFieldHelper<uint8_t>::getBit  (metaType, 6);

	if (bigEndian_)
	{
		BinIO::decode<BigEndian, uint32_t, 3>(frameSize_blk_, dataIn);
        BinIO::decode<BigEndian, uint8_t>(dataSource_, dataIn);
		BinIO::decode<BigEndian, uint16_t>(frameType_, dataIn);
		BinIO::decode<BigEndian, uint8_t>(revision_, dataIn);
	}
	else
	{
		BinIO::decode<LittleEndian, uint32_t, 3>(frameSize_blk_, dataIn);
        BinIO::decode<LittleEndian, uint8_t>(dataSource_, dataIn);
		BinIO::decode<LittleEndian, uint16_t>(frameType_, dataIn);
		BinIO::decode<LittleEndian, uint8_t>(revision_, dataIn);
	}
}

/**
 * Reads from an input stream to decode a frame primary header.
 * @param dataIn Input stream.
 * @return Returns the primary header created.
 */
std::auto_ptr<PrimaryHeader> PrimaryHeader::decodePrimaryHeader(std::istream& dataIn)
{
	std::auto_ptr<PrimaryHeader> phPtr(new PrimaryHeader);
	phPtr->decode(dataIn);
	return phPtr;
}

/**
 * Encodes PrimaryHeader elements into serialized data stream according to MFM2.x
 * specifications.
 * @param dataOut Output data stream.
 */
void PrimaryHeader::encode(std::ostream& dataOut) const
{
	// Make sure that output stream is in exception mode
	dataOut.exceptions(ios_base::eofbit | ios_base::failbit | ios_base::badbit);

	// Prepare meta-information in first byte (see MFM2.x specification)
	uint8_t metaType = 0;
	BitFieldHelper<uint8_t>::setBit  (metaType, 7, not bigEndian_);
	BitFieldHelper<uint8_t>::setField(metaType, 0, 4, logpow2<uint8_t>(blkSize_B_));
	BitFieldHelper<uint8_t>::setBit  (metaType, 6, blobNess_);

	BinIO::encode<uint8_t>(metaType, dataOut);

	if (bigEndian_)
	{
		BinIO::encode<BigEndian, uint32_t, 3>(frameSize_blk_, dataOut);
		BinIO::encode<BigEndian, uint8_t>(dataSource_, dataOut);
		BinIO::encode<BigEndian, uint16_t>(frameType_, dataOut);
		BinIO::encode<BigEndian, uint8_t>(revision_, dataOut);
	}
	else
	{
		BinIO::encode<LittleEndian, uint32_t, 3>(frameSize_blk_, dataOut);
        BinIO::encode<LittleEndian, uint8_t>(dataSource_, dataOut);
        BinIO::encode<LittleEndian, uint16_t>(frameType_, dataOut);
        BinIO::encode<LittleEndian, uint8_t>(revision_, dataOut);
	}
}
//______________________________________________________________________
/**
 * Creates a primary header adapted to a given frame format description.
 * @param format Frame format description.
 * @return Returns the primary header created.
 */
std::auto_ptr<PrimaryHeader> PrimaryHeader::createPrimaryHeader(FrameFormat const & format)
{
	std::auto_ptr<PrimaryHeader> phPtr(new PrimaryHeader);
	phPtr->set_bigEndian(format.endianness() == BigEndian);
	phPtr->set_blkSize_B(format.blkSize_B());
	phPtr->set_blobNess(format.frameKind() == BLOB);
	if (format.frameKind() == BLOB)
	{
		phPtr->set_frameSize_blk(format.frameSize_blk());
	}
	else
	{
		phPtr->set_frameSize_blk(format.headerSize_blk());
	}
	phPtr->set_dataSource(0u);
	phPtr->set_frameType(format.frameType());
	phPtr->set_revision(format.revision());
	return phPtr;
}
//______________________________________________________________________
/**
 * Prints the contents of this primary header to standard output.
 */
void PrimaryHeader::print() const
{
	std::cout << "PrimaryHeader[bigEndian = "       << bigEndian_
		      << ", blobNess = "      << blobNess_
		      << ", blkSize_B = "     << blkSize_B_
		      << ", frameSize_blk = " << frameSize_blk_
		      << ", dataSource = "    << (short) dataSource_
		      << ", frameType = "     << frameType_
		      << ", revision = "      << (short) revision_
		      << "]";
}
//______________________________________________________________________
/**
 * Computes smallest block size that can be used to store a frame of given size, considering the frame size in blocks is coded on 3 Bytes.
 * @param frameSize_B Frame size, in Bytes.
 * @return Optimal block size, in Bytes.
 */
uint16_t PrimaryHeader::minimal_blkSize_B(const uint64_t frameSize_B)
{
	if (frameSize_B <= PrimaryHeader::SPEC_SIZE_B)
		return 1u;

	uint8_t log2BlkSize_B = std::max<double>(0., std::ceil(std::log(frameSize_B)/std::log(2) - 24.));
	return pow2<uint16_t>(log2BlkSize_B);
}
//______________________________________________________________________
} // namespace mfm
