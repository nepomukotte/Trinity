/* =============================================================================
 * PrimaryHeader.h
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
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#ifndef mfm_PrimaryHeader_h_INCLUDED
#define mfm_PrimaryHeader_h_INCLUDED

#include "mfm/Common.h"
#include <utl/Rounding.hpp>
#include <utl/BinIO.hpp>
#include <boost/cstdint.hpp>
#include <iostream>
#include <memory>

namespace mfm
{
//______________________________________________________________________
class FrameFormat;
//______________________________________________________________________
/**
 * This class represents the minimal header that is common to all types of
 * MFM2.2 frames.
 */
class PrimaryHeader
{
public:
	/**
	 * Size in bytes of serialized header as defined by MFM2.x
	 */
	static const unsigned SPEC_SIZE_B = 8;
public:
	PrimaryHeader(utl::Endianness byteOrder = utl::BigEndian);
	PrimaryHeader(PrimaryHeader const & primaryHeader);
	PrimaryHeader(std::istream& dataIn);
	virtual ~PrimaryHeader();

	///@{ Getters: see corresponding attribute descriptions
    bool     isBigEndian()   const { return bigEndian_;     }
    utl::Endianness endianness()  const { return bigEndian_?utl::BigEndian:utl::LittleEndian;}
	bool     isBlobFrame()   const { return blobNess_;      }
	uint16_t blkSize_B()     const { return blkSize_B_;     }
    uint32_t frameSize_blk() const { return frameSize_blk_; }
	uint64_t frameSize_B()   const { return ((uint64_t) frameSize_blk()) * blkSize_B(); }
    uint8_t  dataSource()    const { return dataSource_;    }
    uint16_t frameType()     const { return frameType_;     }
    uint8_t  revision()      const { return revision_;      }
	///@}

	///@{ Setters: see corresponding attribute descriptions
    void set_bigEndian(bool _bigEndian)             { bigEndian_ = _bigEndian; }
    void set_endianness(utl::Endianness _endianness){ set_bigEndian(_endianness == utl::BigEndian); }
	void set_blobNess(bool _blobNess)               { blobNess_ = _blobNess; }
    virtual void set_blkSize_B(uint16_t _blkSize_B)         {  blkSize_B_ = _blkSize_B; }
	virtual void set_frameSize_blk(uint32_t _frameSize_blk) { frameSize_blk_ = _frameSize_blk; }
	virtual void set_frameSize_B(uint64_t _frameSize_B); // See inlined further
	void set_dataSource(uint8_t _dataSource)        { dataSource_ = _dataSource; }
	void set_frameType(uint16_t _frameType)         { frameType_ = _frameType; }
	void set_revision(uint8_t _revision)            { revision_ = _revision; }
	///@}

	///@{ Comparison operators
    virtual bool operator==(const PrimaryHeader& rvalue) const;
    virtual bool operator!=(const PrimaryHeader& rvalue) const;
	///@}

	virtual void decode(std::istream& dataIn);
	static std::auto_ptr<PrimaryHeader> decodePrimaryHeader(std::istream& dataIn);
	static std::auto_ptr<PrimaryHeader> createPrimaryHeader(FrameFormat const & format);
	virtual void encode(std::ostream& dataOut) const;
	virtual void print() const;
	static uint16_t minimal_blkSize_B(const uint64_t frameSize_B);
private:
	bool     bigEndian_;       ///< True if frame format is big endian.
	bool     blobNess_;        ///< Whether the frame is a blob or not.
	uint16_t blkSize_B_;       ///< Size of a unit block in bytes.
	uint32_t frameSize_blk_;   ///< Frame size in blocks.
	uint8_t  dataSource_;      ///< Data source.
	uint16_t frameType_;       ///< Actual frame type (defined by application).
	uint8_t  revision_;        ///< Revision number of that particular type.
};

/**
 * Setter for the total frame size with an argument in Bytes.
 * Resizes the frame so that it is made of an integral number of units blocks,
 * and its size in Bytes is greater or equal to the given size @param _frameSize_B .
 */
inline void PrimaryHeader::set_frameSize_B(uint64_t _frameSize_B)
{
	set_frameSize_blk(utl::divide<utl::RoundUp>(_frameSize_B, (uint64_t) blkSize_B()));
}
//______________________________________________________________________
}

#endif /* mfm_PrimaryHeader_h_INCLUDED */
