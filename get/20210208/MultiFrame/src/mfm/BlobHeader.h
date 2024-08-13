/* =============================================================================
 * $Id: BlobHeader.h 716 2012-04-10 15:17:58Z psizun $
 * @file BlobHeader.h
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

#ifndef mfm_BlobHeader_h_INCLUDED
#define mfm_BlobHeader_h_INCLUDED

#include <mfm/Header.h>

namespace mfm {

/**
 * An MFM frame header adapted to a Blob frame.
 */
class BlobHeader: public mfm::Header
{
public:
    BlobHeader(utl::Endianness _byteOrder = utl::BigEndian);
    BlobHeader(PrimaryHeader const &);
    BlobHeader(BlobHeader const &);
    virtual BlobHeader* clone() const;
    virtual ~BlobHeader() { ; }

    ///@{ Getters
    FrameKind metaType() const { return BLOB; }
    bool isBlobFrame() const { return true; }
    uint16_t headerSize_blk() const;
    uint32_t headerSize_B() const;
    uint16_t itemSize_B() const;
    uint32_t itemCount()  const;
    uint32_t headerReserveSize_B() const;
    uint32_t dataReserveSize_B() const;
    ///@}

    ///@{ Setters
    void set_headerSize_blk(uint16_t);
    void set_headerSize_B(uint32_t);
    void set_itemSize_B(uint16_t);
    void set_itemCount(uint32_t);
	void set_dataSize_blk(uint32_t _dataSize_blk);
	void set_dataSize_B(uint64_t _dataSize_B);
    ///@}
};

} /* namespace mfm */
#endif /* mfm_BlobHeader_h_INCLUDED */
