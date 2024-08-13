/*
 * $Id: FrameFormat.h 1800 2015-10-15 13:17:30Z psizun $
 * @file FrameFormat.h
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

#ifndef mfm_FrameFormat_h_INCLUDED
#define mfm_FrameFormat_h_INCLUDED

#include <mfm/Exception.h>
#include <mfm/Common.h>
#include <mfm/PrimaryHeader.h>
#include <CCfg/Io/Document.h>
#include <CCfg/Object.h>
#include <CCfg/CConfig.h>
#include <utl/BinIO.hpp> // utl::Endianness
#include <map>
#include <string>
#include <memory> // auto_ptr
#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>


namespace mfm {
//______________________________________________________________________
class Frame;
//______________________________________________________________________
/**
 * Exception thrown when a frame format is invalid.
 */
class InvalidFrameFormat : public ::mfm::Exception
{
public:
	InvalidFrameFormat();
	~InvalidFrameFormat() throw();
private:
	std::string formatMessage() const;
};
//______________________________________________________________________
/**
 * This class gives access to the contents of the CompoundConfig describing a frame format.
 */
class FrameFormat
{
public:
	FrameFormat();
	FrameFormat(const std::string & path, const std::string & formatName, int revision);
	virtual ~FrameFormat();
	std::auto_ptr<Frame> createFrame() const;
	bool isValid() const;
	uint16_t frameType() const;
	std::string formatName() const;
	size_t revision() const;
	::mfm::FrameKind frameKind() const;
	::utl::Endianness endianness() const;
	size_t blkSize_B() const;
	void set_blkSize_B(size_t);
	size_t headerSize_blk() const;
	size_t itemSize_B() const;
	uint32_t frameSize_blk() const;
	void findItemField(const std::string & name, size_t & offset_B, size_t & size_B) const;
	std::string findItemField(size_t const & offset_B, size_t const & size_B) const;
	void findHeaderField(const std::string & name, size_t & offset_B, size_t & size_B) const;
	std::string findHeaderField(size_t const & offset_B, size_t const & size_B) const;
	void findItemBitField(std::string const & fieldName, std::string const & bitFieldName, size_t & pos_b, size_t & width_b) const;
	std::string findItemBitField(std::string const & fieldName, size_t const & pos_b, size_t const & width_b) const;
	void findHeaderBitField(std::string const & fieldName, std::string const & name, size_t & pos_b, size_t & width_b) const;
private:
	void clear();
	void loadFromCConfig(const std::string & path);
	struct BitFieldInfo
	{
		size_t offset_b;
		size_t width_b;
	};
	typedef std::map< std::string, BitFieldInfo >  BitFieldMap;
	struct FieldInfo
	{
		size_t offset_B;
		size_t size_B;
		BitFieldMap bitFields;
	};
	typedef std::map< std::string, FieldInfo >  FieldMap;
private:
	std::string formatName_; ///< String identifier of the frame format.
	PrimaryHeader primaryHeader_;
	uint16_t headerSize_blk_;
	uint16_t itemSize_B_;
	FieldMap headerFields;
	FieldMap itemFields;
	bool isValid_;
};
//______________________________________________________________________
} /* namespace mfm */
#endif /* mfm_FrameFormat_h_INCLUDED */
