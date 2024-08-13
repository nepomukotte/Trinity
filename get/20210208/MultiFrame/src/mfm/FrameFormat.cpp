/*
 * $Id: FrameFormat.cpp 1800 2015-10-15 13:17:30Z psizun $
 * @file FrameFormat.cpp
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

#include "mfm/FrameFormat.h"
#include "mfm/Common.h"
#include "mfm/Frame.h"
#include <CCfg/Object.h>
#include <CCfg/View/ClassIterator.h>
#include <CCfg/CConfig.h>
#include <utl/BinIO.hpp>
#include <utl/Rounding.hpp>
using namespace CCfg;
#include <sstream>
#include <iostream>

namespace mfm {
//______________________________________________________________________
/**
 * Default constructor.
 * Constructs an invalid frame format.
 */
FrameFormat::FrameFormat() :
		formatName_(), headerSize_blk_(0), itemSize_B_(0), isValid_(false)
{
	;
}
//______________________________________________________________________
/**
 * Constructs a format based on the filepath of its description file.
 */
FrameFormat::FrameFormat(const std::string & path, const std::string & formatName, int revision) :
		formatName_(formatName), headerSize_blk_(0), itemSize_B_(0), isValid_(false)
{
	try
	{
		primaryHeader_.set_revision(revision);
		loadFromCConfig(path);
	}
	catch (CCfg::Io::Document::IoError & e)
	{
		throw InvalidFrameFormat();
	}
}
//______________________________________________________________________
FrameFormat::~FrameFormat()
{
	;
}
//______________________________________________________________________
void FrameFormat::clear()
{
	headerFields.clear();
	itemFields.clear();
}
//______________________________________________________________________
void FrameFormat::loadFromCConfig(const std::string & path)
{
	clear();
	isValid_ = false;
	try
	{
		Io::Document doc;
		CConfig config(&doc.loadFromFile(path));
		CConfig formatCfg = config("Frame", formatName_);
		CConfig revisionCfg = formatCfg("Revision", revision());

		// Get primary header info
		primaryHeader_.set_frameType(formatCfg("typeCode").getIntValue());
		CConfig metaTypeCfg = revisionCfg("Header")("Field", "metaType");
		primaryHeader_.set_blobNess(metaTypeCfg("BitField", "ISBLOB")("value").getIntValue());
		primaryHeader_.set_bigEndian(not metaTypeCfg("BitField", "ISLEND")("value").getIntValue());
		primaryHeader_.set_blkSize_B(pow2<uint16_t>(metaTypeCfg("BitField", "P2BLCK")("value").getIntValue()));
		try
		{
			primaryHeader_.set_frameSize_blk(metaTypeCfg("BitField", "P2BLCK")("value").getIntValue()); // FIXME: ??!!
		}
		catch (const CCfg::CConfig::NotFound &)
		{
		}
		headerSize_blk_ = revisionCfg("Header")("Field", "headerSize")("value").getIntValue();
		itemSize_B_ = revisionCfg("Header")("Field", "itemSize")("value").getIntValue();

		// Get header fields
		ConfigList fields = revisionCfg("Header").getAttributes("Field");
		for (size_t i = 0; i < fields.size(); ++i)
		{
			std::string fieldName = fields[i].getIndex();
			//std::cout << "Found header field '" << fieldName << '\'' << std::endl;
			headerFields[fieldName].offset_B = fields[i]("offset").getIntValue();
			headerFields[fieldName].size_B = fields[i]("size").getIntValue();
			// Get bit fields
			ConfigList bitFields = fields[i].getAttributes("BitField");
			for (size_t j = 0; j < bitFields.size(); ++j)
			{
				std::string bitFieldName = bitFields[j].getIndex();
				//std::cout << "Found header bitfield '" << fieldName << "'.'" << bitFieldName << '\'' << std::endl;
				headerFields[fieldName].bitFields[bitFieldName].offset_b = bitFields[j]("offset").getIntValue();
				headerFields[fieldName].bitFields[bitFieldName].width_b = bitFields[j]("width").getIntValue();
			}
		}

		// Get items
		fields = revisionCfg("Item", "*").getAttributes("Field");
		for (size_t i = 0; i < fields.size(); ++i)
		{
			std::string fieldName = fields[i].getIndex();
			//std::cout << "Found item field '" << fieldName << '\'' << std::endl;
			itemFields[fieldName].offset_B = fields[i]("offset").getIntValue();
			itemFields[fieldName].size_B = fields[i]("size").getIntValue();
			// Get bit fields
			ConfigList bitFields = fields[i].getAttributes("BitField");
			for (size_t j = 0; j < bitFields.size(); ++j)
			{
				std::string bitFieldName = bitFields[j].getIndex();
				//std::cout << "Found item bitfield '" << fieldName << "'.'" << bitFieldName << '\'' << std::endl;
				itemFields[fieldName].bitFields[bitFieldName].offset_b = bitFields[j]("offset").getIntValue();
				itemFields[fieldName].bitFields[bitFieldName].width_b = bitFields[j]("width").getIntValue();
			}
		}

		isValid_ = true;
	}
	catch (CCfg::Io::Document::IoError & e)
	{
		throw InvalidFrameFormat();
	}
	catch (CCfg::Attribute::NotFound & e)
	{
		throw InvalidFrameFormat();
	}
}
//______________________________________________________________________
/**
 * Creates a minimal frame (e.g. without items) corresponding to this format.
 * @return Returns a new Frame respecting this format.
 */
std::auto_ptr<Frame> FrameFormat::createFrame() const
{
    return Frame::create(*this);
}
//______________________________________________________________________
/**
 * Checks whether this format is valid.
 * @return Whether the format is valid.
 */
bool FrameFormat::isValid() const
{
	return isValid_;
}
//______________________________________________________________________
/**
 * Returns the value of the frameType header field.
 * @return Value of the field according to the format description.
 */
uint16_t FrameFormat::frameType() const
{
	return primaryHeader_.frameType();
}
//______________________________________________________________________
/**
 * Returns the value of the frameSize header field.
 * @return Value of the field according to the format description.
 */
uint32_t FrameFormat::frameSize_blk() const
{
	return primaryHeader_.frameSize_blk();
}
//______________________________________________________________________
/**
 * Returns the string identifier of the frame format.
 * @return Name of the format according to the format description.
 */
std::string FrameFormat::formatName() const
{
	return formatName_;
}
//______________________________________________________________________
/**
 * Returns the value of the revision header field.
 * @return Value of the field according to the format description.
 */
size_t FrameFormat::revision() const
{
	return primaryHeader_.revision();
}
//______________________________________________________________________
/**
 * Returns the kind of frame this format corresponds to.
 * @return Kind of frame (basic, layered or blob) this format describes.
 */
FrameKind FrameFormat::frameKind() const
{
	if (primaryHeader_.isBlobFrame())
		return BLOB;
	if (0 == itemSize_B())
		return LAYERED;
	return BASIC;
}
//______________________________________________________________________
/**
 * Returns the default endianness to use for this format.
 * FIXME: Not sure it should be part of a frame format description.
 * @return Byte order to use when creating new frames of this format.
 */
utl::Endianness FrameFormat::endianness() const
{
	return (primaryHeader_.isBigEndian()?utl::BigEndian:utl::LittleEndian);
}
//______________________________________________________________________
/**
 * @return Size of a unit block [Bytes].
 */
size_t FrameFormat::blkSize_B() const
{
	return primaryHeader_.blkSize_B();
}
//______________________________________________________________________
void FrameFormat::set_blkSize_B(size_t new_blkSize_B)
{
	headerSize_blk_ = utl::divide<utl::RoundUp>((uint64_t) headerSize_blk_ * primaryHeader_.blkSize_B(), (uint64_t) new_blkSize_B);
	primaryHeader_.set_blkSize_B(new_blkSize_B);
}
//______________________________________________________________________
/**
 * Returns the default header size, as an integer number of unit blocks.
 * This should only be used to create new frames from scratch.
 * @return Size of a frame header, in unit block units.
 */
size_t FrameFormat::headerSize_blk() const
{
	return headerSize_blk_;
}
//______________________________________________________________________
/**
 * Returns the size of an item, in Bytes.
 * @return Size of an item, in Bytes.
 */
size_t FrameFormat::itemSize_B() const
{
	return itemSize_B_;
}
//______________________________________________________________________
/**
 * Searches for field with given name within the item description.
 * @param fieldName Name of the field to search for.
 * @param[out] offset_B Offset of the field, in Bytes,  w.r.t to the beginning of the item.
 * @param[out] size_B Size of the field, in Bytes.
 * @throws mfm::FieldNotFound Exception thrown if field is not found.
 */
void FrameFormat::findItemField(const std::string & fieldName, size_t & offset_B,
		size_t & size_B) const
{
	//std::cout << "findItemField(" << fieldName << ')' << std::endl;
	FieldMap::const_iterator fieldIter = itemFields.find(fieldName);
	if (fieldIter == itemFields.end())
	{
		throw mfm::FieldNotFound(fieldName);
	}

	offset_B = fieldIter->second.offset_B;
	size_B = fieldIter->second.size_B;
}
//______________________________________________________________________
/**
 * Searches for name of item field with given offset and size.
 * @return Returns the name of the field.
 * @param[in] offset_B Offset of the field, in Bytes,  w.r.t to the beginning of the item.
 * @param[in] size_B Size of the field, in Bytes.
 * @throws mfm::FieldNotFound Exception thrown if field is not found.
 */
std::string FrameFormat::findItemField(size_t const & offset_B,
		size_t const & size_B) const
{
	FieldMap::const_iterator fieldIter;
	for (fieldIter = itemFields.begin(); fieldIter != itemFields.end(); ++fieldIter)
	{
		if (fieldIter->second.offset_B == offset_B and
				fieldIter->second.size_B == size_B)
		return fieldIter->first;
	}

	std::ostringstream oss;
	oss << "with offset=" << offset_B << " and size=" << size_B;
	throw mfm::FieldNotFound(oss.str());
}
//______________________________________________________________________
/**
 * Searches for field with given name within the header description.
 * @param fieldName Name of the field to search for.
 * @param[out] offset_B Offset of the field, in Bytes,  w.r.t to the beginning of the header.
 * @param[out] size_B Size of the field, in Bytes.
 * @throws mfm::FieldNotFound Exception thrown if field is not found.
 */
void FrameFormat::findHeaderField(const std::string & fieldName, size_t & offset_B,
		size_t & size_B) const
{
	//std::cout << "findHeaderField(" << fieldName << ')' << std::endl;
	FieldMap::const_iterator fieldIter = headerFields.find(fieldName);
	if (fieldIter == headerFields.end())
	{
		throw mfm::FieldNotFound(fieldName);
	}

	offset_B = fieldIter->second.offset_B;
	size_B = fieldIter->second.size_B;
}
//______________________________________________________________________
/**
 * Searches for name of header field with given offset and size.
 * @return Returns the name of the field.
 * @param[in] offset_B Offset of the field, in Bytes,  w.r.t to the beginning of the item.
 * @param[in] size_B Size of the field, in Bytes.
 * @throws mfm::FieldNotFound Exception thrown if field is not found.
 */
std::string FrameFormat::findHeaderField(size_t const & offset_B,
		size_t const & size_B) const
{
	FieldMap::const_iterator fieldIter;
	for (fieldIter = headerFields.begin(); fieldIter != headerFields.end(); ++fieldIter)
	{
		if (fieldIter->second.offset_B == offset_B and
				fieldIter->second.size_B == size_B)
		return fieldIter->first;
	}

	std::ostringstream oss;
	oss << "with offset=" << offset_B << " and size=" << size_B;
	throw mfm::FieldNotFound(oss.str());
}
//______________________________________________________________________
/**
 * Searches for bit field with given name within the item field with given name.
 * @param fieldName Name of the field.
 * @param bitFieldName Name of the bit field to search for.
 * @param[out] pos_b Offset of the bit field, in bits.
 * @param[out] width_b Size of the bit field, in bits.
 * @throws mfm::BitFieldNotFound Exception thrown if bit field is not found.
 */
void FrameFormat::findItemBitField(std::string const & fieldName, std::string const & bitFieldName,
		size_t & pos_b, size_t & width_b) const
{
	//std::cout << "findItemBitField(" << fieldName << ',' << bitFieldName << ')' << std::endl;
	// Find field
	FieldMap::const_iterator fieldIter = itemFields.find(fieldName);
	if (fieldIter == itemFields.end())
	{
		throw mfm::BitFieldNotFound(fieldName, bitFieldName);
	}

	// Find bitfield
	const BitFieldMap & bitFields =  fieldIter->second.bitFields;
	BitFieldMap::const_iterator bitFieldIter = bitFields.find(bitFieldName);
	if (bitFieldIter == bitFields.end())
	{
		throw mfm::BitFieldNotFound(fieldName, bitFieldName);
	}

	pos_b = bitFieldIter->second.offset_b;
	width_b = bitFieldIter->second.width_b;
}
//______________________________________________________________________
/**
 * Searches for name of bit field with given position and width, within field with given name.
 * @return Returns the name of the bit field.
 * @param[in] fieldName Name of the item field.
 * @param[in] pos_b Position of the bit field, in bits, within the field.
 * @param[in] width_b Width of the bit field, in bits.
 * @throws mfm::BitFieldNotFound Exception thrown if bit field is not found.
 */
std::string FrameFormat::findItemBitField(std::string const & fieldName, size_t const & pos_b,
		size_t const & width_b) const
{
	// Find field
	FieldMap::const_iterator fieldIter = itemFields.find(fieldName);
	if (fieldIter == itemFields.end())
	{
		throw mfm::FieldNotFound(fieldName);
	}

	const BitFieldMap & bitFields =  fieldIter->second.bitFields;
	BitFieldMap::const_iterator bitFieldIter;
	for (bitFieldIter = bitFields.begin(); bitFieldIter != bitFields.end(); ++bitFieldIter)
	{
		if (bitFieldIter->second.offset_b == pos_b and
			bitFieldIter->second.width_b == width_b)
		return bitFieldIter->first;
	}

	std::ostringstream oss;
	oss << "with offset=" << pos_b << " and width=" << width_b;
	throw mfm::BitFieldNotFound(fieldName, oss.str());
}
//______________________________________________________________________
/**
 * Searches for bit field with given name within the header field with given name.
 * @param fieldName Name of the field.
 * @param bitFieldName Name of the bit field to search for.
 * @param[out] pos_b Offset of the bit field, in bits.
 * @param[out] width_b Size of the bit field, in bits.
 * @throws mfm::BitFieldNotFound Exception thrown if bit field is not found.
 */
void FrameFormat::findHeaderBitField(std::string const & fieldName,
		std::string const & bitFieldName, size_t & pos_b, size_t & width_b) const
{
	//std::cout << "findHeaderBitField(" << fieldName << ',' << bitFieldName << ')' << std::endl;
	// Find field
	FieldMap::const_iterator fieldIter = headerFields.find(fieldName);
	if (fieldIter == headerFields.end())
	{
		throw mfm::BitFieldNotFound(fieldName, bitFieldName);
	}

	// Find bitfield
	const BitFieldMap & bitFields =  fieldIter->second.bitFields;
	BitFieldMap::const_iterator bitFieldIter = bitFields.find(bitFieldName);
	if (bitFieldIter == bitFields.end())
	{
		throw mfm::BitFieldNotFound(fieldName, bitFieldName);
	}

	pos_b = bitFieldIter->second.offset_b;
	width_b = bitFieldIter->second.width_b;
}
//______________________________________________________________________
/**
 * Constructor.
 */
InvalidFrameFormat::InvalidFrameFormat()
{
	;
}
//______________________________________________________________________
/**
 * Destructor.
 */
InvalidFrameFormat::~InvalidFrameFormat() throw ()
{
}
//______________________________________________________________________
/**
 * Format an error message for the actual exception based on its attributes.
 * @return Returns the generated error message.
 */
std::string InvalidFrameFormat::formatMessage() const
{
	std::ostringstream message;
	message << "mfm::InvalidFrameFormat";
	return message.str();
}
//______________________________________________________________________
} /* namespace mfm */
