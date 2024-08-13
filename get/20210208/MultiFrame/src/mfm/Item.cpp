/*
 * $Id: Item.cpp 1039 2013-02-13 08:34:27Z psizun $
 * @file Item.cpp
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

#include "mfm/Item.h"
#include "mfm/Exception.h"
#include "mfm/Frame.h"
#include "mfm/FrameFormat.h"
#include <iostream>
#include <sstream>

namespace mfm {
//______________________________________________________________________
/**
 * Default constructor.
 * Constructs an invalid item.
 */
Item::Item() : framePtr_(0), itemIndex_(0)
{
}
//______________________________________________________________________
Item::Item(Frame* _frame, size_t _itemIndex)
	: framePtr_(_frame), itemIndex_(_itemIndex),
	  serializer_(_frame->serializer(), framePtr_->itemSize_B(itemIndex_), framePtr_->itemOffset(itemIndex_))
{
	;
}
//______________________________________________________________________
Item::Item(Frame* frame, size_t _itemIndex, const uint64_t & itemSize_B, const uint64_t itemOffset_B)
	: framePtr_(frame), itemIndex_(_itemIndex), serializer_(frame->serializer(), itemSize_B, itemOffset_B)
{
	;
}
//______________________________________________________________________
Item::~Item()
{
	;
}
//______________________________________________________________________
bool Item::isValid() const
{
	return (framePtr_ != 0 and itemIndex_ < frame().itemCount());
}
//______________________________________________________________________
Header const & Item::header() const
{
	return framePtr_->header();
}
//______________________________________________________________________
Frame const & Item::frame() const
{
	return *framePtr_;
}
//______________________________________________________________________
size_t Item::index() const { return itemIndex_; }
//______________________________________________________________________
/*
 * (see inherited doc)
 */
std::string Item::fieldName(size_t const & fieldAbsOffset_B, size_t const & fieldSize_B) const
{
	// Compute field offset w.r.t to beginning of item
	size_t fieldOffsetB = fieldAbsOffset_B - offset_B();
	return frame().findFormat().findItemField(fieldOffsetB, fieldSize_B);
}
//______________________________________________________________________
/*
 * (see inherited doc)
 */
void Item::findBitField(std::string const & fieldName, std::string const & bitFieldName,
		size_t & pos_b, size_t & width_b) const
{
	return frame().findFormat().findItemBitField(fieldName, bitFieldName, pos_b, width_b);
}
//______________________________________________________________________
/**
 * @return Size of the item, in Bytes.
 */
uint64_t Item::size_B() const
{
	return serializer_.size_B();
}
//______________________________________________________________________
/**
 * Returns the offset of this item w.r.t. the beginning of its frame.
 * @return Offset of the item, in Bytes.
 */
uint64_t Item::offset_B() const
{
	return serializer_.offset_B();
}
//______________________________________________________________________
/**
 * Returns the field inside this item starting at given offset and with given size.
 * @param fieldOffset_B Offset of the field, in Bytes, w.r.t. the beginning of this item.
 * @param fieldSize_B Size of the field, in Bytes.
 * @return Returns a manipulator for the field.
 */
Field Item::field(size_t const & fieldOffset_B, size_t const & fieldSize_B)
{
	if (fieldOffset_B + fieldSize_B > size_B())
	{
		std::ostringstream oss;
		oss << "accessing field with (offset=" << fieldOffset_B << " B, size=" << fieldSize_B << " B) in item of " << size_B() << " B";
		throw mfm::OutOfRangeError(oss.str());
	}
	return framePtr_->field(this, serializer_.offset_B() + fieldOffset_B, fieldSize_B);
}
//______________________________________________________________________
/**
 * Returns the field inside this item with given name.
 * @param name Name of the field.
 * @return Returns a manipulator for the field.
 */
Field Item::field(std::string const & name)
{
	size_t fieldPos, fieldSize;
	framePtr_->findItemField(name, itemIndex_, fieldPos, fieldSize);
	return framePtr_->field(this, fieldPos, fieldSize);
}
//______________________________________________________________________
/**
 * Returns a field within the item with same offset and size as the input field.
 * It may be used to access the same field within another item.
 */
Field Item::field(const Field & f)
{
	return field(f.relativeOffset_B(), f.size_B());
}
//______________________________________________________________________
} /* namespace mfm */
