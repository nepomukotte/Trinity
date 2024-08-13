/*
 * $Id: BitField.cpp 1039 2013-02-13 08:34:27Z psizun $
 * @file BitField.cpp
 * @created 28 mars 2012
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

#include "mfm/BitField.h"
#include "mfm/Field.h"
#include "mfm/Exception.h"

namespace mfm {
//______________________________________________________________________
/**
 * Default constructor.
 */
BitField::BitField() : field_(), pos_b_(0), width_b_(0)
{
}
//______________________________________________________________________
/**
 * Constructs a bit field.
 */
BitField::BitField(Field const & parentField, const size_t & pos_b, const size_t & width_b)
	: field_(parentField), pos_b_(pos_b), width_b_(width_b)
{
	;
}
//______________________________________________________________________
/**
 * Copy constructor.
 */
BitField::BitField(BitField const & r)
	: field_(r.field_), pos_b_(r.pos_b_), width_b_(r.width_b_)
{
	;
}
//______________________________________________________________________
/**
 * Assignment operator.
 */
BitField& BitField::operator=(const BitField & r)
{
	field_ = r.field_;
	pos_b_ = r.pos_b_;
	width_b_ = r.width_b_;
	return *this;
}
//______________________________________________________________________
/**
 * Destructor.
 */
BitField::~BitField()
{
	;
}
//______________________________________________________________________
template <>
DynamicBitset BitField::value() const
{
	// Create bit set for bitfield, which all bits equal  to zero
	DynamicBitset bitfield(width_b_);

	// Iterate through bits of full field
	DynamicBitset pattern = field_.value<DynamicBitset>();
	DynamicBitset::size_type index = pattern.find_first();
	while (index != DynamicBitset::npos and index < pos_b_ + width_b_)
	{
		if (index >= pos_b_)
			bitfield.set(index - pos_b_, true);
		index = pattern.find_next(index);
	}
	return bitfield;
}
//______________________________________________________________________
/**
 * Sets a bit field, using a mfm::DynamicBitset.
 * @param bits Value of the bit field.
 */
template <>
void BitField::setValue(const DynamicBitset & bits)
{
	// Check consistency
	if (bits.size() != width_b_)
	{
		throw mfm::OutOfRangeError();
	}

	// Get current bits for full field
	DynamicBitset pattern = field_.value<DynamicBitset>();

	// Modify bits according to bitfield
	for (DynamicBitset::size_type offset=0; offset < bits.size(); ++offset)
	{
		pattern.set(pos_b_ + offset, bits.test(offset));
	}

	// Re-encode field
	field_.setValue(pattern);
}
//______________________________________________________________________
} /* namespace mfm */
