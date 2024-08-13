/*
 * $Id: Item.h 1039 2013-02-13 08:34:27Z psizun $
 * @file Item.h
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

#ifndef mfm_Item_h_INCLUDED
#define mfm_Item_h_INCLUDED

#include <mfm/Field.h>
#include <mfm/Common.h>
#include <mfm/AbstractFieldContainer.h>

namespace mfm {
//______________________________________________________________________
class Frame;
/**
 * This class represents a (basic) Item within a MFM frame.
 */
class Item : public AbstractFieldContainer
{
public:
	Item();
	Item(Frame* frame, size_t _itemIndex=0);
	Item(Frame* frame, size_t _itemIndex, const uint64_t & itemSize_B, const uint64_t itemOffset_B);
	virtual ~Item();
	bool isValid() const;
	Header const & header() const;
	Frame const & frame() const;
	size_t index() const;
	std::string fieldName(size_t const & offset_B, size_t const & size_B) const;
	void findBitField(std::string const & fieldName, std::string const & bitFieldName,
			size_t & pos_b, size_t & width_b) const;
	uint64_t size_B() const;
	uint64_t offset_B() const;
	Field field(size_t const & pos, size_t const & size);
	Field field(std::string const & name="");
	Field field(const Field & f);
private:
	Frame* framePtr_; ///< Frame this item belongs to.
	size_t itemIndex_; ///< Index of this item within its frame.
	Serializer serializer_; ///< Access to buffer.
};
//______________________________________________________________________
} /* namespace mfm */

#endif /* mfm_Item_h_INCLUDED */
