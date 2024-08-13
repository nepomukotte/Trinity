/*
 * $Id: ItemIterator.h 1039 2013-02-13 08:34:27Z psizun $
 * @file ItemIterator.h
 * @created 8 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: patrick.sizun@cea.fr
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mfm_ItemIterator_h_INCLUDED
#define mfm_ItemIterator_h_INCLUDED

#include "mfm/Item.h"

namespace mfm {
// Forward declarations
class Frame;
//_________________________________________________________________________________________________
/**
 * A class to iterate through the items of a Frame.
 */
class ItemIterator
{
private:
	class ItemIteratorImpl
	{
	public:
		ItemIteratorImpl(Frame & frame);
		virtual ItemIteratorImpl* clone() = 0;
		virtual bool hasNext() const;
		virtual ItemIteratorImpl & next() = 0;
		Item & item() { return currentItem_; }
		Item const & item() const { return currentItem_; }
	protected:
		Frame & frame_; ///< Frame through which to iterate.
		Item currentItem_; ///< Current item.
		bool isInitialized_;
	};
	class BasicFrameItemIteratorImpl : public ItemIteratorImpl
	{
	public:
		BasicFrameItemIteratorImpl(Frame & frame);
		BasicFrameItemIteratorImpl* clone();
		BasicFrameItemIteratorImpl & next();
	};
	class LayeredFrameItemIteratorImpl : public ItemIteratorImpl
	{
	public:
		LayeredFrameItemIteratorImpl(Frame & frame);
		LayeredFrameItemIteratorImpl* clone();
		LayeredFrameItemIteratorImpl & next();
	};
	class BlobFrameItemIteratorImpl : public ItemIteratorImpl
	{
	public:
		BlobFrameItemIteratorImpl(Frame & frame);
		BlobFrameItemIteratorImpl* clone();
		BlobFrameItemIteratorImpl & next();
	};
public:
	ItemIterator();
	ItemIterator(Frame & frame);
	ItemIterator(const ItemIterator& iter);
	virtual ~ItemIterator();
	ItemIterator& operator=(const ItemIterator& iter);
	bool isValid() const;

	virtual bool hasNext() const;
	ItemIterator& next();
	ItemIterator& operator ++();

	Item& operator *();
	Item* operator ->();
	const Item& operator *() const;
	const Item* operator ->() const;
private:
	std::auto_ptr< ItemIteratorImpl > iterator_; ///< Iterator implementation.
};
//_________________________________________________________________________________________________
} /* namespace mfm */
#endif /* mfm_ItemIterator_h_INCLUDED */
