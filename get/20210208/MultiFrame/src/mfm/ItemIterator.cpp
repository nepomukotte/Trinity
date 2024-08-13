/*
 * $Id: ItemIterator.cpp 1039 2013-02-13 08:34:27Z psizun $
 * @file ItemIterator.cpp
 * @created 8 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors:
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

#include "mfm/ItemIterator.h"
#include "mfm/Frame.h"

namespace mfm {
//_________________________________________________________________________________________________
/** Default constructor.
 * Constructs an invalid iterator.
 */
ItemIterator::ItemIterator()
	: iterator_(0)
{
}
//_________________________________________________________________________________________________
/** Default constructor.
 * Constructs a iterator to iterate through items of a frame.
 * @param frame The frame to iterate through.
 */
ItemIterator::ItemIterator(Frame & frame)
{
	if (frame.header().isBasicFrame())
	{
		iterator_.reset(new BasicFrameItemIteratorImpl(frame));
	}
	else if (frame.header().isLayeredFrame())
	{
		iterator_.reset();
	}
	else
	{
		iterator_.reset();
	}
}
//_________________________________________________________________________________________________
/**
 * Copy constructor.
 * @param value The iterator to copy.
 */
ItemIterator::ItemIterator(const ItemIterator& value)
	: iterator_(value.iterator_->clone())
{
}
//_________________________________________________________________________________________________
/**
 * Assignment operator.
 * @param value The iterator to copy.
 */
ItemIterator& ItemIterator::operator =(const ItemIterator& value)
{
	iterator_.reset(value.iterator_->clone());
	return *this;
}
//_________________________________________________________________________________________________
/**
 * Destructor.
 */
ItemIterator::~ItemIterator()
{
}
//_________________________________________________________________________________________________
/**
 * Indicates whether the iterator is valid.
 */
bool ItemIterator::isValid() const
{
	return (iterator_.get() != 0);
}
//_________________________________________________________________________________________________
/**
 * Indicates whether the iterator has reached the last item in the frame.
 */
bool ItemIterator::hasNext() const
{
	return (isValid()?iterator_->hasNext():false);
}
//_________________________________________________________________________________________________
/**
 * Moves to the next item.
 * @return Returns a reference to the iterator.
 */
ItemIterator& ItemIterator::next()
{
	assert(isValid());
	iterator_->next();
	return *this;
}
//_________________________________________________________________________________________________
/**
 * Moves to the next item.
 * @return Returns a reference to the iterator.
 */
ItemIterator& ItemIterator::operator ++()
{
	return next();
}
//_________________________________________________________________________________________________
/**
 * Retrieves the current item in the iterator.
 * \note
 * The function "next()" must have been called at least once before.
 * @return Returns the current item in the iterator.
 */
Item& ItemIterator::operator *()
{
	assert(isValid());
	return iterator_->item();
}
//_________________________________________________________________________________________________
/**
 * Allows to directly call methods on the current item in the frame.
 * @return Returns a pointer to the current item in the frame.
 */
Item* ItemIterator::operator ->()
{
	assert(isValid());
	return &iterator_->item();
}
//_________________________________________________________________________________________________
/**
 * Retrieves the current item in the iterator (const version).
 * @return Returns the current item in the iterator.
 */
const Item& ItemIterator::operator *() const
{
	assert(isValid());
	return iterator_->item();
}
//_________________________________________________________________________________________________
/**
 * Allows to directly call methods on the current item in the frame (const version).
 * @return Returns a pointer to the current object in the view.
 */
const Item* ItemIterator::operator ->() const
{
	assert(isValid());
	return &iterator_->item();
}
//_________________________________________________________________________________________________
ItemIterator::ItemIteratorImpl::ItemIteratorImpl(Frame & frame)
	: frame_(frame), currentItem_(), isInitialized_(false)
{
}
//_________________________________________________________________________________________________
bool ItemIterator::ItemIteratorImpl::hasNext() const
{
	return (not isInitialized_ or currentItem_.index() + 1 < frame_.header().itemCount());
}
//_________________________________________________________________________________________________
ItemIterator::BasicFrameItemIteratorImpl::BasicFrameItemIteratorImpl(Frame & frame)
	: ItemIteratorImpl(frame)
{
}
//_________________________________________________________________________________________________
ItemIterator::BasicFrameItemIteratorImpl* ItemIterator::BasicFrameItemIteratorImpl::clone()
{
	return new BasicFrameItemIteratorImpl(*this);
}
//_________________________________________________________________________________________________
ItemIterator::BasicFrameItemIteratorImpl & ItemIterator::BasicFrameItemIteratorImpl::next()
{
	//assert(hasNext());

	if (not isInitialized_)
	{
		isInitialized_ = true;
		currentItem_ = Item(&frame_, 0u, frame_.header().itemSize_B(), frame_.header().headerSize_B());
	}
	else
	{
		const size_t newIndex = currentItem_.index() + 1;
		const size_t newOffset_B =  currentItem_.offset_B() + frame_.header().itemSize_B();
		currentItem_ = Item(&frame_, newIndex, frame_.header().itemSize_B(), newOffset_B);
	}

	return *this;
}
//_________________________________________________________________________________________________
ItemIterator::LayeredFrameItemIteratorImpl::LayeredFrameItemIteratorImpl(Frame & frame)
	: ItemIteratorImpl(frame)
{
}
//_________________________________________________________________________________________________
ItemIterator::LayeredFrameItemIteratorImpl* ItemIterator::LayeredFrameItemIteratorImpl::clone()
{
	return new LayeredFrameItemIteratorImpl(*this);
}
//_________________________________________________________________________________________________
ItemIterator::LayeredFrameItemIteratorImpl & ItemIterator::LayeredFrameItemIteratorImpl::next()
{
	//assert(hasNext());

	if (not isInitialized_)
	{
		isInitialized_ = true;
		currentItem_ = Item(&frame_, 0u, frame_.header().itemSize_B(), frame_.header().headerSize_B());
	}
	else
	{
		const size_t newIndex = currentItem_.index() + 1;
		const size_t newOffset_B = currentItem_.offset_B() + currentItem_.size_B();
		currentItem_ = Item(&frame_, newIndex, frame_.frameSize_B(newIndex, newOffset_B), newOffset_B);
	}

	return *this;
}
//_________________________________________________________________________________________________
ItemIterator::BlobFrameItemIteratorImpl::BlobFrameItemIteratorImpl(Frame & frame)
	: ItemIteratorImpl(frame)
{
}
//_________________________________________________________________________________________________
ItemIterator::BlobFrameItemIteratorImpl* ItemIterator::BlobFrameItemIteratorImpl::clone()
{
	return new BlobFrameItemIteratorImpl(*this);
}
//_________________________________________________________________________________________________
ItemIterator::BlobFrameItemIteratorImpl & ItemIterator::BlobFrameItemIteratorImpl::next()
{
	//assert(hasNext());

	isInitialized_ = true;
	currentItem_ = Item(&frame_, 0u, frame_.header().frameSize_B() - frame_.header().headerSize_B(), frame_.header().headerSize_B());

	return *this;
}
//_________________________________________________________________________________________________
} /* namespace mfm */
