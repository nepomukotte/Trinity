/**
 * @file RtBitSet.h
 * @created Nov 8, 2012
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the Utilities software project.
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

#ifndef utl_RtBitSet_h_INCLUDED
#define utl_RtBitSet_h_INCLUDED

#include "utl/WordIterator.hpp"
#include "utl/MemOperator.hpp"
#include "utl/BitFieldHelper.hpp"
#include <boost/type_traits/make_signed.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#include <boost/mpl/if.hpp>
#include <boost/cstdint.hpp>
#include <boost/math/common_factor_rt.hpp>
#include <iomanip>
#include <climits>
#include <typeinfo>

#ifdef VXWORKS
#ifdef min
#define UNDEF_VX_MINMAX
#undef min
#undef max
#endif
#endif

namespace utl
{

/**
 * This template structure includes all the meta info we need on a word type.
 * @param W The word type
 */

template <typename W>
struct WordInfo
{
	static const size_t BYTE_COUNT = sizeof(W);
	static const size_t BIT_COUNT = sizeof(W) * CHAR_BIT;
};

struct AtomOperator
{
	struct copy
	{
		template < typename AtomType >
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom  = sourceAtom; return targetAtom; }
	};

	struct bitAnd
	{
		template < typename AtomType >
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom &= sourceAtom; return targetAtom; }
	};

	struct bitOr
	{
		template < typename AtomType >
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom |= sourceAtom; return targetAtom; }
	};

	struct bitXor
	{
		template < typename AtomType >
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom ^= sourceAtom; return targetAtom; }
	};

	struct bitNot
	{
		template < typename AtomType >
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom = ~sourceAtom; return targetAtom; }
	};
};


/**
 * Forward declarations for friend statements and helpers
 */
template < class W > class BitFieldHelper;
template < typename W > struct make_RtBitSet;

/**
 * Template class representing a "word" made of an arbitrary number of bits set at run time (run time bitset)
 * @tparam A Specifies word type on which size the total buffer size of the bitset must be aligned.
 */
template < typename A = unsigned int >
class RtBitSet
{
public:
	/**
	 * Depending on its signed or unsigned character, this bitset will be made of a vector of signed or unsigned words.
	 * Atom refers to the type of this elementary word.
	 */
	typedef A Atom;
	typedef typename boost::make_unsigned <Atom>::type unsigned_Atom; ///< Unsigned version of atom word.
	typedef typename boost::make_signed   <Atom>::type   signed_Atom; ///< Signed version of atom word.
	typedef typename boost::mpl::if_c< boost::is_signed<Atom>::value, signed char, unsigned char >::type Byte;
	typedef WordIterator<Atom> iter_Atom;

	friend class BitFieldHelper< std::vector<Atom> >;

public:
	///@name Constructors, initializer & normalizer
	///@{
	/// Constructor for BitSet over external buffer (no fill)
	RtBitSet(size_t _bitCount, void* x_buffer_ = 0) : i_buffer(0) { resetBuffer(_bitCount, x_buffer_); }

	/// Copy constructor: you need to code it separately from other template constructors
	RtBitSet(const RtBitSet& aBitSet) : i_buffer(0) { resetBuffer(aBitSet.bitCount()); setValue(aBitSet); }

	/// Assignment op: you need to code it separately from other template assignments
	RtBitSet& operator = (const RtBitSet& bitSet) { setValue(bitSet); return *this; }

	template < typename W2 >
	RtBitSet(const W2& initWord) : i_buffer(0) { resetBuffer(WordInfo<W2>::BIT_COUNT); setValue(initWord); }

	template < typename W2 >
	RtBitSet(const std::vector<W2>& initVect) : i_buffer(0) { resetBuffer(WordInfo<W2>::BIT_COUNT*initVect.size()); setValue(initVect); }

	template < typename A2 >
	RtBitSet(const RtBitSet<A2>& initBitSet) : i_buffer(0) { resetBuffer(initBitSet.bitCount()); setValue(initBitSet); }

	virtual ~RtBitSet();

	RtBitSet& resetBuffer(size_t _bitCount, void* _x_buffer = 0);
	RtBitSet& fill(const Atom& initWord);
	RtBitSet& normalize();
	///@}

	///@name Value getters and setters
	///@{
	template < typename W2 >
	void setValue(const W2& w) { apply<AtomOperator::copy, W2>(w); }
	template < typename W2 >
	void setValue(const std::vector<W2>& sourceVect) { apply<AtomOperator::copy, W2>(sourceVect); }
	template < typename A2 >
	void setValue(const RtBitSet<A2>& bitSet) { apply<AtomOperator::copy, A2>(bitSet); }
	template < typename AtomType >
	void setValue(const AtomType* sourceArray, size_t arraySize) { apply<AtomOperator::copy, AtomType>(sourceArray, arraySize); }

	template < typename W2 >
	void getValue(W2& w) const { apply<AtomOperator::copy, W2>(w); }
	template < typename W2 >
	void getValue(std::vector<W2>& targetVect) const { apply<AtomOperator::copy, W2>(targetVect); }
	template < typename A2 >
	void getValue(RtBitSet<A2>& bitSet) const { bitSet.setValue(*this); }
	template < typename AtomType >
	void getValue(AtomType* targetArray, size_t arraySize) const { apply<AtomOperator::copy, AtomType>(targetArray, arraySize); }
	///@}

	///@name Bitwise operations
	///@{
	template < typename W2 >
	void bitAnd(const W2& w) { apply<AtomOperator::bitAnd, W2>(w); }
	void bitAnd(const unsigned_Atom& w) { apply<AtomOperator::bitAnd>(w); }
	void bitAnd(const signed_Atom& w) { apply<AtomOperator::bitAnd>(w); }
	template < typename A2 >
	void bitAnd(const RtBitSet<A2>& bitSet) { apply<AtomOperator::bitAnd, A2>(bitSet); }

	template < typename W2 >
	void bitOr(const W2& w) { apply<AtomOperator::bitOr, W2>(w); }
	void bitOr(const unsigned_Atom& w) { apply<AtomOperator::bitOr>(w); }
	void bitOr(const signed_Atom& w) { apply<AtomOperator::bitOr>(w); }
	template < typename A2 >
	void bitOr(const RtBitSet<A2>& bitSet) { apply<AtomOperator::bitOr, A2>(bitSet); }

	template < typename W2 >
	void bitXor(const W2& w) { apply<AtomOperator::bitXor, W2>(w); }
	void bitXor(const unsigned_Atom& w) { apply<AtomOperator::bitXor>(w); }
	void bitXor(const signed_Atom& w) { apply<AtomOperator::bitXor>(w); }
	template < typename A2 >
	void bitXor(const RtBitSet<A2>& bitSet) { apply<AtomOperator::bitXor, A2>(bitSet); }

	void bitNot();
	void shiftLeft  (size_t nBits);
	void shiftRight (size_t nBits);

	///@name Assignment / transform operators
	///@{
	template < typename W2 >
	RtBitSet& operator   = (const W2& w) { setValue(w); return *this; }
	RtBitSet& operator <<= (size_t nBits) { shiftLeft (nBits); return *this; }
	RtBitSet& operator >>= (size_t nBits) { shiftRight(nBits); return *this; }
	RtBitSet  operator   ~ () const { RtBitSet bitSetTmp(*this); bitSetTmp.bitNot(); return bitSetTmp; }
	RtBitSet  operator  << (size_t nBits) const { RtBitSet bitSetTmp(*this); bitSetTmp <<= nBits; return bitSetTmp; }
	RtBitSet  operator  >> (size_t nBits) const { RtBitSet bitSetTmp(*this); bitSetTmp >>= nBits; return bitSetTmp; }
	template < typename W2 >
	RtBitSet& operator  &= (const W2& w) { bitAnd(w); return *this; }
	template < typename W2 >
	RtBitSet& operator  |= (const W2& w) { bitOr (w); return *this; }
	template < typename W2 >
	RtBitSet& operator  ^= (const W2& w) { bitXor(w); return *this; }
	///@}

	///@name Composition / cast operations
	///@{
	template < typename W2 >
	RtBitSet  operator   & (const W2& w) const { RtBitSet bitSetTmp(*this); bitSetTmp &= w; return bitSetTmp; }
	template < typename W2 >
	RtBitSet  operator   | (const W2& w) const { RtBitSet bitSetTmp(*this); bitSetTmp |= w; return bitSetTmp; }
	template < typename W2 >
	RtBitSet  operator   ^ (const W2& w) const { RtBitSet bitSetTmp(*this); bitSetTmp ^= w; return bitSetTmp; }
	template < typename W2 >
	operator W2 () const { return toScalar<W2>(); }
	///@}

	///@name Comparison operations
	///@{
	template < typename W2 >
	int compare(const W2& w) const;
	template < typename A2 >
	int compare(const RtBitSet<A2>& bitSet) const;
	template < typename W2 >
	bool isEqualTo(const W2& w) const { return compare(w) == 0; }
	template < typename W2 >
	bool isGreaterThan(const W2& w) const { return compare(w) > 0; }
	template < typename W2 >
	bool operator == (const W2& w) const { return isEqualTo(w); }
	template < typename W2 >
	bool operator != (const W2& w) const { return not isEqualTo(w); }
	template < typename W2 >
	bool operator > (const W2& w) const { return isGreaterThan(w); }
	template < typename W2 >
	bool operator >= (const W2& w) const { return isGreaterThan(w) or isEqualTo(w); }
	template < typename W2 >
	bool operator <= (const W2& w) const { return not isGreaterThan(w); }
	template < typename W2 >
	bool operator < (const W2& w) const { return not isGreaterThan(w) and not isEqualTo(w); }
	///@}

	///@name Elementary word (atom) access & iteration operations
	///@{
	iter_Atom itAtomBegin() const { return iter_Atom::getBegin(p_buffer, atomCount_); }
	iter_Atom itAtomEnd  () const { return iter_Atom::getEnd  (p_buffer, atomCount_); }
	iter_Atom itAtomLast () const { return itAtomEnd() - 1; }

	      Atom& atomAt(std::ptrdiff_t iAtom)       { return itAtomBegin()[iAtom]; }
    const Atom& atomAt(std::ptrdiff_t iAtom) const { return itAtomBegin()[iAtom]; }

	      Atom& atomFirst()       { return atomAt(0); }
	const Atom& atomFirst() const { return atomAt(0); }

	      Atom& atomLast ()       { return atomAt(atomCount()-1); }
	const Atom& atomLast () const { return atomAt(atomCount()-1); }
	///@}

	///@name Miscellaneous operations
	///@{
	static unsigned atomBitCount() { return WordInfo<Atom>::BIT_COUNT; }
	static unsigned atomByteCount() { return WordInfo<Atom>::BYTE_COUNT; }

	size_t bitCount () const { return bitCount_; }
	size_t atomCount() const { return atomCount_; }
	size_t byteCount() const { return atomCount() * atomByteCount(); }
	size_t minByteCount() const { return bitCount()%CHAR_BIT == 0 ? bitCount()/CHAR_BIT : 1 + bitCount()/CHAR_BIT; }

	bool isPositive() { return bitCount() == 0 or utl::isPositive(atomLast()); }
	bool isNegative() { return not isPositive(); }

	static bool isBigEndian() { return iter_Atom::isBigEndian(); }

	template < typename W2 >
	W2* bufferPtr() { return reinterpret_cast<W2*>(p_buffer); }

	template < typename W2 >
    const W2* bufferPtr() const { return reinterpret_cast<const W2*>(p_buffer); }

	Atom* buffer() { return bufferPtr<Atom>(); }
	const Atom* buffer() const { return bufferPtr<Atom>(); }

	void format(std::ostream& os) const;

	template <typename W2>
	W2 toScalar() const { W2 w; getValue(w); return w; }
	///@}

	///@name Bit field access operations
	///@{
	bool getBit(unsigned pos) const;
	void setBit(unsigned pos, bool bitVal);

	template < typename value_t >
	void getField(unsigned pos, unsigned wid, value_t& value) const;

	template < typename value_t >
	value_t getField(unsigned pos, unsigned wid) const;

	template < typename value_t >
	void setField(unsigned pos, unsigned wid, const value_t& value);
	///@}

protected:
	///@name Raster (array) operations
	///@{
	template < class AtomOp, typename AtomType >
	static void apply(AtomType* pTarget, size_t maxAtoms, const AtomType* pSource, size_t nAtoms);

	template < class AtomOp, typename AtomType >
	void apply(const AtomType* sourceArray, size_t arraySize);

	template < class AtomOp, typename AtomType >
	void apply(AtomType* targetArray, size_t arraySize) const;

	template < class AtomOp, typename AtomType >
	void apply(const AtomType& sourceAtom) { apply<AtomOp, AtomType>(&sourceAtom, 1u); }

	template < class AtomOp, typename AtomType >
	void apply(AtomType& targetAtom) const { apply<AtomOp, AtomType>(&targetAtom, 1u); }

	template < class AtomOp, typename SourceWord >
	void apply(const std::vector<SourceWord>& sourceVect);

	template < class AtomOp, typename TargetWord >
	void apply(std::vector<TargetWord>& targetVect) const;

	template < class AtomOp >
	void apply(const unsigned_Atom& sourceAtom);

	template < class AtomOp >
	void apply(unsigned_Atom& targetAtom) const;

	template < class AtomOp >
	void apply(const signed_Atom& sourceAtom);

	template < class AtomOp >
	void apply(signed_Atom& targetAtom) const;

	template < class AtomOp, typename A2 >
	void apply(const RtBitSet<A2>& sourceRtBitSet);
	///@}

private:
	size_t bitCount_;
	size_t atomCount_;
	Atom* p_buffer; ///< Actual buffer on which this bitset is working.
	Atom* i_buffer; ///< Internally allocated buffer.
	void* x_buffer; ///< Externally defined buffer.
};

/**
 * Run time bitset destructor: de-allocates locally allocated buffer (internal buffer).
 */
template < typename A >
RtBitSet<A>::~RtBitSet()
{
	if (i_buffer != 0) delete [] i_buffer;
}

template < typename A >
RtBitSet<A>& RtBitSet<A>::resetBuffer(size_t _bitCount, void* _x_buffer)
{
	bitCount_ = _bitCount;
	atomCount_ = bitCount_/atomBitCount() + (bitCount_ % atomBitCount() == 0 ? 0 : 1);

	x_buffer = _x_buffer; // Set external buffer (might be null)
	if (i_buffer != 0)
	{
		delete [] i_buffer; // Free possible former allocated buffer
		i_buffer = 0;
	}

	if (x_buffer == 0) // If external buffer is null
	{
		i_buffer = new Atom[atomCount_]; // Allocate internal buffer
		p_buffer = i_buffer; // Set actual BitSet buffer to allocated buffer
	}
	else
	{
		i_buffer = 0;
		p_buffer = reinterpret_cast<Atom*>(x_buffer); // Set actual BitSet buffer to external buffer
	}
	return *this;
}

/**
 * Fills bitset with specified word value.
 * @param initWord The value which must be copied to the bitset as many times as needed to fill all its words.
 * @return
 */
template < typename A >
RtBitSet<A>& RtBitSet<A>::fill(const Atom& initWord)
{
	for (size_t i=0u; i < atomCount(); ++i)
	{
		buffer()[i] = initWord;
	}
	normalize();
	return *this;
}

/**
 * Normalizes the RtBitSet format, namely by setting the sign bit properly.
 */
template < typename A >
RtBitSet<A>& RtBitSet<A>::normalize()
{
	const size_t bitsTail = bitCount() % atomBitCount();
	if (bitsTail > 0)
	{
		const size_t deadBits = atomBitCount() - bitsTail;
		atomLast() <<= deadBits;
		atomLast() >>= deadBits;
	}
	return *this;
}

template < typename A >
void RtBitSet<A>::bitNot()
{
	apply<AtomOperator::bitNot, Atom>(bufferPtr<Atom>(), atomCount(), bufferPtr<Atom>(), atomCount());
}

template < typename A >
void RtBitSet<A>::shiftLeft(size_t nBits)
{
	if (nBits >= bitCount())
	{
		setValue(0);
		return;
	}

	iter_Atom itBegin(itAtomBegin());

	size_t nAtomJumps = nBits / atomBitCount();
	size_t nBitJumpsLeft  =  nBits - atomBitCount() * nAtomJumps;
	size_t nBitJumpsRight =  atomBitCount() - nBitJumpsLeft;

	iter_Atom itTarget(itAtomLast());
	iter_Atom itSource(itTarget - nAtomJumps);

	while (itSource > itBegin)
	{
		*itTarget  = *itSource << nBitJumpsLeft;
		*itTarget |= *--itSource >> nBitJumpsRight;
		--itTarget;
	}

	*itTarget = *itSource << nBitJumpsLeft;
	while(--itTarget >= itBegin)
	{
		*itTarget = 0;
	}
	normalize();
}

template < typename A >
void RtBitSet<A>::shiftRight(size_t nBits)
{
	if (nBits >= bitCount())
	{
		setValue(isPositive() ? 0 : -1);
		return;
	}

	iter_Atom itLast(itAtomLast());

	size_t nAtomJumps = nBits / atomBitCount();
	size_t nBitJumpsRight =  nBits - atomBitCount() * nAtomJumps;
	size_t nBitJumpsLeft  =  atomBitCount() - nBitJumpsRight;

	iter_Atom itTarget(itAtomBegin());
	iter_Atom itSource (itTarget + nAtomJumps);

	while (itSource < itLast)
	{
		*itTarget  = static_cast<unsigned_Atom>(*itSource) >> nBitJumpsRight;
		*itTarget |= *++itSource << nBitJumpsLeft;
		++itTarget;
	}

	*itTarget  = *itSource >> nBitJumpsRight;
	Atom padding = *itLast < 0 ? ~0 : 0;
	while(++itTarget <= itLast)
	{
		*itTarget = padding;
	}
	normalize();
}

/**
 * Returns 0 if this bitset is numerically equal to a specified word, 1 if it is numerically greater, -1 otherwise.
 * @tparam W2 Type of the word this bitset is compared to.
 * @param w The word this bitset is compared to.
 */
template < typename A >
template < typename W2 >
int RtBitSet<A>::compare(const W2& w) const
{
	if (WordInfo<W2>::BIT_COUNT > bitCount())
	{
		typedef typename boost::mpl::if_c
		<
			boost::is_signed<A>::value,
			typename boost::make_signed <W2>::type,
			typename boost::make_unsigned <W2>::type
		>
		::type ThisWord;

		ThisWord thisWord;
		getValue(thisWord);
		     if (thisWord > w) return  1;
		else if (thisWord < w) return -1;
		                       return  0;
	}
	else
	{
		typename make_RtBitSet<W2>::type bitSet(bitCount(), w);
		return compare(bitSet);
	}
}

/**
 * Returns 0 if this bitset is numerically equal to another specified bitset, 1 if it is numerically greater, -1 otherwise.
 * @tparam W2 Signedness of the bitset this one is compared to.
 * @param w The bitset this one is compared to.
 */
template < typename A >
template < typename A2 >
int RtBitSet<A>::compare(const RtBitSet<A2>& bitSet) const
{
	if (bitCount() == bitSet.bitCount())
	{
		iter_Atom itBegin = itAtomBegin();
		iter_Atom it = itAtomEnd();
		WordIterator<typename RtBitSet<A2>::Atom> jt = bitSet.itAtomEnd();
		do
		{
			--it; --jt;
			     if (*it > *jt) return  1;
			else if (*it < *jt) return -1;

		} while (it != itBegin);
		return 0;
	}
	else
	{
		if (bitSet.bitCount() > bitCount())
		{
			RtBitSet<A> largerThis(bitSet.bitCount(), *this);
			return largerThis.compare(bitSet); // recursive call
		}
		else
		{
			RtBitSet<A2> largerThat(bitCount(), bitSet);
			return compare(largerThat); // recursive call
		}
	}
}

/**
 * This static method executes in a loop over the source and target array words the same elementary operation (e.g. copy
 * operation). In the case where target array is of a greater size than source array, excess target words are "wiped"
 * according to the signedness of the array word type (if AtomType is signed, the sign bit is extended).
 * @tparam AtomType The type of the elementary element of the arrays.
 * @tparam AtomOp The elementary operation applied (as functor).
 * @param pTarget Pointer to the head of the target array.
 * @param maxAtoms Size of the target array.
 * @param pSource Pointer to the head of the source array.
 * @param nAtoms Size of the source array.
 */
template < typename A >
template < class AtomOp, typename AtomType >
void RtBitSet<A>::apply(AtomType* pTarget, size_t maxAtoms, const AtomType* pSource, size_t nAtoms)
{
	AtomType wipeValue = utl::isPositive(*(WordIterator<AtomType>::getEnd(pSource, nAtoms) - 1)) ? 0 : -1;
	WordIterator<AtomType> targetIt = WordIterator<AtomType>::getBegin(pTarget, maxAtoms);
	WordIterator<AtomType> sourceIt = WordIterator<AtomType>::getBegin(pSource, nAtoms);

	std::ptrdiff_t delta = nAtoms - maxAtoms;
	size_t nOps = nAtoms;
	if (delta >= 0)
	{
		nOps = maxAtoms;
	}

	while (nOps-- > 0)
	{
		AtomOp()(*targetIt, *sourceIt);
		++targetIt;
		++sourceIt;
	}

	while (delta++ < 0)
	{
		AtomOp()(*targetIt, wipeValue);
		++targetIt;
	}
}

/**
 * Specialized apply method from array of words to bitset.
 * @tparam AtomOp Elementary raster operation (e.g. copy, bitwise or, etc.)
 * @tparam AtomType Elementary word type of the array.
 * @param sourceArray The atom array whose value is operated with this bitset.
 * @param arraySize The size of the array.
 */
template < typename A >
template < class AtomOp, typename AtomType >
void RtBitSet<A>::apply(const AtomType* sourceArray, size_t arraySize)
{
	// Define Byte with same sign as source, i.e. AtomType template parameter (here)
	typedef typename boost::mpl::if_c
		<
			boost::is_signed<AtomType>::value,
			signed char,
			unsigned char
		>
		::type SourceByte;

	apply< AtomOp, SourceByte >
	(
		bufferPtr<SourceByte>(), byteCount(), reinterpret_cast<const SourceByte*>(sourceArray), arraySize*sizeof(AtomType)
	);
	normalize();
}

/**
 * Specialized apply method from bitset to array of words.
 * @param targetArray The atom array whose value is operated with this bitset.
 */
template < typename A >
template < class AtomOp, typename AtomType >
void RtBitSet<A>::apply(AtomType* targetArray, size_t arraySize) const
{
	apply<AtomOp, Byte>
	(
		reinterpret_cast<Byte*>(targetArray), arraySize*sizeof(AtomType), bufferPtr<Byte>(), byteCount()
	);
}

/**
 * This method executes the specified elementary operation between a source vector and this bitset (e.g. copy
 * operation). In this case this bitset is modified but not the source vector.
 * @tparam SourceWord The type of the source vector's elementary word.
 * @tparam AtomOp The elementary operation applied (as functor).
 */
template < typename A >
template < class AtomOp, typename SourceWord >
void RtBitSet<A>::apply(const std::vector<SourceWord>& sourceVect)
{
	typedef typename boost::mpl::if_c
		<
			boost::is_signed<SourceWord>::value,
			signed char,
			unsigned char
		>
		::type SourceByte;

	if (sourceVect.size() > 0)
	{
		apply<AtomOp, SourceByte>
		(
			bufferPtr<SourceByte>(), byteCount(), reinterpret_cast<const SourceByte*>(&sourceVect[0]), sourceVect.size()*sizeof(SourceWord)
		);
		normalize();
	}
}

/**
 * This method executes the specified elementary operation between this bitset a target vector (e.g. copy operation). In
 * this case the target vector is modified but not this bitset. The size of the target vector is adjusted so that it
 * can contain the whole bitset.
 * @tparam TargetWord The type of the target vector's elementary word.
 * @tparam AtomOp The elementary operation applied (as functor).
 */
template < typename A >
template < class AtomOp, typename TargetWord >
void RtBitSet<A>::apply(std::vector<TargetWord>& targetVect) const
{
	size_t lcmBytes = boost::math::lcm(sizeof(TargetWord), byteCount());
	targetVect.resize(lcmBytes / sizeof(TargetWord));
	for (size_t i = 0; i < targetVect.size(); ++i)
	{
		targetVect[i] = 0;
	}
	size_t maxTargetAtoms = lcmBytes / atomByteCount();
	if (targetVect.size() > 0)
	{
		apply<AtomOp, Atom>
		(
			reinterpret_cast<Atom*>(&targetVect[0]), maxTargetAtoms, bufferPtr<Atom>(), atomCount()
		);
	}
}

/**
 * Specialized apply method from word to bitset method for words of type Atom in its unsigned version.
 * This specialization optimizes the performance of the raster operation by using words bigger than bytes.
 * @param sourceWord The word whose value is operated with this bitset.
 */
template < typename A >
template < class AtomOp >
void RtBitSet<A>::apply(const unsigned_Atom& sourceWord)
{
	apply<AtomOp, unsigned_Atom>(bufferPtr<unsigned_Atom>(), atomCount(), &sourceWord, 1u);
	normalize();
}

/**
 * Specialized apply method from bitset to word for words of type Atom in its unsigned version.
 * This specialization optimizes the performance of the raster operation by using words bigger than bytes.
 * @param targetWord The word whose value is operated with this bitset.
 */
template < typename A >
template < class AtomOp >
void RtBitSet<A>::apply(unsigned_Atom& targetWord) const
{
	apply<AtomOp, unsigned_Atom>(&targetWord, 1u, bufferPtr<unsigned_Atom>(), atomCount());
}

/**
 * Specialized apply method from word to bitset method for words of type Atom in its signed version.
 * This specialization optimizes the performance of the raster operation by using words bigger than bytes.
 * @param sourceWord The word whose value is operated with this bitset.
 */
template < typename A >
template < class AtomOp >
void RtBitSet<A>::apply(const signed_Atom& sourceWord)
{
	apply<AtomOp, signed_Atom>(bufferPtr<signed_Atom>(), atomCount(), &sourceWord, 1u);
	normalize();
}

/**
 * Specialized apply method from bitset to word for words of type Atom in its signed version.
 * This specialization optimizes the performance of the raster operation by using words bigger than bytes.
 * @param sourceWord The word whose value is operated with this bitset.
 */
template < typename A >
template < class AtomOp >
void RtBitSet<A>::apply(signed_Atom& targetWord) const
{
	apply<AtomOp, signed_Atom>(&targetWord, 1u, bufferPtr<signed_Atom>(), atomCount());
}

/**
 * Specialized apply method from bitset to bitset.
 * This specialization optimizes the performance of the raster operation by using words bigger than bytes.
 * @param sourceRtBitSet The bitset whose value is operated with this bitset.
 * @param arraySize The size of the array.
 */
template < typename A >
template < class AtomOp, typename A2 >
void RtBitSet<A>::apply(const RtBitSet<A2>& sourceRtBitSet)
{
	typedef typename RtBitSet<A2>::Atom AtomType;
	apply<AtomOp, AtomType>(sourceRtBitSet.template bufferPtr<AtomType>(),  sourceRtBitSet.atomCount());
	normalize();
}

template < typename A >
void RtBitSet<A>::format(std::ostream& os) const
{
	os << (boost::is_signed<A>::value ? "(s)" : "(u)");
	iter_Atom it(itAtomEnd());
	iter_Atom itBegin(itAtomBegin());
	os << std::hex << std::setfill('0');
	size_t iBit = atomCount() * sizeof(Atom) * CHAR_BIT - 1;
	while (it != itBegin)
	{
		--it;
		for (int i = atomBitCount()-1; i >= 0; --i)
		{
			if (iBit < bitCount())
			{
				os << ((*it >> i) & 1);
			}
			else
			{
				os << '.';
			}
			if (iBit != 0)
			{
				if (iBit % atomBitCount() == 0) os << '|';
				else if(iBit % 8 == 0) os << '\'';
				else if(iBit % 4 == 0) os << ' ';
			}
			--iBit;
		}
	}
	os << std::dec;
}

template < typename A >
bool RtBitSet<A>::getBit(unsigned pos) const
{
	return BitFieldHelper<Atom>::getBit(atomAt(pos / atomBitCount()), pos % atomBitCount());

}

template < typename A >
void RtBitSet<A>::setBit(unsigned pos, bool bitVal)
{
	BitFieldHelper<Atom>::setBit(atomAt(pos / atomBitCount()), pos % atomBitCount(), bitVal);
	normalize();
}

template < typename A >
template < typename value_t >
void RtBitSet<A>::getField(unsigned pos, unsigned wid, value_t& value) const
{
	std::size_t iAtom = pos / atomBitCount();
	iter_Atom it = itAtomBegin() + iAtom;

	unsigned posInAtom = pos - iAtom * atomBitCount();
	unsigned widInAtom = std::min(wid, atomBitCount() - posInAtom);
	if (widInAtom == wid) // if this is the last part, extraction must have same signedness as value
	{
		value = BitFieldHelper<Atom>::getField(*it++, posInAtom, widInAtom);
	}
	else  // if this is NOT the last part, extraction must be unsigned (to avoid risking sign extension if value_t is signed)
	{
		typename boost::make_unsigned<value_t>::type xtract;
		BitFieldHelper<Atom>::getField(*it++, posInAtom, widInAtom, xtract);
		value = xtract;
		for (unsigned widSum = widInAtom; widSum < wid; widSum += widInAtom)
		{
			widInAtom = std::min(wid - widSum, atomBitCount());
			if (widSum + widInAtom < wid)
			{
				typename boost::make_unsigned<value_t>::type xtract;
				BitFieldHelper<Atom>::getField(*it++, 0, widInAtom, xtract);
				value |= xtract << widSum;
			}
			else // this is the last part => extraction must have same signedness as value
			{
				value_t xtract;
				BitFieldHelper<Atom>::getField(*it++, 0, widInAtom, xtract);
				value |= xtract << widSum;
			}
		}
	}
}

template < typename A >
template < typename value_t >
value_t RtBitSet<A>::getField(unsigned pos, unsigned wid) const
{
	value_t val;
	getField(pos, wid, val);
	return val;
}

template < typename A >
template < typename value_t >
void RtBitSet<A>::setField(unsigned pos, unsigned wid, const value_t& value)
{
	std::size_t iAtom = pos / atomBitCount();
	iter_Atom it = itAtomBegin() + iAtom;

	unsigned posInAtom = pos - iAtom * atomBitCount();
	unsigned widInAtom = std::min(wid, atomBitCount() - posInAtom);

	BitFieldHelper<Atom>::setField(*it++, posInAtom, widInAtom, value);
	for (unsigned widSum = widInAtom; widSum < wid; widSum += widInAtom)
	{
		widInAtom = std::min(wid - widSum, static_cast<unsigned>(atomBitCount()));
		BitFieldHelper<Atom>::setField(*it++, 0, widInAtom, value >> widSum);
	}
	normalize();
}

template < typename A >
std::ostream& operator << (std::ostream& os, const ::utl::RtBitSet<A>& bitSet)
{
	bitSet.format(os);
	return os;
}

template < typename WordType, typename A >
WordType operator | (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	utl::RtBitSet<WordType> bitSet(wordLeft);
	bitSet |= bitSetRight;
	WordType w;
	bitSet.getValue(w);
	return w;
}

template < typename AL, typename AR >
::utl::RtBitSet<AL> operator | (const ::utl::RtBitSet<AL>& bitSetLeft, const ::utl::RtBitSet<AR>& bitSetRight)
{
	::utl::RtBitSet<AL> bitSet(bitSetLeft);
	bitSet |= bitSetRight;
	return bitSet;
}

template < typename WordType, typename A >
WordType operator & (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	typename ::utl::make_RtBitSet<WordType>::type bitSet(wordLeft);
	bitSet &= bitSetRight;
	WordType w;
	bitSet.getValue(w);
	return w;
}

template < typename AL, typename AR >
::utl::RtBitSet<AL> operator & (const ::utl::RtBitSet<AL>& bitSetLeft, const ::utl::RtBitSet<AR>& bitSetRight)
{
	::utl::RtBitSet<AL> bitSet(bitSetLeft);
	bitSet ^= bitSetRight;
	return bitSet;
}

template < typename WordType, typename A >
WordType operator ^ (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	typename ::utl::make_RtBitSet<WordType>::type bitSet(wordLeft);
	bitSet ^= bitSetRight;
	WordType w;
	bitSet.getValue(w);
	return w;
}

template < typename AL, typename AR >
::utl::RtBitSet<AL> operator ^ (const ::utl::RtBitSet<AL>& bitSetLeft, const ::utl::RtBitSet<AR>& bitSetRight)
{
	::utl::RtBitSet<AL> bitSet(bitSetLeft);
	bitSet ^= bitSetRight;
	return bitSet;
}

template < typename WordType, typename A >
bool operator > (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	return bitSetRight < wordLeft;
}

template < typename WordType, typename A >
bool operator < (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	return bitSetRight > wordLeft;
}

template < typename WordType, typename A >
bool operator >= (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	return bitSetRight <= wordLeft;
}

template < typename WordType, typename A >
bool operator <= (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	return bitSetRight >= wordLeft;
}

template < typename WordType, typename A >
bool operator == (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	return bitSetRight == wordLeft;
}

template < typename WordType, typename A >
bool operator != (const WordType& wordLeft, const ::utl::RtBitSet<A>& bitSetRight)
{
	return bitSetRight != wordLeft;
}

} // namespace utl

#ifdef UNDEF_VX_MINMAX
// Defining back VxWorks min/max macros
#define max(x, y)	(((x) < (y)) ? (y) : (x))
#define min(x, y)	(((x) < (y)) ? (x) : (y))
#endif

#endif /* utl_RtBitSet_h_INCLUDED */
