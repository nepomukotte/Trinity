/**
 * @file BitFieldHelperVector.hpp
 * @created Nov 12, 2012
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


#ifndef utl_BitFieldHelperVector_hpp_INCLUDED
#define utl_BitFieldHelperVector_hpp_INCLUDED

#include "config.h"
#include "utl/BitFieldHelperBitSet.hpp"
#include <boost/type_traits/make_signed.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/conditional.hpp>
#include <vector>
#include <algorithm>
#include <climits>
#include <iostream>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace utl
{
/**
 * Template specialization of BitFieldHelper for vector objects.
 */
template < class W >
class BitFieldHelper < std::vector<W> >
{
public:
	typedef W WordType;
	typedef std::vector < WordType > VectorType;
	typedef WordIterator<WordType> WordIter;
	enum { WORD_BITCOUNT = sizeof(WordType) * CHAR_BIT };

public:
	/// @name Compile-time bit-field manipulation.
	//@{
	template < unsigned pos >
	static bool getBit(const VectorType& vect);

	template < unsigned pos, unsigned width, typename ValueType >
	static void getField(const VectorType& vect, ValueType& value);

	template < unsigned pos, unsigned width, typename ValueType >
	static ValueType getField(const VectorType& vect);

	template < unsigned pos >
	static void setBit(VectorType& vect, bool askSet);

	template < unsigned pos, unsigned width, typename ValueType >
	static void setField(VectorType& vect, const ValueType& value);

	/// @name Run-time bit-field manipulation.
	//@{
	static bool getBit(const VectorType& vect, unsigned pos);
	static void setBit(VectorType& vect, unsigned pos, bool askSet);

	template < typename ValueType >
	static void getField(const VectorType& vect, unsigned pos, unsigned width, ValueType& value);

	template < typename ValueType >
	static ValueType getField(const VectorType& vect, unsigned pos, unsigned width);

	template < typename ValueType >
	static void setField(VectorType& vect, unsigned pos, unsigned width, const ValueType& value);

	static void shiftLeft(VectorType& vect, size_t nBits);
	static void shiftRight(VectorType& vect, size_t nBits);
	//@}

private:
	BitFieldHelper() {}
	~BitFieldHelper() {}
};

// =============================================================================
// utl::BitFieldHelper inline methods body (compile-time versions)
// -----------------------------------------------------------------------------

/**
 * Returns the value of a bit in a vect.
 */
template < class W >
template < unsigned pos >
inline bool BitFieldHelper< std::vector<W> >::getBit(const VectorType& vect)
{
	return BitFieldHelper< WordType >::getBit< pos % WORD_BITCOUNT >(WordIter::getBegin(vect)[pos/WORD_BITCOUNT]);
}

/**
 * Sets out parameter to the value of a bit field in a vect.
 */
template < class W >
template<unsigned pos, unsigned width, typename ValueType>
inline void BitFieldHelper< std::vector<W> >::getField(const VectorType& vect, ValueType& value)
{
	getField(vect, pos, width, value);
}

/**
 * Returns the value of a bit field in a vect.
 */
template < class W >
template<unsigned pos, unsigned width, typename ValueType>
inline ValueType BitFieldHelper< std::vector<W> >::getField(const VectorType& vect)
{
	ValueType value;
	getField(vect, pos, width, value);
	return value;
}

/**
 * Sets or clears a bit of a variable.
 */
template < class W >
template<unsigned pos>
inline void BitFieldHelper< std::vector<W> >::setBit(VectorType& vect, bool askSet)
{
	BitFieldHelper< WordType >::setBit< pos % WORD_BITCOUNT >(WordIter::getBegin(vect)[pos/WORD_BITCOUNT], askSet);
}

/**
 * Sets the value of a bit field of a variable.
 */
template < class W >
template<unsigned pos, unsigned width, typename ValueType>
inline void BitFieldHelper< std::vector<W> >::setField(VectorType& vect, const ValueType& value)
{
	setField(vect, pos, width, value);
}

// =============================================================================
// utl::BitFieldHelper inline methods body (run-time versions)
// -----------------------------------------------------------------------------

/**
 * Returns the value of a bit  in a vect (whose position is not known at compile-time).
 */
template < class W >
inline bool BitFieldHelper< std::vector<W> >::getBit(const VectorType& vect, unsigned pos)
{
	return BitFieldHelper< WordType >::getBit(WordIter::getBegin(vect)[pos/WORD_BITCOUNT], pos % WORD_BITCOUNT);
}

/**
 * Sets or clears a bit of a variable (whose position is not known at compile-time).
 */
template < class W >
inline void BitFieldHelper< std::vector<W> >::setBit(VectorType& vect, unsigned pos, bool askSet)
{
	BitFieldHelper< WordType >::setBit(WordIter::getBegin(vect)[pos/WORD_BITCOUNT], pos % WORD_BITCOUNT, askSet);
}

/**
 * Sets out parameter to the value of a bit field in a vect.
 */
template < class W >
template<typename ValueType>
inline void BitFieldHelper< std::vector<W> >::getField(const VectorType& vect, unsigned pos, unsigned width, ValueType& value)
{
	typename ::utl::make_BitSet<ValueType>::type bitSet(vect);
	bitSet.getField(pos, width, value);
}

/**
 * Returns the value of a bit field in a vect.
 */
template < class W >
template<typename ValueType>
inline ValueType BitFieldHelper< std::vector<W> >::getField(const VectorType& vect, unsigned pos, unsigned width)
{
	ValueType value;
	getField(vect, pos, width, value);
	return value;
}

/**
 * Sets the value of a bit field of a variable.
 */
template < class W >
template<typename ValueType>
inline void BitFieldHelper< std::vector<W> >::setField(VectorType& vect, unsigned pos, unsigned width, const ValueType& value)
{
	typename ::utl::make_BitSet<ValueType>::type bitSet(vect);
	bitSet.setField(pos, width, value);
	bitSet.getValue(vect);
}

template < class W >
void BitFieldHelper< std::vector<W> >::shiftLeft(VectorType& vect, size_t nBits)
{
	typename ::utl::make_BitSet<W>::type bitSet(vect);
	bitSet.shiftLeft(nBits);
	bitSet.getValue(vect);
}


template < class W >
inline void BitFieldHelper< std::vector<W> >::shiftRight(VectorType& vect, size_t nBits)
{
	typename ::utl::make_BitSet<W>::type bitSet(vect);
	bitSet.shiftRight(nBits);
	bitSet.getValue(vect);
}

template < class W >
std::vector<W>& operator >>= (std::vector<W>& vect, size_t nBits)
{
	::utl::BitFieldHelper< std::vector<W> >::shiftRight(vect, nBits);
	return vect;
}

template < class W >
std::vector<W>& operator <<= (std::vector<W>& vect, size_t nBits)
{
	::utl::BitFieldHelper< std::vector<W> >::shiftLeft(vect, nBits);
	return vect;
}

template < class W >
std::vector<W> operator >> (const std::vector<W>& vect, size_t nBits)
{
	std::vector<W> outVect(vect);
	outVect >>= nBits;
	return outVect;
}

template < class W >
std::vector<W> operator << (const std::vector<W>& vect, size_t nBits)
{
	std::vector<W> outVect(vect);
	outVect <<= nBits;
	return outVect;
}

} // namespace utl


#endif /* utl_BitFieldHelperVector_hpp_INCLUDED */
