// =============================================================================
// BitFieldHelper.hpp
// -----------------------------------------------------------------------------
// class utl::BitFieldHelper
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef utl_BitFieldHelper_h_INCLUDED
#define utl_BitFieldHelper_h_INCLUDED

#include "utl/Logging.h"
#include <boost/type_traits/make_signed.hpp>
#include <boost/type_traits/make_unsigned.hpp>
#include <boost/type_traits/is_signed.hpp>
#include <boost/type_traits/remove_cv.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>
#include <climits>

namespace utl
{

/**
 * This class provides template methods to manipulate bit-fields easily and efficiently.
 */
template<class WordType>
class BitFieldHelper
{
public:
	typedef typename boost::remove_cv<WordType>::type WordBasicType;
	typedef typename boost::make_signed<WordBasicType>::type WordSignedType;
	typedef typename boost::make_unsigned<WordBasicType>::type WordUnsignedType;

	/**
	 * Compile-time mapping that returns the signed or unsigned version of WordType
	 * depending on the sign of the specified template parameter.
	 */
	template<typename ValueType>
	struct WordTypeMatchingSignOf
		: boost::mpl::if_c<
			boost::is_signed<ValueType>::value,
			WordSignedType,
			WordUnsignedType>
	{
	};

	static const WordBasicType ZERO = 0;
	static const WordBasicType ONE = 1;
	static const WordBasicType ALL = ~0;
	static const unsigned SIZEOF = sizeof(WordBasicType);
	static const unsigned BITCOUNT = SIZEOF * CHAR_BIT;

	///@name Compile-time bit-field manipulation.
	///@{
	template<WordUnsignedType pos>
	static bool getBit(WordBasicType var);

	template<WordUnsignedType pos>
	static void clearBit(WordType& var);

	template<WordUnsignedType pos>
	static void setBit(WordType& var);

	template<WordUnsignedType pos>
	static void setBit(WordType& var, bool askSet);

	template<WordUnsignedType pos, WordUnsignedType width>
	static WordBasicType getField(WordBasicType var);

	template<WordUnsignedType pos, WordUnsignedType width, typename FieldType>
	static FieldType getField(WordBasicType var);

	template<WordUnsignedType pos, WordUnsignedType width, typename FieldType>
	static void getField(WordBasicType var, FieldType& result);

	template<WordUnsignedType pos, WordUnsignedType width>
	static void setField(WordType& var, WordBasicType value);

	template<WordUnsignedType pos, WordUnsignedType width>
	static WordBasicType getReadMask();

	template<WordUnsignedType pos, WordUnsignedType width>
	static WordBasicType getClearMask();
	///@}

	///@name Run-time bit-field manipulation.
	///@{
	static bool getBit(WordBasicType var, WordUnsignedType pos);

	static void setBit(WordType& var, WordUnsignedType pos, bool askSet);

	static WordBasicType getField(WordBasicType var, WordUnsignedType pos, WordUnsignedType width);

	template<typename FieldType>
	static FieldType getField(WordBasicType var, WordUnsignedType pos, WordUnsignedType width);

	template<typename FieldType>
	static void getField(WordBasicType var, WordUnsignedType pos, WordUnsignedType width, FieldType& result);

	static void setField(WordType& var, WordUnsignedType pos, WordUnsignedType width, WordBasicType value);

	static WordBasicType getReadMask(WordUnsignedType pos, WordUnsignedType width);

	static WordBasicType getClearMask(WordUnsignedType pos, WordUnsignedType width);
	///@}

private:
	BitFieldHelper() {}
	~BitFieldHelper() {}
};

// =============================================================================
// Wrappers functions for some built-in bitwise operators.
// They allow stricter type coherency (no automatic integer promotion),
// and allow alternative implementations.
// Overloading operators for primitive types would propagate out of control, so
// this was not an option.
// -----------------------------------------------------------------------------

/**
 * Performs a bitwise NOT, or bit-toggling.
 * @param value
 * @return
 */
template<typename T>
T toggle(T value)
{
	return ~value;
}

/**
 * Computes an integer exponential of base 2 (which essentially have the nth bit at 1 and the rest at 0).
 * @param n
 * @return
 */
template<typename T>
T exp2(T n)
{
	return 1 << n;
}

#if __cplusplus >= 201103L
/**
 * Binary left shift for unsigned values (C++11 compliant).
 * @param value the value to shift
 * @param nbits the number of bits to shift
 * @return Returns value left-shifted by nbits.
 */
template<typename T>
T lhs(T value, typename boost::enable_if<boost::is_unsigned<T>, T>::type nbits)
{
	return value << nbits;
}

/**
 * Binary left shift for signed values (C++11 compliant).
 * This is a workaround for C++11 considering the left shift of a signed number as
 * an undefined behaviour.
 * @param value the value to shift
 * @param nbits the number of bits to shift
 * @return Returns value left-shifted by nbits.
 */
template<typename T>
T lhs(T value, typename boost::enable_if<boost::is_signed<T>, typename boost::make_unsigned<T>::type>::type nbits)
{
	return value * exp2(nbits);
}
#else
/**
 * Binary left shift.
 * @param value the value to shift
 * @param nbits the number of bits to shift
 * @return
 */
template<typename T>
T lhs(T value, typename boost::make_unsigned<T>::type nbits)
{
	return value << nbits;
}
#endif

/**
 * Binary right shift.
 * @param value the value to shift
 * @param nbits the number of bits to shift
 * @return Returns value right-shifted by nbits.
 */
template<typename T>
T rhs(T value, typename boost::make_unsigned<T>::type nbits)
{
	return value >> nbits;
}

// =============================================================================
// utl::BitFieldHelper inline methods body
// -----------------------------------------------------------------------------

/**
 * Returns the value of a bit of an integer variable (compile-time version).
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos>
inline bool BitFieldHelper<WordType>::getBit(WordBasicType var)
{
	return (var & exp2(pos)) != 0;
}

/**
 * Returns the value of a bit of an integer variable.
 */
template<class WordType>
inline bool BitFieldHelper<WordType>::getBit(WordBasicType var, WordUnsignedType pos)
{
	return (var & exp2(pos)) != 0;
}

/**
 * Clears the bit at the specified position.
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos>
inline void BitFieldHelper<WordType>::clearBit(WordType& var)
{
	var &= toggle(exp2(pos));
}

/**
 * Sets the bit at the specified position.
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos>
inline void BitFieldHelper<WordType>::setBit(WordType& var)
{
	var |= exp2(pos);
}

/**
 * Sets or clears a bit of a variable (compile-time version).
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos>
inline void BitFieldHelper<WordType>::setBit(WordType& var, bool askSet)
{
	if(askSet) setBit<pos>(var);
	else       clearBit<pos>(var);
}

/**
 * Sets or clears a bit of a variable.
 */
template<class WordType>
inline void BitFieldHelper<WordType>::setBit(WordType& var, WordUnsignedType pos, bool askSet)
{
	if(askSet) var |= exp2(pos);
	else       var &= toggle(exp2(pos));
}

/**
 * Returns the value of a bit-field, with the same type as the input (compile-time version).
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos, typename BitFieldHelper<WordType>::WordUnsignedType width>
inline typename BitFieldHelper<WordType>::WordBasicType
BitFieldHelper<WordType>::getField(WordBasicType var)
{
	WordBasicType result;
	getField<pos, width>(var, result);
	return result;
}
/**
 * Returns the value of a bit-field, with the same type as the input.
 */
template<class WordType>
inline typename BitFieldHelper<WordType>::WordBasicType
BitFieldHelper<WordType>::getField(WordBasicType var, WordUnsignedType pos, WordUnsignedType width)
{
	WordBasicType result;
	getField(var, pos, width, result);
	return result;
}

/**
 * Returns the value of a bit-field, allowing the field to have a different type (compile-time version).
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos, typename BitFieldHelper<WordType>::WordUnsignedType width, typename FieldType>
inline FieldType BitFieldHelper<WordType>::getField(WordBasicType var)
{
	FieldType result;
	getField<pos, width>(var, result);
	return result;
}

/**
 * Returns the value of a bit-field, allowing the field to have a different type.
 */
template<class WordType>
template<typename FieldType>
inline FieldType BitFieldHelper<WordType>::getField(WordBasicType var, WordUnsignedType pos, WordUnsignedType width)
{
	FieldType result;
	getField(var, pos, width, result);
	return result;
}

/**
 * Retrieves the value of a bit-field (compile-time version).
 * The field does not need to have the same type (or same signed-ness) as the input.
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos, typename BitFieldHelper<WordType>::WordUnsignedType width, typename FieldType>
inline void BitFieldHelper<WordType>::getField(WordBasicType var, FieldType& result)
{
	typedef typename WordTypeMatchingSignOf<FieldType>::type SignCorrectWordType;

	result = static_cast<FieldType>(rhs(static_cast<SignCorrectWordType>(lhs(var, BITCOUNT - width - pos)), BITCOUNT - width));
}

/**
 * Retrieves the value of a bit-field.
 * The field does not need to have the same type (or same signed-ness) as the input.
 */
template<class WordType>
template<typename FieldType>
inline void BitFieldHelper<WordType>::getField(WordBasicType var, WordUnsignedType pos, WordUnsignedType width, FieldType& result)
{
	typedef typename WordTypeMatchingSignOf<FieldType>::type SignCorrectWordType;

	result = static_cast<FieldType>(rhs(static_cast<SignCorrectWordType>(lhs(var, BITCOUNT - width - pos)), BITCOUNT - width));
}

/**
 * Sets the value of a bit-field of a variable (compile-time version).
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos, typename BitFieldHelper<WordType>::WordUnsignedType width>
inline void BitFieldHelper<WordType>::setField(WordType& var, WordBasicType value)
{
	var = (var & getClearMask<pos, width>()) | (lhs(value, pos) & getReadMask<pos, width>());
}

/**
 * Sets the value of a bit-field of a variable.
 */
template<class WordType>
inline void BitFieldHelper<WordType>::setField(WordType& var, WordUnsignedType pos, WordUnsignedType width, WordBasicType value)
{
	var = (var & getClearMask(pos, width)) | (lhs(value, pos) & getReadMask(pos, width));
}

/**
 * Returns the mask used to isolate the bit-field range, in order to read it (compile-time version).
 * Protection against the machine-dependent case (and seen at run-time)
 * for which width == size of the word type
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos, typename BitFieldHelper<WordType>::WordUnsignedType width>
inline typename BitFieldHelper<WordType>::WordBasicType
BitFieldHelper<WordType>::getReadMask()
{
	return (width >= BITCOUNT) ? ALL : lhs(toggle(lhs(ALL, width)), pos);
}

/**
 * Returns the mask used to isolate the bit-field range, in order to read it.
 * Protection against the machine-dependent case (and seen at run-time)
 * for which width == size of the word type
 */
template<class WordType>
inline typename BitFieldHelper<WordType>::WordBasicType
BitFieldHelper<WordType>::getReadMask(WordUnsignedType pos, WordUnsignedType width)
{
	return (width >= BITCOUNT) ? ALL : lhs(toggle(lhs(ALL, width)), pos);
}

/**
 * Returns the mask used to exclude the bit-field range, in order to clear it (compile-time version).
 */
template<class WordType>
template<typename BitFieldHelper<WordType>::WordUnsignedType pos, typename BitFieldHelper<WordType>::WordUnsignedType width>
inline typename BitFieldHelper<WordType>::WordBasicType
BitFieldHelper<WordType>::getClearMask()
{
	return toggle(getReadMask<pos, width>());
}

/**
 * Returns the mask used to exclude the bit-field range, in order to clear it.
 */
template<class WordType>
inline typename BitFieldHelper<WordType>::WordBasicType
BitFieldHelper<WordType>::getClearMask(WordUnsignedType pos, WordUnsignedType width)
{
	return toggle(getReadMask(pos, width));
}

}

#endif // utl_BitFieldHelper_h_INCLUDED
