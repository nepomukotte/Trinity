// =============================================================================
// Enum.h
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#ifndef CCfg_Enum_h_INCLUDED
#define CCfg_Enum_h_INCLUDED

#include "Range.h"
#include "Value.h"
#include "Exception.h"
#include <set>
#include <ios>
#include <cassert>

namespace CCfg
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class defines an enumeration range.
 * It is a list of values of several type (see CCfg::Value) that contains all
 * valid values for a configuration attribute.
 *
 * @author fchateau
 */
class Enum : public Range
{
public:
	typedef std::set<Value> Set;

	// inner classes
	class IncompatibleType;
	template<typename T> class Inserter;
	template<typename T> struct LinearList;

	// friends
	template<typename T> friend class Inserter;
	friend std::ostream& operator<<(std::ostream& stream, const Enum& enumeration);

	// static constants
	static const char OPEN_CHAR = '{';
	static const char CLOSE_CHAR = '}';
	static const char SEPARATOR = ',';
	static const char RESERVED_CHARACTERS[]; // concatenation of previous characters

	// static methods
	template<typename T> static LinearList<T> linear(T first, T last, T step = 1);




	// @{ construction/destruction
	Enum();
	Enum(const Enum& right);
	virtual ~Enum();
	Enum& operator=(const Enum& right);
	// @}

	// @{ base class overriden methods
	virtual bool contains(const Value& value) const;
	virtual bool equals(const Range& rvalue) const;
	virtual RangePtr clone() const;
	virtual std::string toString() const;
	virtual void fromString(const std::string& str);
	virtual void write(std::ostream& stream) const;
	// @}

	// @{ accessors
	const Set& getValues() const;
	// @}

	// @{ insertion methods
	template<typename T> void insert(T val);
	template<typename T> void insert(const LinearList<T>& lin);
	// @}

	// @{  insertion operators
	template<typename T> Inserter<T> operator=(T val);
	template<typename T> Inserter<T> operator+=(T val);
	template<typename T> Inserter<T> operator=(const LinearList<T>& lin);
	template<typename T> Inserter<T> operator+=(const LinearList<T>& lin);
	// @}

private:
	template<typename T> void unchecked_insert(T val);
	template<typename T> void unchecked_insert(const LinearList<T>& lin);
	void checkTypeHomogeneity(const Value& v);

	Set values;	///< the list of values in the enumeration.
};

std::ostream& operator<<(std::ostream& stream, const Enum& enumeration);


/**
 * This exception reports that a method attempted to insert a Value that would
 * have made the Enum heterogenous (have more than a unique type of values).
 *
 * @author fchateau
 */
class Enum::IncompatibleType : public Exception
{
	Enum reporterEnum;	///< a copy of the enum that reported the error.
	Value invalidValue;	///< the value that could not be inserted in the enum because of its type.
public:
	// @{ construction/destruction
	IncompatibleType(const Enum& reporterEnum, const Value& invalidValue);
	IncompatibleType(const IncompatibleType& rvalue);
	virtual ~IncompatibleType() throw();
	// @}

	// @{ accessors
	const Enum& getReporterEnum() const;
	ValueTypeId getEnumType() const;
	const Value& getInvalidValue() const;
	// @}

private:
	std::string formatMessage() const;
	IncompatibleType& operator=(const IncompatibleType& right);
};

/**
 * This is a temporary proxy object allowing to initialize a CCfg::Enum from a
 * list of values separated by commas.
 * \verbatim
 * Enum e;
 * e = 5, 10, 15, 20;
 * \endverbatim
 */
template<typename T>
class Enum::Inserter
{
public:
	Enum& enumeration;

	Inserter(Enum& enumeration);

	Inserter<T>& insert(T val);
	Inserter<T>& insert(const LinearList<T>& lin);

	Inserter<T>& operator,(T val);
	Inserter<T>& operator,(const LinearList<T>& lin);

	operator Enum&();
};

/**
 * This is a temporary object that represent a linear list of values to use in
 * CCfg::Enum initialization.
 */
template<typename T>
struct Enum::LinearList
{
	T first, last, step;

	LinearList(T first, T last, T step);
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Enum inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Builds a linear list of values.
 * Its better to call this method than calling the constructor because template
 * methods have automatic template type deduction.
 */
template<typename T>
inline Enum::LinearList<T> Enum::linear(T first, T last, T step)
{
	return Enum::LinearList<T>(first, last, step);
}

/**
 * Returns the sets of values.
 */
inline const Enum::Set& Enum::getValues() const
{
	return values;
}

/**
 * Inserts a new value in the enumeration.
 */
template<typename T>
inline void Enum::insert(T val)
{
	Value v(val);
	checkTypeHomogeneity(v);
	values.insert(v);
}

/**
 * Inserts a new value in the enumeration.
 * This is the specialization for CCfg::Value.
 * @throws Enum::checkTypeHomogeneity
 */
template<>
inline void Enum::insert(Value val)
{
	checkTypeHomogeneity(val);
	values.insert(val);
}

/**
 * Inserts a linear list in the enumeration.
 * @throws Enum::checkTypeHomogeneity
 */
template<typename T> void Enum::insert(const LinearList<T>& lin)
{
	Value firstVal(lin.first);
	checkTypeHomogeneity(firstVal);
	for(T val = lin.first; val <= lin.last; val += lin.step)
	{
		Value v(val);
		values.insert(v);
	}
}

/**
 * Inserts a new value in the enumeration.
 * Type homogeneity is not checked.
 */
template<typename T>
inline void Enum::unchecked_insert(T val)
{
	Value v(val);
	values.insert(v);
}

/**
 * Inserts a new value in the enumeration.
 * This is the specialization for CCfg::Value. Type homogeneity is not checked.
 */
template<>
inline void Enum::unchecked_insert(Value val)
{
	values.insert(val);
}

/**
 * Inserts a linear list in the enumeration.
 * Type homogeneity is not checked.
 */
template<typename T> void Enum::unchecked_insert(const LinearList<T>& lin)
{
	for(T val = lin.first; val <= lin.last; val += lin.step)
	{
		unchecked_insert(val);
	}
}

/**
 * Overwrites an enumeration with a new value and begins an initialization list.
 * The template type can be CCfg::Value or any type accepted by its constructor.
 */
template<typename T>
inline Enum::Inserter<T> Enum::operator=(T val)
{
	values.clear();
	unchecked_insert(val);
	return Inserter<T>(*this);
}

/**
 * Inserts an new value in the enumeration and begins an initialization list.
 * The template type can be CCfg::Value or any type accepted by its constructor.
 * @throws Enum::checkTypeHomogeneity
 */
template<typename T>
inline Enum::Inserter<T> Enum::operator+=(T val)
{
	insert(val);
	return Inserter<T>(*this);
}

/**
 * Overwrites an enumeration with a linear list and begins an initialization list.
 * The template type can be CCfg::Value or any type accepted by its constructor.
 */
template<typename T>
inline Enum::Inserter<T> Enum::operator=(const LinearList<T>& lin)
{
	values.clear();
	unchecked_insert(lin);
	return Inserter<T>(*this);
}

/**
 * Inserts a linear list in the enumeration and begins an initialization list.
 * @throws Enum::checkTypeHomogeneity
 */
template<typename T>
inline Enum::Inserter<T> Enum::operator+=(const LinearList<T>& lin)
{
	insert(lin);
	return Inserter<T>(*this);
}

/**
 * Stream insertion operator.
 */
inline std::ostream& operator<<(std::ostream& stream, const Enum& enumeration)
{
	enumeration.write(stream);
	return stream;
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Enum::Inserter inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
template<typename T>
inline Enum::Inserter<T>::Inserter(Enum& enumeration)
	: enumeration(enumeration)
{
}

/**
 * Inserts a value in the underlying enumeration.
 * @throws Enum::checkTypeHomogeneity
 */
template<typename T>
inline Enum::Inserter<T>& Enum::Inserter<T>::insert(T val)
{
	enumeration.insert(val);
	return *this;
}

/**
 * Inserts a linear list in the underlying enumeration.
 * @throws Enum::checkTypeHomogeneity
 */
template<typename T>
inline Enum::Inserter<T>& Enum::Inserter<T>::insert(const LinearList<T>& lin)
{
	enumeration.insert(lin);
	return *this;
}

/**
 * Inserts the value after each comma in an initialization list.
 */
template<typename T>
inline Enum::Inserter<T>& Enum::Inserter<T>::operator,(T val)
{
	enumeration.unchecked_insert(val);
	return *this;
}

/**
 * Inserts the linear list after each comma in an initialization list.
 */
template<typename T>
inline Enum::Inserter<T>& Enum::Inserter<T>::operator,(const LinearList<T>& lin)
{
	enumeration.unchecked_insert(lin);
	return *this;
}

/**
 * Conversion operator.
 * This operator is needed to keep the assignment return object compatible with
 * standard assignment rules.
 */
template<typename T>
inline Enum::Inserter<T>::operator Enum&()
{
	return enumeration;
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Enumm::LinearList class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
template<typename T>
inline Enum::LinearList<T>::LinearList(T first, T last, T step)
	: first(first), last(last), step(step)
{
}

}

#endif // CCfg_Enum_h_INCLUDED
