/**
 * @file InplaceParser.hpp
 * @date 5 juin 2013
 * @author fchateau
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau
 *
 * This file is part of the Utilities software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#ifndef utl_io_InplaceParser_hpp_INCLUDED
#define utl_io_InplaceParser_hpp_INCLUDED

#include "Buffer.hpp"

namespace utl
{
namespace io
{
/**
 * This class encapsulate a non-owning pointer to a single data object.
 * The object pointer may be set or not.
 */
template<class T>
class ValueRef
{
public:
	typedef T DataField;

	ValueRef() : ref_(0) {}

	bool isSet() const { return ref_ != 0; }
	const DataField& get() const { return *ref_; }

	operator const DataField&() const { return *ref_; }
	const DataField* operator->() const { return ref_; }

	const DataField*& get() { return ref_; }
	void clear() { ref_ = 0; }
private:
	const DataField* ref_;
};

/**
 * This class encapsulate a non-owning pointer to an array of contiguous objects.
 * The number of objects is arbitrary and the pointer can be unset or the array empty.
 */
template<class T>
class ArrayRef
{
public:
	typedef T DataField;

	ArrayRef() : count_(0), ref_(0) {}

	bool isSet() const { return ref_ != 0; }
	unsigned count() const { return count_; }
	const DataField* get() const { return ref_; }
	const DataField& operator[](unsigned i) const { return ref_[i]; }

	unsigned& count() { return count_; }
	const DataField*& get() { return ref_; }
	void clear() { count_ = 0; ref_ = 0; }

private:
	unsigned count_;
	const DataField* ref_;
};

/**
 * This class enables in-place parsing of data, which means decoding some data
 * in a buffer and binding POD data on it, so that no copy of the data is necessary.
 * Binding is achieved through ValueRef and ArrayRef classes.
 *
 * For more complex data layout, where variable-length or optional data are
 * contained into arrays, direct binding with a predefined number of ValueRef and
 * ArrayRef might not be possible: defining vectors of them is necessary.
 *
 * In order to enable parsing of user-defined binding fields (structures made of
 * several Value/ArrayRef, the parse_1(BindingField&) method is here to be specialized.
 * The specializations must reside inside a utl::io namespace block.
 *
 * Parsing a DataField can be either unchecked or checked (the data is examined to
 * determine if it matches the DataField). for a DataField to be checked a
 * bool isDataField(const InplaceParser&, const DataField*) must exist in the utl::io
 * namespace. A user-defined overload can be defined for any type (or group of
 * type, if using a template method and SFINAE). It has to be declared prior to
 * including InplaceParser.ipp.
 */
class InplaceParser
{
public:
	template<class T>
	struct has_checker;

	InplaceParser();
	InplaceParser(const ConstBuffer& buf);
	InplaceParser(const InplaceParser& p);
	InplaceParser& operator=(const InplaceParser& val);

	template<class BindingField>
	bool parse_1(BindingField& f);

	template<class DataField>
	bool parse_1(ValueRef<DataField>& dest);

	template<class DataField>
	bool parse_N(ArrayRef<DataField>& dest, unsigned count);

	template<class DataField>
	bool parse_0_1(ValueRef<DataField>& dest);

	template<class DataField>
	bool parse_0_N(ArrayRef<DataField>& dest);

	template<class DataField>
	bool parse_1_N(ArrayRef<DataField>& dest);

	template<class BindingField>
	bool parse_0_N(std::vector<BindingField>& dest);

	template<class BindingField>
	bool parse_1_N(std::vector<BindingField>& dest);

	template<class T>
	bool bind(const T*& data, size_t n = 1);

	template<class T>
	bool tryBind(const T*& data, size_t n = 1) const;

	template<class T>
	bool ignore(size_t n = 1);

	bool commit(InplaceParser& dest) const;

	const void* cur() const { return bufCurPtr; }
	const void* end() const { return bufEndPtr; }

	bool isEndOfBuffer() const;

private:
	const char* bufCurPtr;
	const char* bufEndPtr;
};


} /* namespace io */
} /* namespace utl */

#endif /* utl_io_InplaceParser_hpp_INCLUDED */
