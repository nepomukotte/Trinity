/* =============================================================================
 * Buffer.hpp
 * -----------------------------------------------------------------------------
 * class utl::io::MutableBuffer
 * class utl::io::ConstBuffer
 * Created on: 12 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This file is part of the Utilities software project.
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
 * =============================================================================
 */

#ifndef utl_io_Buffer_hpp_INCLUDED
#define utl_io_Buffer_hpp_INCLUDED

#include <boost/array.hpp>
#include <boost/cstdint.hpp>
#include <vector>
#include <string>

#ifdef min
#undef min
#endif

namespace utl
{
namespace io
{

struct MutableBuffer
{
	MutableBuffer() : data(), maxBytes() {}
	MutableBuffer(void* data, size_t maxBytes) : data(data), maxBytes(maxBytes) {}

	uint8_t* begin() const { return reinterpret_cast<uint8_t*>(data); }
	uint8_t* end() const { return begin() + maxBytes; }

	void*  data;
	size_t maxBytes;
};

struct ConstBuffer
{
	ConstBuffer() : data(), maxBytes() {}
	ConstBuffer(const void* data, size_t maxBytes) : data(data), maxBytes(maxBytes) {}
	ConstBuffer(const MutableBuffer& b) : data(b.data), maxBytes(b.maxBytes) {}

	const uint8_t* begin() const { return reinterpret_cast<const uint8_t*>(data); }
	const uint8_t* end() const { return begin() + maxBytes; }

	const void* data;
	size_t maxBytes;
};

// =============================================================================
// Buffer over pointer
// -----------------------------------------------------------------------------

template<class T>
ConstBuffer buffer(const T* ptr, size_t nElems) { return ConstBuffer(ptr, nElems * sizeof(T)); }

template<class T>
MutableBuffer buffer(T* ptr, size_t nElems) { return MutableBuffer(ptr, nElems * sizeof(T)); }

// =============================================================================
// Buffer over pointer range
// -----------------------------------------------------------------------------

template<class T>
ConstBuffer buffer(const T* begin, const T* end) { return ConstBuffer(begin, (end - begin) * sizeof(T)); }

template<>
inline ConstBuffer buffer(const void* begin, const void* end) { return buffer(reinterpret_cast<const uint8_t*>(begin), reinterpret_cast<const uint8_t*>(end)); }

template<class T>
MutableBuffer buffer(T* begin, T* end) { return MutableBuffer(begin, (end - begin) * sizeof(T)); }

template<>
inline MutableBuffer buffer(void* begin, void* end) { return buffer(reinterpret_cast<uint8_t*>(begin), reinterpret_cast<uint8_t*>(end)); }

// =============================================================================
// Buffer over C array
// -----------------------------------------------------------------------------

template<class T, size_t N>
ConstBuffer buffer(const T (&array)[N]) { return ConstBuffer(&array[0], N * sizeof(T)); }

template<class T, size_t N>
MutableBuffer buffer(T (&array)[N]) { return MutableBuffer(&array[0], N * sizeof(T)); }

template<class T, size_t N>
ConstBuffer buffer(const T (&array)[N], size_t maxBytes) { return ConstBuffer(&array[0], std::min(maxBytes, N * sizeof(T))); }

template<class T, size_t N>
MutableBuffer buffer(T (&array)[N], size_t maxBytes) { return MutableBuffer(&array[0], std::min(maxBytes, N * sizeof(T))); }

// =============================================================================
// Buffer over boost array
// -----------------------------------------------------------------------------

template<class T, size_t N>
ConstBuffer buffer(const boost::array<T, N>& v) { return ConstBuffer(&v[0], v.size() * sizeof(T)); }

template<class T, size_t N>
MutableBuffer buffer(boost::array<T, N>& v) { return MutableBuffer(&v[0], v.size() * sizeof(T)); }

template<class T, size_t N>
ConstBuffer buffer(const boost::array<T, N>& v, size_t maxBytes) { return ConstBuffer(&v[0], std::min(maxBytes, v.size() * sizeof(T))); }

template<class T, size_t N>
MutableBuffer buffer(boost::array<T, N>& v, size_t maxBytes) { return MutableBuffer(&v[0], std::min(maxBytes, v.size() * sizeof(T))); }

// =============================================================================
// Buffer over std::vector
// -----------------------------------------------------------------------------

template<class T, class A>
ConstBuffer buffer(const std::vector<T, A>& v)  { return ConstBuffer(&v[0], v.size() * sizeof(T)); }

template<class T, class A>
MutableBuffer buffer(std::vector<T, A>& v) { return MutableBuffer(&v[0], v.size() * sizeof(T)); }

template<class T, class A>
ConstBuffer buffer(const std::vector<T, A>& v, size_t maxBytes)  { return ConstBuffer(&v[0], std::min(maxBytes, v.size() * sizeof(T))); }

template<class T, class A>
MutableBuffer buffer(std::vector<T, A>& v, size_t maxBytes) { return MutableBuffer(&v[0], std::min(maxBytes, v.size() * sizeof(T))); }

// =============================================================================
// Buffer over std::string
// -----------------------------------------------------------------------------

template<class E, class T, class A>
ConstBuffer buffer(const std::basic_string<E, T, A>& s) { return ConstBuffer(&s[0], s.size() * sizeof(E)); }

template<class E, class T, class A>
MutableBuffer buffer(std::basic_string<E, T, A>& s) { return MutableBuffer(&s[0], s.size() * sizeof(E)); }

template<class E, class T, class A>
ConstBuffer buffer(const std::basic_string<E, T, A>& s, size_t maxBytes) { return ConstBuffer(&s[0], std::min(maxBytes, s.size() * sizeof(E))); }

template<class E, class T, class A>
MutableBuffer buffer(std::basic_string<E, T, A>& s, size_t maxBytes) { return MutableBuffer(&s[0], std::min(maxBytes, s.size() * sizeof(E))); }

// =============================================================================
// Buffer over POD and primitive types
// -----------------------------------------------------------------------------

template<class T>
ConstBuffer buffer(const T& var) { return ConstBuffer(&var, sizeof(T)); }

template<class T>
MutableBuffer buffer(T& var) { return MutableBuffer(&var, sizeof(T)); }

} /* namespace io */
} /* namespace utl */

#endif /* utl_io_Buffer_hpp_INCLUDED */
