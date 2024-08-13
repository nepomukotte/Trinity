// =============================================================================
// PointerValueStreamer.hpp
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

#ifndef Utils_Io_PointerValueStreamer_hpp_INCLUDED
#define Utils_Io_PointerValueStreamer_hpp_INCLUDED

#include <ios>
#include <ostream>

namespace Utils
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 *
 * @author fchateau
 */
template<typename T>
struct PointerValueStreamer
{
	const T* ptr;
	explicit PointerValueStreamer(const T* ptr);
	PointerValueStreamer(const PointerValueStreamer& rvalue);
	~PointerValueStreamer();

private:
	PointerValueStreamer();
	PointerValueStreamer& operator=(const PointerValueStreamer& rvalue);
};

template<typename T>
PointerValueStreamer<T> writeValue(const T* ptr);

template<typename T>
std::ostream& operator<<(std::ostream& stream, const PointerValueStreamer<T>& ptrStreamer);

////////////////////////////////////////////////////////////////////////////////
// Utils::Io::PointerValueStreamer inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
template<typename T>
inline PointerValueStreamer<T>::PointerValueStreamer(const T* ptr)
	: ptr(ptr)
{
}

/**
 * Copy constructor.
 */
template<typename T>
inline PointerValueStreamer<T>::PointerValueStreamer(const PointerValueStreamer& rvalue)
	: ptr(rvalue.ptr)
{
}

/**
 * Destructor.
 */
template<typename T>
inline PointerValueStreamer<T>::~PointerValueStreamer()
{
}

////////////////////////////////////////////////////////////////////////////////
// Namespace inline methods body
////////////////////////////////////////////////////////////////////////////////

template<typename T>
inline PointerValueStreamer<T> writeValue(const T* ptr)
{
	return PointerValueStreamer<T>(ptr);
}

/**
 * Stream insertion operator.
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& stream, const PointerValueStreamer<T>& ptrStreamer)
{
	if(ptrStreamer.ptr == 0)
	{
		stream << "null";
	}
	else
	{
		stream << *ptrStreamer.ptr;
	}
	return stream;
}


}
}

#endif // Utils_Io_PointerValueStreamer_hpp_INCLUDED
