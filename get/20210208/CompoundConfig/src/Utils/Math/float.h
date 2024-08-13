// =============================================================================
// float.h
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

#ifndef Utils_Math_Float_h_INCLUDED
#define Utils_Math_Float_h_INCLUDED

namespace Utils
{
namespace Math
{
////////////////////////////////////////////////////////////////////////////////
// Namespace methods declarations
////////////////////////////////////////////////////////////////////////////////

template<typename R> bool equals(R x, R y, R epsilon);
template<typename R> R leastAbove(R x);
template<typename R> R greatestBelow(R x);

}
}
#endif // Utils_Math_Float_h_INCLUDED
