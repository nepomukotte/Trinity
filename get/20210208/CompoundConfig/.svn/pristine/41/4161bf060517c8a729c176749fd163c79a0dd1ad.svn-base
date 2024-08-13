// =============================================================================
// float.cpp
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

#include "float.h"
#include <limits>
#include <cmath>

namespace Utils
{
namespace Math
{
/**
 * Compares two real numbers for equality like x == y, using an epsilon.
 * Remark: when comparing 0 to another number, even 1e-300, it always returns
 * false because it assumes numbers of different magnitude are different.
 * @param x left operand
 * @param y right operand
 * @param epsilon a number so that equals(z, z*(1+epsilon), epsilon) == true.
 */
template<typename R> bool equals(R x, R y, R epsilon)
{
	if(x == y)
	{
		return true;
	}

	int ex = 0, ey = 0;
	R mx = frexp(x, &ex);
	R my = frexp(y, &ey);
	if(ex != ey) // different exponent
	{
		return false;
	}

	return fabs(mx - my) < epsilon;
}

/**
 * Returns the smallest number that compares greater than the specified number.
 */
template<typename R> R leastAbove(R x)
{
	int exponent;
	R mantissa = frexp(x, &exponent);
	mantissa += std::numeric_limits<R>::epsilon();
	return ldexp(mantissa, exponent);
}

/**
 * Return the greatest number that compares lesser than than the specified number.
 */
template<typename R> R greatestBelow(R x)
{
	int exponent;
	R mantissa = frexp(x, &exponent);
	mantissa -= std::numeric_limits<R>::epsilon();
	return ldexp(mantissa, exponent);
}

template bool equals(float, float, float);
template bool equals(double, double, double);
template bool equals(long double, long double, long double);

template float leastAbove(float);
template double leastAbove(double);
template long double leastAbove(long double);

template float greatestBelow(float);
template double greatestBelow(double);
template long double greatestBelow(long double);
}
}
