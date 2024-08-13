/* =============================================================================
 * Rounding.h
 * -----------------------------------------------------------------------------
 * Created on: 19 August 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of general utilities package.
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

#ifndef utl_Rounding_INCLUDED
#define utl_Rounding_INCLUDED

#include <boost/type_traits/make_signed.hpp>
#include <boost/type_traits/make_unsigned.hpp>

namespace utl
{
// enums
enum RoundMode {RoundDown, RoundNear, RoundUp};

template<RoundMode round, typename T>
struct do_divide
{
	static T apply(T numerator, T denominator);
};

template<typename T>
struct do_divide<RoundDown, T>
{
	static T apply(T numerator, T denominator)
	{
		return numerator / denominator;
	}
};

template<typename T>
struct do_divide<RoundNear, T>
{
	typedef typename boost::make_signed<T>::type Signed;
	typedef typename boost::make_unsigned<T>::type Unsigned;

	static Unsigned apply(Unsigned numerator, Unsigned denominator)
	{
		return (numerator + denominator / 2) / denominator;
	}
	static Signed apply(Signed numerator, Signed denominator)
	{
		Signed result = numerator;
		if((numerator < 0) xor (denominator < 0))
		{
			result -= denominator / 2;
		}
		else
		{
			result += denominator / 2;
		}
		return result / denominator;
	}
};

template<typename T>
struct do_divide<RoundUp, T>
{
	typedef typename boost::make_signed<T>::type Signed;
	typedef typename boost::make_unsigned<T>::type Unsigned;

	static Unsigned apply(Unsigned numerator, Unsigned denominator)
	{
		return (numerator + denominator - 1) / denominator;
	}
	static Signed apply(Signed numerator, Signed denominator)
	{
		Signed result = numerator;
		if(numerator < 0)
		{
			if(denominator < 0)
				result += denominator + 1;
			else
				result += -denominator + 1;
		}
		else
		{
			if(denominator < 0)
				result += -denominator - 1;
			else
				result += denominator - 1;
		}
		return result / denominator;
	}
};

/**
 * Performs an integer division and rounds the result according to the template parameter.
 */
template<RoundMode round, typename T>
T divide(T numerator, T denominator)
{
	return do_divide<round, T>::apply(numerator, denominator);
}

/**
 * Aligns a value on the the upper/lower/nearest multiple of the boundary, according to the RoundMode.
 */
template<RoundMode round, typename T>
T align(T value, T boundary)
{
	return divide<round>(value, boundary) * boundary;
}

}

#endif // utl_Rounding_INCLUDED
