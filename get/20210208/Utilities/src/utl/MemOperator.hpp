/**
 * @file MemOperator.hpp
 * @created Oct 8, 2013
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


#ifndef utl_MemOperator_hpp_INCLUDED
#define utl_MemOperator_hpp_INCLUDED

#include "utl/WordIterator.hpp"
#include <boost/integer.hpp>

namespace utl
{

/**
 * General helper function testing positiveness of a scalar word.
 * @tparam W Type of the word to be tested.
 * @param w The word to be tested.
 * @return true if w is positive or null.
 */
template < typename W >
inline bool isPositive(const W& w) { return w >= 0; }

template <>
inline bool isPositive(const unsigned char&) { return true; }

template <>
inline bool isPositive(const unsigned short&) { return true; }

template <>
inline bool isPositive(const unsigned int&) { return true; }

template <>
inline bool isPositive(const unsigned long&) { return true; }

/**
 * General helper function testing strict negativeness of a scalar word.
 * @tparam W Type of the word to be tested.
 * @param w The word to be tested.
 * @return true if w is strictly negative.
 */
template < typename W >
inline bool isNegative(const W& w) { return not isPositive(w); }


///**
// * This structure contains nested structures corresponding to the different elementary operations that can be applied to
// * a single word (this single word type being represented by the Atom type).
// * @tparam A The single word scalar type (the Atom type)
// */
//template < typename A >
//struct AtomOp
//{
//	typedef A AtomType;
//
//	struct copy
//	{
//		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom  = sourceAtom; return targetAtom; }
//	};
//
//	struct bitAnd
//	{
//		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom &= sourceAtom; return targetAtom; }
//	};
//
//	struct bitOr
//	{
//		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom |= sourceAtom; return targetAtom; }
//	};
//
//	struct bitXor
//	{
//		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom ^= sourceAtom; return targetAtom; }
//	};
//
//	struct bitNot
//	{
//		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom = ~sourceAtom; return targetAtom; }
//	};
//};

/**
 * This structure implements a general function applying each of the AtomOp single word operations to a series of
 * consecutive words in memory (burst operation).
 */
template < class A >
struct MemOperator
{
	typedef A AtomType;
	struct copy
	{
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom  = sourceAtom; return targetAtom; }
	};

	struct bitAnd
	{
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom &= sourceAtom; return targetAtom; }
	};

	struct bitOr
	{
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom |= sourceAtom; return targetAtom; }
	};

	struct bitXor
	{
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom ^= sourceAtom; return targetAtom; }
	};

	struct bitNot
	{
		AtomType& operator() (AtomType& targetAtom, AtomType sourceAtom) { targetAtom = ~sourceAtom; return targetAtom; }
	};

	template < class Op >
	static void apply(AtomType* pTarget, size_t maxAtoms, const AtomType* pSource, size_t nAtoms);
};

/**
 * This static method executes in a loop over the source and target word arrays the same elementary operation (e.g. copy
 * operation). In the case where target array is of a greater size than source array, excess target words are "wiped"
 * according to the signedness of the array word type: if AtomType is signed, the sign bit is extended. The operation is
 * a succession of the same operation over an atom reflecting the endianness of the compilation platform, i.e the
 * operation burst is done in the same byte order as it would be done on a single word; you might think of the word
 * array as one big word.
 * @tparam A The type of the elementary element of the arrays.
 * @tparam AtomOp The elementary operation applied (as functor).
 * @param pTarget Pointer to the head of the target array.
 * @param maxAtoms Size of the target array.
 * @param pSource Pointer to the head of the source array.
 * @param nAtoms Size of the source array.
 */
template < class A >
template < class Op >
void MemOperator<A>::apply(AtomType* pTarget, size_t maxAtoms, const AtomType* pSource, size_t nAtoms)
{
	AtomType wipeValue = isPositive <AtomType>(*(WordIterator<AtomType>::getEnd(pSource, nAtoms) - 1)) ? 0 : -1;
	std::ptrdiff_t delta = nAtoms - maxAtoms;

	if (delta > 0)
	{
		nAtoms = maxAtoms;
#ifdef WORDS_BIG_ENDIAN
		pSource += delta;
#endif
	}
	else
	{
#ifdef WORDS_BIG_ENDIAN
		AtomType* pWipe = pTarget;
#else
		AtomType* pWipe = pTarget + nAtoms;
#endif
		while (delta++ < 0)
		{
			Op()(*pWipe, wipeValue);
			++pWipe;
		}
	}
	while (nAtoms-- > 0)
	{
		Op()(*pTarget, *pSource);
		++pTarget;
		++pSource;
	}
}

} // namespace utl



#endif /* utl_MemOperator_hpp_INCLUDED */
