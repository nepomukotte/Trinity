/**
 * @file BitFieldHelperBitSet.hpp
 * @created Nov 23, 2012
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

#ifndef utl_BitFieldHelperBitSet_hpp_INCLUDED
#define utl_BitFieldHelperBitSet_hpp_INCLUDED

#include "config.h"
#include "utl/BitSet.hpp"
#include <algorithm>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

namespace utl
{

/**
 * Template specialization of BitFieldHelper for BitSet objects.
 */
template < size_t N, typename A >
class BitFieldHelper < BitSet<N,A> >
{
public:
	typedef BitSet<N,A> bitset_t;
	typedef typename bitset_t::Atom Atom;
	typedef WordIterator<Atom> iter_Atom;

public:
	/// @name Compile-time bit-field manipulation.
	//@{
	template < unsigned pos >
	static bool getBit(const bitset_t& bitSet) { return bitSet.getBit(pos); }

	template < unsigned pos >
	static void setBit(bitset_t& bitSet, bool askSet) { return bitSet.setBit(pos, askSet); }

	template < unsigned pos, unsigned wid, typename value_t >
	static void getField(const bitset_t& bitSet, value_t& value) { bitSet.getField(pos, wid, value); }

	template < unsigned pos, unsigned wid, typename value_t >
	static value_t getField(const bitset_t& bitSet) { return bitSet.getField(pos, wid); }

	template < unsigned pos, unsigned wid, typename value_t >
	static void setField(bitset_t& bitSet, const value_t& value) { bitSet.setField(pos, wid, value); }

	/// @name Run-time bit-field manipulation.
	//@{
	static bool getBit(const bitset_t& bitSet, unsigned pos) { return bitSet.getBit(pos); }

	static void setBit(bitset_t& bitSet, unsigned pos, bool askSet) { bitSet.setBit(pos, askSet); }

	template < typename value_t >
	static void getField(const bitset_t& bitSet, unsigned pos, unsigned wid, value_t& value) { bitSet.getField(pos, wid, value); }

	template < typename value_t >
	static value_t getField(const bitset_t& bitSet, unsigned pos, unsigned wid) { return bitSet.getField(pos, wid); }

	template < typename value_t >
	static void setField(bitset_t& bitSet, unsigned pos, unsigned wid, const value_t& value) { bitSet.setField(pos, wid, value); }
	//@}

private:
	BitFieldHelper() {}
	~BitFieldHelper() {}
};

} // namespace utl

#endif /* utl_BitFieldHelperBitSet_hpp_INCLUDED */
