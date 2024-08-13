/**
 * @file BitSet.hpp
 * @created Nov 27, 2012
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

#ifndef utl_BitSet_h_INCLUDED
#define utl_BitSet_h_INCLUDED

#include "utl/RtBitSet.hpp"

namespace utl
{

/**
 * Template class representing a "word" made of an arbitrary number of bits set at compile time.
 * @tparam NBITS Fixed number of bits of this bitset.
 * @tparam A Specifies word type on which size the total buffer size of the bitset must be aligned.
 */
template<size_t NBITS, typename A = unsigned int>
class BitSet : public RtBitSet<A>
{
public:
	typedef RtBitSet<A> Base;
	typedef typename Base::Atom Atom;
	typedef typename Base::unsigned_Atom unsigned_Atom; ///< Unsigned version of atom word.
	typedef typename Base::signed_Atom signed_Atom; ///< Signed version of atom word.
	typedef typename Base::Byte Byte; ///< Byte type (signed or unsigned according to bitset sign)

public:
	BitSet(void* x_buffer_ = 0) : Base(NBITS, x_buffer_) { }
	BitSet(Atom* x_buffer_) : Base(NBITS, x_buffer_) { }


	/// Copy constructor: you need to code it separately from other template constructors
	BitSet(const BitSet& aBitSet) : Base(NBITS) { setValue(aBitSet); }

	/// Assignment op: you need to code it separately from other template assignments
	BitSet& operator = (const BitSet& bitSet) { setValue(bitSet); return *this; }

	template < typename W2 >
	BitSet(const W2& initWord) : Base(NBITS) { setValue(initWord); }

	template < typename W2 >
	BitSet(const std::vector<W2>& initVect) : Base(NBITS) { setValue(initVect); }

	template < typename A2 >
	BitSet(const RtBitSet<A2>& initBitSet) : Base(NBITS) { setValue(initBitSet); }

	template < size_t N2, typename A2 >
	BitSet(const BitSet<N2,A2>& initBitSet) : Base(NBITS) { setValue(initBitSet); }

	BitSet& resetBuffer(void* _x_buffer = 0) { return static_cast<BitSet&>(Base::resetBuffer(NBITS, _x_buffer)); }

	static size_t bitCount () { return NBITS; }
	static size_t byteCount() { return atomCount() * Base::atomByteCount(); }
	static size_t atomCount() { return NBITS%Base::atomBitCount() == 0 ? NBITS/Base::atomBitCount() : 1+NBITS/Base::atomBitCount(); }
	static size_t minByteCount() { return bitCount()%CHAR_BIT == 0 ? bitCount()/CHAR_BIT : 1 + bitCount()/CHAR_BIT; }

	///@name Value getters and setters
	///@{
	template < typename W2 >
	void setValue(const W2& w) { Base::setValue(w); }
	template < typename W2 >
	void setValue(const std::vector<W2>& sourceVect) { Base::setValue(sourceVect); }
	void setValue(const unsigned_Atom& w) { Base::setValue(w); }
	void setValue(const signed_Atom& w) { Base::setValue(w); }
	template < typename A2 >
	void setValue(const RtBitSet<A2>& bitSet) { Base::setValue(bitSet); }
	template < size_t N2, typename A2 >
	void setValue(const BitSet<N2,A2>& bitSet) { Base::setValue(static_cast< const RtBitSet<A2>& >(bitSet)); }
	template < typename AtomType >
	void setValue(const AtomType* sourceArray, size_t arraySize) { Base::setValue(sourceArray, arraySize); }

	template < typename W2 >
	void getValue(W2& w) const { Base::getValue(w); }
	template < typename W2 >
	void getValue(std::vector<W2>& targetVect) const { Base::getValue(targetVect); }
	void getValue(unsigned_Atom& w) const { Base::getValue(w); }
	void getValue(signed_Atom& w) const { Base::getValue(w); }
	template < typename A2 >
	void getValue(RtBitSet<A2>& bitSet) const { Base::getValue(bitSet); }
	template < size_t N2, typename A2 >
	void getValue(BitSet<N2,A2>& bitSet) const { Base::getValue(static_cast< RtBitSet<A2>& >(bitSet)); }
	template < typename AtomType >
	void getValue(AtomType* targetArray, size_t arraySize) const { Base::getValue(targetArray, arraySize); }
	///@}

	///@name Bitwise operations
	///@{
	template < typename W2 >
	void bitAnd(const W2& w) { Base::bitAnd(w); }
	void bitAnd(const unsigned_Atom& w) { Base::bitAnd(w); }
	void bitAnd(const signed_Atom& w) { Base::bitAnd(w); }
	template < size_t N2, typename A2 >
	void bitAnd(const RtBitSet<A2>& bitSet) { Base::bitAnd(bitSet); }
	template < size_t N2, typename A2 >
	void bitAnd(const BitSet<N2,A2>& bitSet) { Base::bitAnd(static_cast< RtBitSet<A2>& >(bitSet)); }

	template < typename W2 >
	void bitOr(const W2& w) { Base::bitOr(w); }
	void bitOr(const unsigned_Atom& w) { Base::bitOr(w); }
	void bitOr(const signed_Atom& w) { Base::bitOr(w); }
	template < size_t N2, typename A2 >
	void bitOr(const RtBitSet<A2>& bitSet) { Base::bitOr(bitSet); }
	template < size_t N2, typename A2 >
	void bitOr(const BitSet<N2,A2>& bitSet) { Base::bitOr(static_cast< RtBitSet<A2>& >(bitSet)); }

	template < typename W2 >
	void bitXor(const W2& w) { Base::bitXor(w); }
	void bitXor(const unsigned_Atom& w) { Base::bitXor(w); }
	void bitXor(const signed_Atom& w) { Base::bitXor(w); }
	template < size_t N2, typename A2 >
	void bitXor(const RtBitSet<A2>& bitSet) { Base::bitXor(bitSet); }
	template < size_t N2, typename A2 >
	void bitXor(const BitSet<N2,A2>& bitSet) { Base::bitXor(static_cast< RtBitSet<A2>& >(bitSet)); }

	void bitNot();
	void shiftLeft  (size_t nBits);
	void shiftRight (size_t nBits);

	template < typename W2 >
	BitSet& operator = (const W2& w) { return static_cast<BitSet&>(Base::operator = (w)); }
	template < typename A2 >
	BitSet& operator = (const RtBitSet<A2>& bitSet) { return static_cast<BitSet&>(Base::operator = (bitSet)); }
	template < size_t N2, typename A2 >
	BitSet& operator = (const BitSet<N2,A2>& bitSet) { return static_cast<BitSet&>(Base::operator = (static_cast<const RtBitSet<A2>&>(bitSet))); }

	BitSet& operator <<= (size_t nBits) { return static_cast<BitSet&>(Base::operator <<= (nBits)); }
	BitSet& operator >>= (size_t nBits) { return static_cast<BitSet&>(Base::operator >>= (nBits)); }
	BitSet  operator   ~ ()  const { BitSet bitSetTmp(*this); bitSetTmp.bitNot(); return bitSetTmp; }
	BitSet  operator  << (size_t nBits) const { BitSet bitSetTmp(*this); bitSetTmp <<= nBits; return bitSetTmp; }
	BitSet  operator  >> (size_t nBits) const { BitSet bitSetTmp(*this); bitSetTmp >>= nBits; return bitSetTmp; }

	template < typename W2 >
	BitSet& operator  &= (const W2& w) { return static_cast<BitSet&>(Base::operator &= (w)); }
	template < typename A2 >
	BitSet& operator  &= (const RtBitSet<A2>& bitSet) { return static_cast<BitSet&>(Base::operator &= (bitSet)); }
	template < size_t N2, typename A2 >
	BitSet& operator  &= (const BitSet<N2,A2>& bitSet) { return static_cast<BitSet&>(Base::operator &= (static_cast<const RtBitSet<A2>&>(bitSet))); }

	template < typename W2 >
	BitSet& operator  |= (const W2& w) { return static_cast<BitSet&>(Base::operator |= (w)); }
	template < typename A2 >
	BitSet& operator  |= (const RtBitSet<A2>& bitSet) { return static_cast<BitSet&>(Base::operator |= (bitSet)); }
	template < size_t N2, typename A2 >
	BitSet& operator  |= (const BitSet<N2,A2>& bitSet) { return static_cast<BitSet&>(Base::operator |= (static_cast<const RtBitSet<A2>&>(bitSet))); }

	template < typename W2 >
	BitSet& operator  ^= (const W2& w) { return static_cast<BitSet&>(Base::operator ^= (w)); }
	template < typename A2 >
	BitSet& operator  ^= (const RtBitSet<A2>& bitSet) { return static_cast<BitSet&>(Base::operator ^= (bitSet)); }
	template < size_t N2, typename A2 >
	BitSet& operator  ^= (const BitSet<N2,A2>& bitSet) { return static_cast<BitSet&>(Base::operator ^= (static_cast<const RtBitSet<A2>&>(bitSet))); }

	template < typename W2 >
	BitSet  operator   & (const W2& w) const { BitSet bitSetTmp(*this); bitSetTmp &= w; return bitSetTmp; }
	template < typename W2 >
	BitSet  operator   | (const W2& w) const { BitSet bitSetTmp(*this); bitSetTmp |= w; return bitSetTmp; }
	template < typename W2 >
	BitSet  operator   ^ (const W2& w) const { BitSet bitSetTmp(*this); bitSetTmp ^= w; return bitSetTmp; }
	template < typename W2 >
	operator W2 () const { return toScalar<W2>(); }

	template <typename W2>
	W2 toScalar() const { W2 w; getValue(w); return w; }
};

} // namespace utl

template < size_t NBITS, typename A >
std::ostream& operator << (std::ostream& os, const ::utl::BitSet<NBITS,A>& bitSet)
{
	return os << static_cast< const ::utl::RtBitSet<A>& >(bitSet);
}

#endif /* utl_BitSet_h_INCLUDED */
