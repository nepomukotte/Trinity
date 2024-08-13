/**
 * @file WordIterator.hpp
 * @created Nov 14, 2012
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


#ifndef utl_WordIterator_hpp_INCLUDED
#define utl_WordIterator_hpp_INCLUDED

#include "config.h"
#include "utl/Exception.h"
#include <iterator>
#include <vector>
#include <cstddef>

namespace utl
{
/**
 * This template class allows to iterate over the words of a contiguous memory block according to endianness: the
 * beginning of the memory block is considered to be its least significant word whereas its most significant word
 * corresponds to the last word.
 * @tparam W The word type from which the memory block is made of.
 */
template < typename W >
class WordIterator : public std::iterator<std::random_access_iterator_tag, W>
{
public:
	typedef W WordType;
	struct MisAlignException;

protected:
	template < typename W2 >
	WordIterator(W2* currentWord = 0) : currentWord(currentWord) {}
	template < typename W2 >
	WordIterator(const W2* currentWord = 0) : currentWord(const_cast<W2*>(currentWord)) {}

public:
	template < typename W2 >
	static WordIterator getBegin(W2* headW2, size_t countW2) { return WordIterator(getBeginPtr(headW2, countW2)); }

	template < typename W2 >
	static WordIterator getEnd(W2* headW2, size_t countW2) { return WordIterator(getEndPtr(headW2, countW2)); }

	template < typename W2 >
	static WordIterator getBegin(std::vector<W2>& v) { return WordIterator(getBegin(v.size()>0 ? &v[0] : 0, v.size())); }

	template < typename W2 >
	static WordIterator getEnd(std::vector<W2>& v) { return WordIterator(getEnd(v.size()>0? &v[0] : 0, v.size())); }

	template < typename W2 >
	static WordIterator getBegin(const std::vector<W2>& v) { return WordIterator(getBegin(v.size()>0 ? &v[0] : 0, v.size())); }

	template < typename W2 >
	static WordIterator getEnd(const std::vector<W2>& v) { return WordIterator(getEnd(v.size()>0? &v[0] : 0, v.size())); }

	      WordType& operator *  ()       { return *currentWord; }
	const WordType& operator *  () const { return *currentWord; }
	      WordType& operator [] (std::ptrdiff_t iWord)       { WordIterator it(*this); it.incr(iWord); return *it; }
	const WordType& operator [] (std::ptrdiff_t iWord) const { WordIterator it(*this); it.incr(iWord); return *it; }
	  WordIterator& operator ++ () { incr(); return *this; }
	  WordIterator& operator -- () { decr(); return *this; }
	  WordIterator  operator ++ (int) { WordIterator it(*this); incr(); return it; }
	  WordIterator  operator -- (int) { WordIterator it(*this); decr(); return it; }
	  WordIterator& operator += (std::ptrdiff_t n) { incr(n); return *this; }
	  WordIterator& operator -= (std::ptrdiff_t n) { decr(n); return *this; }
	  WordIterator  operator +  (std::ptrdiff_t n) const { WordIterator it(*this); it += n; return it; }
	  WordIterator  operator -  (std::ptrdiff_t n) const { WordIterator it(*this); it -= n; return it; }
	          bool  operator == (const WordIterator& it) const { return currentWord == it.currentWord; }
	          bool  operator != (const WordIterator& it) const { return currentWord != it.currentWord; }
	          bool  operator >  (const WordIterator& it) const { return isGreaterThan(it); }
	          bool  operator >= (const WordIterator& it) const { return (*this > it) or (*this == it); }
	          bool  operator <  (const WordIterator& it) const { return not (*this >= it); }
	          bool  operator <= (const WordIterator& it) const { return not (*this >  it); }
	std::ptrdiff_t  operator -  (const WordIterator& it) const { return this->minus(it); }

	template < typename W2 >
	static void assertCompat(size_t wordCount) { if ((wordCount * sizeof(W2)) % sizeof(WordType) != 0) throw MisAlignException(); }

#ifdef WORDS_LITTLE_ENDIAN

	static bool isBigEndian() { return false; }

	template < typename W2 >
	static WordType* getBeginPtr(W2* headW2, size_t countW2) { assertCompat<W2>(countW2); return reinterpret_cast<WordType*>(headW2); }
	template < typename W2 >
	static const WordType* getBeginPtr(const W2* headW2, size_t countW2) { assertCompat<W2>(countW2); return reinterpret_cast<const WordType*>(headW2); }
	template < typename W2 >
	static WordType* getEndPtr(W2* headW2, size_t countW2) { return getBeginPtr(headW2 + countW2, countW2); }
	template < typename W2 >
	static const WordType* getEndPtr(const W2* headW2, size_t countW2) { return getBeginPtr(headW2 + countW2, countW2); }

	void incr() { ++currentWord; }
	void decr() { --currentWord; }
	void incr(std::ptrdiff_t n) { currentWord += n; }
	void decr(std::ptrdiff_t n) { currentWord -= n; }
	bool isGreaterThan(const WordIterator& it) const { return currentWord >  it.currentWord; }
	std::ptrdiff_t minus(const WordIterator& it) const { return currentWord -  it.currentWord; }

#else

	static bool isBigEndian() { return true; }

	template < typename W2 >
	static WordType* getBeginPtr(W2* headW2, size_t countW2) { return getEndPtr(headW2, countW2) + countW2*sizeof(W2)/sizeof(WordType); }
	template < typename W2 >
	static const WordType* getBeginPtr(const W2* headW2, size_t countW2) { return getEndPtr(headW2, countW2) + countW2*sizeof(W2)/sizeof(WordType); }
	template < typename W2 >
	static WordType* getEndPtr(W2* headW2, size_t countW2) { assertCompat<W2>(countW2); return reinterpret_cast<WordType*>(headW2) - 1; }
	template < typename W2 >
	static const WordType* getEndPtr(const W2* headW2, size_t countW2) { assertCompat<W2>(countW2); return reinterpret_cast<const WordType*>(headW2) - 1; }

	void incr() { --currentWord; }
	void decr() { ++currentWord; }
	void incr(std::ptrdiff_t n) { currentWord -= n; }
	void decr(std::ptrdiff_t n) { currentWord += n; }
	bool isGreaterThan(const WordIterator& it) const { return currentWord <  it.currentWord; }
	std::ptrdiff_t minus(const WordIterator& it) const { return it.currentWord - currentWord; }

#endif

private:
	mutable WordType* currentWord;
};

template < typename W >
struct WordIterator<W>::MisAlignException : public Exception
{
	MisAlignException() {}
	virtual ~MisAlignException() throw() {}
	virtual void writeMessage(std::ostream& stream) const { stream << "Iterator wordcount incompatible with word type"; }
};

} // namespace utl

#endif /* utl_WordIterator_hpp_INCLUDED */
