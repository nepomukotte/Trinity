/* =============================================================================
 * BinIO.hpp
 * -----------------------------------------------------------------------------
 * class utl::BinIO
 * Created on: 17 March 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the Utilities package of the Mordicus Real-Time
 * Software Development Framework. It implements, mostly in the form of
 * templates, utility functions needed to address endianness issues in output
 * operations.
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
 * commercial purposes. For thism you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#ifndef utl_BinIO_hpp_INCLUDED
#define utl_BinIO_hpp_INCLUDED

#include "config.h"
#include <boost/cstdint.hpp>
#include <boost/version.hpp>
#include <boost/math/special_functions/sign.hpp> // boost::math::signbit
#include <iostream>
#include <cstddef>
#include <limits> // numeric_limits::is_signed

// Make sure host endianness has been properly defined (typically in config.h
// file). If not, compilation is aborted.
#if not defined(WORDS_BIG_ENDIAN) and not defined(WORDS_LITTLE_ENDIAN)
#error "Endianness has not been defined! Maybe you forgot to #include config.h?"
#endif

namespace utl
{

enum Endianness { BigEndian, LittleEndian };

class BinIO
{
public:
	/**
	 * Compile-time word byte swap.
	 * @tparam WordType The type of the word to be byte-swapped.
	 * @param w The word to be byte swapped.
	 */
	template<typename WordType>
	static inline WordType& swap(WordType& w)
	{
		WordType x = w;
		char* ptrx = reinterpret_cast<char*>(&x);
		char* ptrw = reinterpret_cast<char*>(&w) + (sizeof(w) - 1);
		for (size_t i=0; i < sizeof(w); ++i)
		{
			*ptrw-- = *ptrx++;
		}
		return w;
	}

	/**
	 * Compile-time byte swap for any data of fixed length. Useful for data with
	 * non-standard byte size (e.g. 24 bits integer).
	 * @tparam NBYTES The data size in bytes
	 * @param ptr Pointer to the data that should be byte swapped
	 */
	template<unsigned NBYTES>
	static inline void* swap(void* data)
	{
		char bridge;
		char* ptrx = reinterpret_cast<char*>(data);
		char* ptrw = reinterpret_cast<char*>(data) + (NBYTES - 1);
		for (unsigned i=0; i < NBYTES / 2; ++i)
		{
			bridge = *ptrw;
			*ptrw-- = *ptrx;
			*ptrx++ = bridge;
		}
		return data;
	}

	/**
	 * Byte swap for any data with runtime-defined byte size.
	 * @param data Pointer to the data that should be byte swapped
	 * @param nBytes The data size in bytes
	 * @return The data pointer.
	 */
	static inline void* swap(void* data, unsigned nBytes)
	{
		char bridge;
		char* ptrx = static_cast<char*>(data);
		char* ptrw = static_cast<char*>(data) + (nBytes - 1);
		for (unsigned i=0; i < nBytes / 2; ++i)
		{
			bridge = *ptrw;
			*ptrw-- = *ptrx;
			*ptrx++ = bridge;
		}
		return data;
	}

	/**
	 * Compile-time word byte swap if host endianness is different from network byte ordering (big endian).
	 * @tparam WordType The type of the word to be possibly byte swapped.
	 * @param w The word to be possibly byte swapped.
	 */
	template<typename WordType>
	static inline WordType& netswap(WordType& w)
	{
#if WORDS_LITTLE_ENDIAN // if host is little endian
		return swap<WordType>(w);
#else
		return w;
#endif
	}

	/**
	 * Run-time byte swap of any data with runtime-defined byte size if host endianness is different from network byte
	 * ordering (big endian).
	 * @param data Pointer to the data that should be byte swapped
	 * @param nBytes The data size in bytes
	 * @return The data pointer.
	 */
	static inline void* netswap(void* data, unsigned nBytes)
	{
#if WORDS_LITTLE_ENDIAN // if host is little endian
		swap(data, nBytes);
#endif
		return data;
	}

	/**
	 * This function decodes N bytes ordered according to big-endianness BigE
	 * from specified input stream into a word of type W.
	 * @tparam BigE Endianness of the decoded bytes.
	 * @tparam W Type of the output word.
	 * @param N Number of bytes to be decoded (can be different from sizeof(W))
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <Endianness E, typename W>
	static void decode(size_t const N, W& w, std::istream& dataIn)
	{
		char* pw = reinterpret_cast<char*>(&w);

		if(E == BigEndian) // if decoded data is big endian
		{
			if (N < sizeof(W))
			{
				w = 0; // reset word to 0 in case N is strictly smaller than sizeof(W)
				pw += sizeof(W) - N;
				dataIn.read(pw, N);
			}
			else
			{
				if (N > sizeof(W))
				{
					dataIn.ignore(N - sizeof(W));
				}
				dataIn.read(pw, sizeof(W));
			}

#if WORDS_LITTLE_ENDIAN // if host is little endian
			swap(w);
#endif
		}
		else // if decoded data is little endian
		{
			if (N < sizeof(W))
			{
				w = 0; // reset word to 0 in case N is strictly smaller than sizeof(W)
				dataIn.read(pw, N);
			}
			else
			{
				dataIn.read(pw, sizeof(W));
				if (N > sizeof(W))
				{
					dataIn.ignore(N - sizeof(W));
				}
			}
#if WORDS_BIG_ENDIAN // if host is big endian
			swap(w);
#endif
		}
	}

	/**
	 * This function decodes N bytes ordered according to big-endianness BigE
	 * from specified input stream into a word of type W.
	 * @tparam BigE Endianness of the decoded bytes.
	 * @tparam W Type of the output word.
	 * @tparam N Number of bytes to be decoded (can be different from sizeof(W))
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <Endianness E, typename W, unsigned N>
	static void decode(W& w, std::istream& dataIn)
	{
		decode<E,W>(N, w, dataIn);
	}

	/**
	 * This function decodes sizeof(W) bytes ordered according to endianness E
	 * from specified input stream into a word of type W.
	 * @tparam E Endianness of the decoded bytes.
	 * @tparam W Type of the output word.
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <Endianness E, typename W>
	static void decode(W& w, std::istream& dataIn)
	{
		decode<E, W, sizeof(W)>(w, dataIn);
	}

	/**
	 * This function decodes N bytes ordered according to endianness e from
	 * specified input stream into a word of type W.
	 * @tparam W Type of the output word.
	 * @tparam N Number of bytes to be decoded (can be different from sizeof(W))
	 * @param e Endianness of the decoded bytes.
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <typename W, unsigned N>
	static void decode(Endianness e, W& w, std::istream& dataIn)
	{
		if (e == BigEndian)
		{
			decode<BigEndian, W, N>(w, dataIn);
		}
		else
		{
			decode<LittleEndian, W, N>(w, dataIn);
		}
	}

	/**
	 * This function decodes n bytes ordered according to endianness e from
	 * specified input stream into a word of type W.
	 * @tparam W Type of the output word.
	 * @param n Number of bytes to be decoded (can be different from sizeof(W))
	 * @param e Endianness of the decoded bytes.
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <typename W>
	static void decode(Endianness e, size_t const n, W& w, std::istream& dataIn)
	{
		if (e == BigEndian)
		{
			decode<BigEndian, W>(n, w, dataIn);
		}
		else
		{
			decode<LittleEndian, W>(n, w, dataIn);
		}
	}

	/**
	 * This function decodes N bytes ordered according to endianness of host
	 * from specified input stream into a word of type W.
	 * @tparam W Type of the output word.
	 * @tparam N Number of bytes to be decoded (can be different from sizeof(W))
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <typename W, unsigned N>
	static void decode(W& w, std::istream& dataIn)
	{
#if WORDS_BIG_ENDIAN // if host is big endian

		decode<BigEndian, W, N>(w, dataIn);

#elif WORDS_LITTLE_ENDIAN // if host is little endian

		decode<LittleEndian, W, N>(w, dataIn);

#endif
	}

	/**
	 * This function decodes sizeof(W) bytes ordered according to host
	 * endianness from specified input stream into a word of type W.
	 * @tparam W Type of the output word.
	 * @param w Output word into which bytes are decoded.
	 * @param dataIn Binary stream from which the N bytes are read out.
	 */
	template <typename W>
	static void decode(W& w, std::istream& dataIn)
	{
#if WORDS_BIG_ENDIAN // if host is big endian

		decode<BigEndian, W, sizeof(W)>(w, dataIn);

#elif WORDS_LITTLE_ENDIAN // if host is little endian

		decode<LittleEndian, W, sizeof(W)>(w, dataIn);

#endif
	}

	/**
	 * This function encodes a word of type W into N bytes ordered according to
	 * endianness E to specified output stream.
	 * @tparam E Endianness of the encoded bytes.
	 * @tparam W Type of the input word.
	 * @param N Number of bytes to be encoded (can be different from sizeof(W))
	 * @param w Output word from which bytes are encoded.
	 * @param dataOut Binary stream into which the N bytes are written.
	 */
	template <Endianness E, typename W>
	static void encode(size_t const N, const W& _w, std::ostream& dataOut)
	{
		W w = _w;
		char* pw = reinterpret_cast<char*>(&w);

		if(E == BigEndian) // if encoded data is big endian
		{
#if WORDS_LITTLE_ENDIAN // if host is little endian
			swap(w);
#endif
			if (N < sizeof(W))
			{
				pw += sizeof(W) - N;
				dataOut.write(pw, N);
			}
			else
			{
				if (N > sizeof(W))
				{
					char filler = '\0';
					if (std::numeric_limits<W>::is_signed
							and boost::math::signbit(_w) != 0)
					{
						filler = ~'\0';
					}
					for (unsigned i=sizeof(W); i < N; ++i)
					{
						dataOut.put(filler); // Pad MSBs with copies of sign bit
					}
				}
				dataOut.write(pw, sizeof(W));
			}
		}
		else // if encoded data is little endian
		{
#if WORDS_BIG_ENDIAN // if host is big endian
			swap(w);
#endif
			if (N < sizeof(W))
			{
				dataOut.write(pw, N);
			}
			else
			{
				dataOut.write(pw, sizeof(W));
				if (N > sizeof(W))
				{
					char filler = '\0';
					if (std::numeric_limits<W>::is_signed
							and boost::math::signbit(_w) != 0)
					{
						filler = ~'\0';
					}
					for (unsigned i=sizeof(W); i < N; ++i)
					{
						dataOut.put(filler); // Pad MSBs with copies of sign bit
					}
				}
			}
		}
	}

	/**
	 * This function encodes a word of type W into N bytes ordered according to
	 * endianness E to specified output stream.
	 * @tparam E Endianness of the encoded bytes.
	 * @tparam W Type of the input word.
	 * @tparam N Number of bytes to be encoded (can be different from sizeof(W))
	 * @param w Output word from which bytes are encoded.
	 * @param dataOut Binary stream into which the N bytes are written.
	 */
	template <Endianness E, typename W, unsigned N>
	static void encode(const W& _w, std::ostream& dataOut)
	{
		encode<E,W>(N, _w, dataOut);
	}

	/**
	 * This function encodes sizeof(W) bytes ordered according to endianness E
	 * into specified output stream from a word of type W.
	 * @tparam E Endianness of the encoded bytes.
	 * @tparam W Type of the input word.
	 * @param w Input word from which bytes are encoded.
	 * @param dataOut Binary stream into which the N bytes are written.
	 */
	template <Endianness E, typename W>
	static void encode(const W& w, std::ostream& dataOut)
	{
		encode<E, W, sizeof(W)>(w, dataOut);
	}

	/**
	 * This function encodes N bytes ordered according to endianness e into
	 * specified output stream from a word of type W.
	 * @tparam W Type of the input word.
	 * @tparam N Number of bytes to be encoded (can be different from sizeof(W))
	 * @param e Endianness of the encoded bytes.
	 * @param w Input word frm which bytes are encoded.
	 * @param dataOut Binary stream into which the N bytes are written.
	 */
	template <typename W, unsigned N>
	static void encode(Endianness e, const W& w, std::ostream& dataOut)
	{
		if (e == BigEndian)
		{
			encode<BigEndian, W, N>(w, dataOut);
		}
		else
		{
			encode<LittleEndian, W, N>(w, dataOut);
		}
	}

	/**
	 * This function encodes n bytes ordered according to endianness e into
	 * specified output stream from a word of type W.
	 * @tparam W Type of the input word.
	 * @param n Number of bytes to be encoded (can be different from sizeof(W))
	 * @param e Endianness of the encoded bytes.
	 * @param w Input word frm which bytes are encoded.
	 * @param dataOut Binary stream into which the N bytes are written.
	 */
	template <typename W>
	static void encode(Endianness e, size_t const n, const W& w, std::ostream& dataOut)
	{
		if (e == BigEndian)
		{
			encode<BigEndian, W>(n, w, dataOut);
		}
		else
		{
			encode<LittleEndian, W>(n, w, dataOut);
		}
	}

	/**
	 * This function encodes sizeof(W) bytes ordered according to host
	 * endianness into specified output stream from a word of type W.
	 * @tparam W Type of the input word.
	 * @param w Input word from which bytes are encoded.
	 * @param dataOut Binary stream into which the N bytes are written.
	 */
	template <typename W>
	static void encode(const W& w, std::ostream& dataOut)
	{
#if WORDS_BIG_ENDIAN // if host is big endian

		encode<BigEndian, W>(w, dataOut);

#elif WORDS_LITTLE_ENDIAN // if host is little endian

		encode<LittleEndian, W>(w, dataOut);

#endif
	}

	static const Endianness HOST_ENDIANNESS =
#if WORDS_BIG_ENDIAN
		BigEndian;
#elif WORDS_LITTLE_ENDIAN
		LittleEndian;
#endif

private:
	BinIO() {}
	~BinIO() {}
};

} // namespace utl

#endif /* utl_BinIO_hpp_INCLUDED */
