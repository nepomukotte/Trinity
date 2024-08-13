/* =============================================================================
 * Common.h
 * -----------------------------------------------------------------------------
 * class mfm::Common
 * -----------------------------------------------------------------------------
 * Created on: 16 mars 2010
 * Author: Shebli Anvar - Irfu/Sedi/Lilas - CEA Saclay - shebli.anvar@cea.fr
 * -----------------------------------------------------------------------------
 * This software is the property of CEA (Commissariat a l'Energie Atomique).
 * You may use it freely for internal and non commercial purposes. Any other
 * use is subject to a written authorization of a legal representative of CEA.
 * =============================================================================
 */

#ifndef mfm_Common_h_INCLUDED
#define mfm_Common_h_INCLUDED

#include <boost/cstdint.hpp>
#include <boost/dynamic_bitset.hpp>
#include <climits>
#include <vector>

/**
 * Namespace of the MultiFrame Metaformat (MFM) API.
 */
namespace mfm
{

/// Data type used to represent a single Byte within a frame.
typedef char Byte;
/// Container used to store sequences of bytes.
typedef ::std::vector< Byte > ByteArray;
/// Container used to provide access to individual bits.
typedef ::boost::dynamic_bitset< unsigned char > DynamicBitset; // 'char' not supported by dynamic_bitset ?

/**
 * As of specification MFM2.1, there are 3 MFM metatypes corresponding to the 3
 * types of data frames which form the base structures of application frames:
 * Blob frames, Basic frames and Layered frames.
 */
enum FrameKind
{
	BLOB,   ///< Value of mfm::FrameKind associated to blob frames, with a data section consisting of an opaque block of data whose internal structure is not described in the header.
	BASIC,  ///< Value of mfm::FrameKind associated to basic frames, whose data section contains fixed size items.
	LAYERED ///< Value of mfm::FrameKind associated to layered frames, whose data section consists of variable size items which are Frames themselves.
};

/**
 * This is a utility template function raising 2 to an integer power. Used to
 * avoid floating point pow() and log() functions. CAUTION: overflow cases are
 * NOT checked so result will be false in case of overflow!
 * @param exponent The power exponent.
 * @return 2^exponent
 */
template<typename ResultType, typename SeedType>
inline ResultType pow2(SeedType exponent)
{
	register ResultType result = 1;
	return result << exponent;
}

/**
 * This is a utility template function giving the base 2 logarithm of a power of
 * 2. Used to avoid floating point pow() and log() functions.
 * CAUTION: pathological cases are not considered. In particular, if the seed is
 * NOT a power of 2, the result will be the same as that of the greater power
 * of 2 that is smaller than the seed (or zero if seed is zero).
 * @param seed The power of 2 for which we seek the base 2 logarithm.
 * @return log2(seed)
 */
template<typename ResultType, typename SeedType>
inline ResultType logpow2(SeedType seed)
{
	register ResultType result = 0;
	register int maxBits = sizeof(SeedType) * CHAR_BIT;
	for (register int i=0; i < maxBits; ++i)
	{
		switch (seed >> result)
		{
		case 1:
			return result;

		case 0:
			result -= maxBits >> i;
			break;

		default:
			result += maxBits >> i;
			break;
		}
	}
	return result;
}

} // namespace mfm

#endif /* mfm_Common_h_INCLUDED */
