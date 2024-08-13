/**
 * $Id: Channel.h 1029 2013-11-26 08:31:22Z psizun $
 * @file Channel.h
 * @date 6 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: patrick.sizun@cea.fr
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of
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

#ifndef get_Channel_h_INCLUDED
#define get_Channel_h_INCLUDED

#include <boost/cstdint.hpp>
#include <map>
#include <vector>

namespace get {
/**
 * A class to store all values of a channel before plotting.
 */
//_________________________________________________________________________________________________
class Channel
{
public:
	Channel(uint8_t chanIdx=0, uint8_t agetIdx=0);
	virtual ~Channel() {};
	void clear();
	void reserve(const size_t cellCount);
	void addValue(const uint16_t & buckIdx, const uint16_t & sampleValue);
	const std::vector< uint16_t > & buckIndexes() const { return buckIndexes_; }
	const std::vector< uint16_t > & sampleValues() const { return sampleValues_; }
	uint8_t agetIdx() const;
	void setAgetIdx(uint8_t agetIdx);
	uint8_t chanIdx() const;
	void setChanIdx(uint8_t chanIdx);
	bool isFpnChannel() const;
	const uint16_t & sampleAt(const uint16_t & buckIndex) const;
	uint16_t maxSampleValue() const;
	uint16_t minSampleValue() const;
	double meanSampleValue() const;
	size_t sampleCount() const;
	const uint16_t* data() const;
	size_t fpnNeighbourChannel() const;
private:
	struct distance_between_channels : std::binary_function< size_t, size_t , size_t>
	{
		size_t operator()(const size_t & l, const size_t & r) const
		{
			return std::max(l,r) - std::min(l,r);
		}
	};
private:
	static const size_t fpnCount = 4;
	static const size_t fpnChannels[fpnCount];
	std::vector< uint16_t > buckIndexes_; ///< Time bucket indices.
	std::vector< uint16_t > sampleValues_; ///< Digitized values of the signal.
	uint8_t chanIdx_; ///< Channel index (0 to 67)
	uint8_t agetIdx_; /// AGET index within its AsAd container (0 to 3).
};
//_________________________________________________________________________________________________
} /* namespace get */
#endif /* get_Channel_h_INCLUDED */
