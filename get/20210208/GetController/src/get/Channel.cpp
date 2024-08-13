/**
 * $Id: Channel.cpp 1029 2013-11-26 08:31:22Z psizun $
 * @file Channel.cpp
 * @date 6 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors:
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

#include "get/Channel.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <numeric>
#include <stdexcept>

namespace get {
//_________________________________________________________________________________________________
const size_t Channel::fpnChannels[fpnCount] = {11, 22, 45, 56};
//_________________________________________________________________________________________________
Channel::Channel(uint8_t chanIdx, uint8_t agetIdx) : chanIdx_(chanIdx), agetIdx_(agetIdx)
{
	reserve(512);
}
//_________________________________________________________________________________________________
void Channel::clear()
{
	buckIndexes_.clear();
	sampleValues_.clear();
}
//_________________________________________________________________________________________________
void Channel::reserve(const size_t cellCount)
{
	buckIndexes_.reserve(cellCount);
	sampleValues_.reserve(cellCount);
}
//_________________________________________________________________________________________________
void Channel::addValue(const uint16_t & buckIdx, const uint16_t & sampleValue)
{
	buckIndexes_.push_back(buckIdx);
	sampleValues_.push_back(sampleValue);
}
//_________________________________________________________________________________________________
uint8_t Channel::agetIdx() const
{
	return agetIdx_;
}
//_________________________________________________________________________________________________
void Channel::setAgetIdx(uint8_t agetIdx)
{
	this->agetIdx_ = agetIdx;
}
//_________________________________________________________________________________________________
uint8_t Channel::chanIdx() const
{
	return chanIdx_;
}
//_________________________________________________________________________________________________
void Channel::setChanIdx(uint8_t chanIdx)
{
	this->chanIdx_ = chanIdx;
}
//_________________________________________________________________________________________________
bool Channel::isFpnChannel() const
{
	return (std::find(fpnChannels, fpnChannels+fpnCount, chanIdx_) != fpnChannels+fpnCount);
}
//_________________________________________________________________________________________________
/**
 * Returns sample value at given bucket index.
 * @param buckIndex Bucket index.
 * @return Sample value, if bucket index exists.
 * @throws std::out_of_range
 */
const uint16_t & Channel::sampleAt(const uint16_t & buckIndex) const
{
	typedef std::vector< uint16_t >::const_iterator Iter;
	Iter pos = std::find(buckIndexes_.begin(), buckIndexes_.end(), buckIndex);
	if (pos == buckIndexes_.end())
	{
		throw std::out_of_range("Bucket index does not exist!");
	}
	return sampleValues_.at(std::distance(buckIndexes_.begin(), pos));
}
//_________________________________________________________________________________________________
uint16_t Channel::maxSampleValue() const
{
	if (sampleValues_.empty()) return 0;
	return *std::max_element(sampleValues_.begin(), sampleValues_.end());
}
//_________________________________________________________________________________________________
uint16_t Channel::minSampleValue() const
{
	if (sampleValues_.empty()) return 0;
	return *std::min_element(sampleValues_.begin(), sampleValues_.end());
}
//_________________________________________________________________________________________________
double Channel::meanSampleValue() const
{
	if (sampleValues_.empty()) return 0;
	double total = std::accumulate(sampleValues_.begin(), sampleValues_.end(), 0.);
	return total/sampleValues_.size();
}
//_________________________________________________________________________________________________
size_t Channel::sampleCount() const
{
	return sampleValues_.size();
}
//_________________________________________________________________________________________________
const uint16_t* Channel::data() const
{
	if (sampleCount() > 0)
		return &sampleValues_.front();
	return 0;
}
//_________________________________________________________________________________________________
/**
 * Returns index of nearest Fixed Pattern Noise channel.
 */
size_t Channel::fpnNeighbourChannel() const
{
	size_t criterion[fpnCount];
	std::transform(fpnChannels, fpnChannels+fpnCount, criterion, std::bind2nd(distance_between_channels(), chanIdx_));
	size_t pos = std::distance(criterion, std::min_element(criterion, criterion+fpnCount));
	return fpnChannels[pos];
}
//_________________________________________________________________________________________________
} /* namespace get */
