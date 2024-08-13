/*
 * $Id: CoBoEvent.cpp 1668 2017-07-06 13:13:19Z psizun $
 * @file CoBoEvent.cpp
 * @created 6 juin 2012
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

#include "get/CoBoEvent.h"
#include "utl/Logging.h"
#include <mfm/Exception.h>
#include <mfm/Frame.h>
#include <mfm/Field.h>
#include <mfm/BitField.h>
#include <utl/Logging.h>
#include <utl/BinIO.hpp>
#include <utl/BitFieldHelper.hpp>
using namespace mfm;
#include <iomanip>
#include <numeric>
#include <sstream>
#include <boost/cstdint.hpp>

namespace get {
//_________________________________________________________________________________________________
CoBoEvent::ChannelNotFound::ChannelNotFound(const uint8_t & chanIdx, const uint8_t & agetIdx)
	: chanIdx(chanIdx), agetIdx(agetIdx)
{
}
//_________________________________________________________________________________________________
void CoBoEvent::ChannelNotFound::writeMessage(std::ostream & stream) const
{
	stream << "Could not find channel (chanIdx=" << (short) chanIdx << ", agetIdx=" << (short) agetIdx << ")";
}
//_________________________________________________________________________________________________
CoBoEvent::CoBoEvent(uint32_t eventIdx, uint8_t coboIdx, uint8_t asadIdx)
	: eventTime_(0), eventIdx_(eventIdx), coboIdx_(coboIdx), asadIdx_(asadIdx),
	  readOffset_(0), status_(0),
	  hitPatterns_(AGET_PER_ASAD, mfm::DynamicBitset(72)),
	  multiplicities_(AGET_PER_ASAD, 0),
	  windowOut_(0),
	  lastCellIdxs_(AGET_PER_ASAD, 0)
{
}
//_________________________________________________________________________________________________
void CoBoEvent::fromMatrix(std::ifstream & in)
{
	const size_t maxAsic = 4u;
	const size_t maxChannels = 68u;
	const size_t maxBuckets = 512u;

	clear();
	double sampleValues[maxBuckets];
	for (size_t asicIdx=0; asicIdx < maxAsic; ++asicIdx)
	{
		for (size_t chanIdx=0; chanIdx < maxChannels; ++chanIdx)
		{
			in.read(reinterpret_cast<char*>(&sampleValues), maxBuckets*sizeof(double));
			for (size_t buckIdx=0; buckIdx < maxBuckets; ++buckIdx)
			{
				channel(chanIdx, asicIdx).addValue(buckIdx, sampleValues[buckIdx]);
			}
		}
	}
}
//_________________________________________________________________________________________________
bool CoBoEvent::fromFrame(mfm::Frame & frame)
{
	clear();

	try {
		// Get meta-data
		setEventTime(frame.headerField("eventTime").value<uint64_t>());
		setEventIdx(frame.headerField("eventIdx").value<uint32_t>());
		setCoboIdx(frame.headerField("coboIdx").value<uint8_t>());
		setAsadIdx(frame.headerField("asadIdx").value<uint8_t>());
		setReadOffset(frame.headerField("readOffset").value<uint16_t>());
		setStatus(frame.headerField("status").value<uint8_t>());

		// Get hit patterns
		for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
		{
			std::ostringstream patternName;
			patternName << "hitPat_" << agetIdx;
			frame.headerField(patternName.str()).value(hitPattern(agetIdx));
		}
		// Get multiplicities
		for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
		{
			std::ostringstream patternName;
			patternName << "multip_" << agetIdx;
			setMultiplicity(agetIdx, frame.headerField(patternName.str()).value<uint16_t>());
		}

		// Get sliding window (revision 3)
		try
		{
			setWindowOut(frame.headerField("windowOut").value<uint32_t>());
		}
		catch (const mfm::FieldNotFound &)
		{
			setWindowOut(0);
		}

		// Get last cell indices (revision 4)
		for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
		{
			try
			{
				std::ostringstream patternName;
				patternName << "lastCell_" << agetIdx;
				setLastCell(agetIdx, frame.headerField(patternName.str()).value<uint16_t>());
			}
			catch (const mfm::FieldNotFound &)
			{
				setLastCell(agetIdx, 0);
			}
		}

		// Decode samples
		decodeSamples(frame);
	}
	catch (const std::ios_base::failure & e)
	{
		LOG_ERROR() << "Error decoding CoBo data frame (" << e.what() << ")";
		return false;
	}
	return true;
}
//_________________________________________________________________________________________________
void CoBoEvent::decodeSamples(mfm::Frame & frame)
{
	const uint16_t frameType = frame.header().frameType();
	const size_t itemCount = frame.itemCount();
	size_t numSamples = itemCount;

	// Check consistency of frame header
	const size_t actualItemCount = (frame.header().frameSize_B() - frame.header().headerSize_B())
			/ frame.header().itemSize_B();
	if (actualItemCount < itemCount) {
		LOG_ERROR() << "Data frame size is not consistent with sample count!";
		numSamples = actualItemCount;
	}

	// Decode uncompressed CoBo frame
	if (1u == frameType)
	{
		mfm::BitField agetIdxField, chanIdxField, buckIdxField, sampleValueField;
		mfm::Field field;
		// Loop through items
		// Decode first item
		if (numSamples > 0)
		{
			mfm::Item item = frame.itemAt(0u);

			field = item.field("");
			agetIdxField = field.bitField("agetIdx");
			chanIdxField = field.bitField("chanIdx");
			buckIdxField = field.bitField("buckIdx");
			sampleValueField = field.bitField("sample");

			const uint32_t agetIdx = agetIdxField.value<uint32_t>();
			const uint32_t chanIdx = chanIdxField.value<uint32_t>();
			const uint32_t buckIdx = buckIdxField.value<uint32_t>();
			const uint32_t sampleValue = sampleValueField.value<uint32_t>();
			channel(chanIdx, agetIdx).addValue(buckIdx, sampleValue);
		}
		// Loop through other items
		for (size_t itemId=1; itemId < numSamples; ++itemId)
		{
			mfm::Item item = frame.itemAt(itemId);

			field = item.field(field);
			agetIdxField = field.bitField(agetIdxField);
			chanIdxField = field.bitField(chanIdxField);
			buckIdxField = field.bitField(buckIdxField);
			sampleValueField = field.bitField(sampleValueField);

			const uint32_t agetIdx = agetIdxField.value<uint32_t>();
			const uint32_t chanIdx = chanIdxField.value<uint32_t>();
			const uint32_t buckIdx = buckIdxField.value<uint32_t>();
			const uint32_t sampleValue = sampleValueField.value<uint32_t>();
			//LOG_DEBUG() << "item=" << std::setw(6) << itemId << " aget=" << agetIdx << " chan=" << std::setw(2) << chanIdx << " bucket=" << std::setw(3) << buckIdx << " sample=" << sampleValue;
			channel(chanIdx, agetIdx).addValue(buckIdx, sampleValue);
		}
	}
	// Decode compact frame (without bucket index and channel index)
	else if (2u == frameType)
	{
		// In full readout mode, every channel and time bucket is read from the AGET and placed in memory in the order received.
		// Since the order the data is the same for each event, the channel index and bucket index are omitted to reduce the item size.
		// In the current firmware, both samples in an item will be from the same AGET, so agetIdx1 and agetIdx2 will always match.
		const size_t numChannels = 68u;
		const size_t numChips = 4u;
		mfm::BitField agetIdxField, sampleValueField;
		mfm::Field field;

		// Loop through items
		std::vector< uint32_t > chanIdx(numChips, 0u);
		std::vector< uint32_t > buckIdx(numChips, 0u);
		// Decode first item
		if (numSamples > 0)
		{
			mfm::Item item = frame.itemAt(0u);

			field = item.field("");
			agetIdxField = field.bitField("agetIdx");
			sampleValueField = field.bitField("sample");

			const uint32_t agetIdx = agetIdxField.value<uint32_t>();
			const uint32_t sampleValue = sampleValueField.value<uint32_t>();
			channel(chanIdx[agetIdx], agetIdx).addValue(buckIdx[agetIdx], sampleValue);
			chanIdx[agetIdx]++;
		}
		// Loop through other items
		for (size_t itemId=1; itemId < numSamples; ++itemId)
		{
			mfm::Item item = frame.itemAt(itemId);

			field = item.field(field);
			agetIdxField = field.bitField(agetIdxField);
			sampleValueField = field.bitField(sampleValueField);

			const uint32_t agetIdx = agetIdxField.value<uint32_t>();
			if (chanIdx[agetIdx] >= numChannels)
			{
				chanIdx[agetIdx] = 0u;
				buckIdx[agetIdx]++;
			}
			const uint32_t sampleValue = sampleValueField.value<uint32_t>();
			channel(chanIdx[agetIdx], agetIdx).addValue(buckIdx[agetIdx], sampleValue);
			chanIdx[agetIdx]++;
		}
	}

}
//_________________________________________________________________________________________________
bool CoBoEvent::fromFrame_hardcoded(mfm::Frame & frame)
{
	clear();

	try {
		// Get meta-data
		setEventTime(frame.headerField(16u, 6u).value<uint64_t>());
		setEventIdx(frame.headerField(22u, 4u).value<uint32_t>());
		setCoboIdx(frame.headerField(26u, 1u).value<uint8_t>());
		setAsadIdx(frame.headerField(27u, 1u).value<uint8_t>());
		setReadOffset(frame.headerField(28u, 2u).value<uint16_t>());
		setStatus(frame.headerField(30u, 1u).value<uint8_t>());

		// Get hit patterns
		for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
		{
			frame.headerField(31u + 9u*agetIdx, 9u).value(hitPattern(agetIdx));
		}
		// Get multiplicities
		for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
		{
			setMultiplicity(agetIdx, frame.headerField(67u + 2u*agetIdx, 2u).value<uint16_t>());
		}

		const uint8_t formatRevision = frame.header().revision();
		if (formatRevision >= 3u)
		{
			// Get sliding window (revision 3)
			setWindowOut(frame.headerField(75u, 4u).value< uint32_t >());
		}
		else
		{
			setWindowOut(0);
		}

		// Get last cell indices (revision 4)
		if (formatRevision >= 4u)
		{
			for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
			{
				setLastCell(agetIdx, frame.headerField(79u + agetIdx*2u, 2u).value< uint16_t >());
			}
		}
		else
		{
			for (short agetIdx = 0; agetIdx < 4; ++agetIdx)
			{
				setLastCell(agetIdx, 0);
			}
		}

		// Decode samples
		decodeSamples_hardcoded(frame);
	}
	catch (const std::ios_base::failure & e)
	{
		LOG_ERROR() << "Error decoding CoBo data frame (" << e.what() << ")";
		return false;
	}
	return true;
}
//_________________________________________________________________________________________________
void CoBoEvent::decodeSamples_hardcoded(mfm::Frame & frame)
{
	const size_t numSamples = frame.itemCount();
	if (numSamples == 0) return;
	const uint16_t frameType = frame.header().frameType();

	const uint32_t dataOffset_B = frame.header().headerSize_B();
	// Get stream
	std::istream & stream = frame.serializer().inputStream(dataOffset_B);

	// Decode uncompressed CoBo frame
	if (1u == frameType)
	{
		// Loop through items
		uint32_t field = 0;
		for (size_t itemId=0; itemId < numSamples; ++itemId)
		{
			// Read item
			stream.read(reinterpret_cast<char *>(&field), 4lu);
#if WORDS_LITTLE_ENDIAN // if host is little endian
			utl::BinIO::swap(field);
#endif

			// Decode bit fields
			const uint32_t agetIdx = utl::BitFieldHelper< uint32_t >::getField(field, 30u, 2u);
			const uint32_t chanIdx = utl::BitFieldHelper< uint32_t >::getField(field, 23u, 7u);
			const uint32_t buckIdx = utl::BitFieldHelper< uint32_t >::getField(field, 14u, 9u);
			const uint32_t sampleValue = utl::BitFieldHelper< uint32_t >::getField(field, 0u, 12u);
			channel(chanIdx, agetIdx).addValue(buckIdx, sampleValue);
		}
	}
	// Decode compact frame (without bucket index and channel index)
	else if (2u == frameType)
	{
		// Loop through items
		uint16_t field = 0;
		std::vector< uint32_t > chanIdx(4u, 0u);
		std::vector< uint32_t > buckIdx(4u, 0u);
		const size_t numChannels = 68u;
		for (size_t itemId=0; itemId < numSamples; ++itemId)
		{
			// Read item
			stream.read(reinterpret_cast<char *>(&field), 2lu);
#if WORDS_LITTLE_ENDIAN // if host is little endian
			utl::BinIO::swap(field);
#endif
			// Decode bit fields
			const uint32_t agetIdx = utl::BitFieldHelper< uint16_t >::getField(field, 14u, 1u);
			if (chanIdx[agetIdx] >= numChannels)
			{
				chanIdx[agetIdx] = 0u;
				buckIdx[agetIdx]++;
			}
			const uint32_t sampleValue = utl::BitFieldHelper< uint16_t >::getField(field, 0u, 12u);
			channel(chanIdx[agetIdx], agetIdx).addValue(buckIdx[agetIdx], sampleValue);
			chanIdx[agetIdx]++;
		}
	}
}
//_________________________________________________________________________________________________
bool CoBoEvent::containsAGetChip(const uint8_t & agetIdx) const
{
	ChannelMap::const_iterator it;
	for (it = channels_.begin(); it != channels_.end(); ++it)
	{
		if ((it->first).second == agetIdx) return true;
	}
	return false;
}
//_________________________________________________________________________________________________
Channel & CoBoEvent::channel(const uint8_t & chanIdx, const uint8_t & agetIdx)
{
	ChannelKey key(chanIdx, agetIdx);
	ChannelMap::iterator it = channels_.find(key);
	if (it == channels_.end())
	{
		channels_[key].reserve(512);
		channels_[key].setChanIdx(chanIdx);
		channels_[key].setAgetIdx(agetIdx);
	}
	return channels_[key];
}
//_________________________________________________________________________________________________
const Channel & CoBoEvent::channel(const uint8_t & chanIdx, const uint8_t & agetIdx) const
{
	ChannelKey key(chanIdx, agetIdx);
	ChannelMap::const_iterator it = channels_.find(key);
	if (it == channels_.end())
	{
		throw ChannelNotFound(chanIdx, agetIdx);
	}
	return channels_.at(key);
}
//_________________________________________________________________________________________________
void CoBoEvent::clear()
{
	channels_.clear();
	for (size_t agetIdx=0; agetIdx < AGET_PER_ASAD; ++agetIdx)
		hitPatterns_[agetIdx].reset();
}
//_________________________________________________________________________________________________
uint8_t CoBoEvent::asadIdx() const
{
	return asadIdx_;
}
//_________________________________________________________________________________________________
void CoBoEvent::setAsadIdx(const uint8_t & asadIdx)
{
	asadIdx_ = asadIdx;
}
//_________________________________________________________________________________________________
uint8_t CoBoEvent::coboIdx() const
{
	return coboIdx_;
}
//_________________________________________________________________________________________________
void CoBoEvent::setCoboIdx(const uint8_t & coboIdx)
{
	coboIdx_ = coboIdx;
}
//_________________________________________________________________________________________________
const uint64_t & CoBoEvent::eventTime() const
{
	return eventTime_;
}
//_________________________________________________________________________________________________
void CoBoEvent::setEventTime(const uint64_t & eventTime)
{
	eventTime_ = eventTime;
}
//_________________________________________________________________________________________________
const uint32_t & CoBoEvent::eventIdx() const
{
	return eventIdx_;
}
//_________________________________________________________________________________________________
void CoBoEvent::setEventIdx(const uint32_t & eventIdx)
{
	eventIdx_ = eventIdx;
}
//_________________________________________________________________________________________________
void CoBoEvent::setMultiplicity(const uint8_t & agetIdx, const uint16_t & mult)
{
	multiplicities_[agetIdx] = mult;
}
//_________________________________________________________________________________________________
const uint16_t & CoBoEvent::multiplicity(const uint8_t & agetIdx) const
{
	return multiplicities_[agetIdx];
}
//_________________________________________________________________________________________________
void CoBoEvent::setLastCell(const uint8_t & agetIdx, const uint16_t & lastCell)
{
	lastCellIdxs_[agetIdx] = lastCell;
}
//_________________________________________________________________________________________________
const uint16_t & CoBoEvent::lastCell(const uint8_t & agetIdx) const
{
	return lastCellIdxs_[agetIdx];
}
//_________________________________________________________________________________________________
size_t CoBoEvent::sampleCount() const
{
	return std::accumulate(channels_.begin(), channels_.end(), 0, add_channel_samples());
}
//_________________________________________________________________________________________________
} /* namespace get */
