/*
 * $Id: CoBoEvent.h 1659 2017-06-23 11:58:26Z psizun $
 * @file CoBoEvent.h
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

#ifndef get_CoBoEvent_h_INCLUDED
#define get_CoBoEvent_h_INCLUDED

#include "get/Channel.h"
#include <mfm/Common.h>
#include "utl/Exception.h"
#include <QMetaType>
#include <boost/cstdint.hpp>
#include <map>
#include <vector>
#include <fstream>

namespace mfm
{
	class Frame;
}

namespace get {
//_________________________________________________________________________________________________
typedef std::pair< uint8_t, uint8_t> ChannelKey; ///< Pair (channel index, AGet chip index).
typedef std::map< ChannelKey, get::Channel > ChannelMap;
//_________________________________________________________________________________________________
/**
 * A class to store all channel values from a CoBo physics event data frame, e.g. pertaining to an AsAd card.
 */
class CoBoEvent
{
public:
	class ChannelNotFound;
public:
	CoBoEvent(uint32_t eventIdx=0, uint8_t coboIdx=0, uint8_t asadIdx=0);
	void fromMatrix(std::ifstream &);
	bool fromFrame(mfm::Frame &);
	bool fromFrame_hardcoded(mfm::Frame & frame);
	Channel & channel(const uint8_t & chanIdx, const uint8_t & agetIdx=0);
	const Channel & channel(const uint8_t & chanIdx, const uint8_t & agetIdx=0) const;
	const ChannelMap & channels() const { return channels_; }
	bool containsAGetChip(const uint8_t & agetIdx=0) const;
	void clear();
	const uint64_t & eventTime() const;
	void setEventTime(const uint64_t & eventTime);
	const uint32_t & eventIdx() const;
	void setEventIdx(const uint32_t & eventIdx);
	uint8_t coboIdx() const;
	void setCoboIdx(const uint8_t & coboIdx);
	uint8_t asadIdx() const;
	void setAsadIdx(const uint8_t & asadIdx);
	const uint16_t & readOffset() const { return readOffset_; }
	void setReadOffset(const uint16_t & offset) { readOffset_ = offset; }
	const uint8_t & status() const { return status_; }
	void setStatus(const uint8_t & status) { status_ = status; }
	mfm::DynamicBitset & hitPattern(const uint8_t & agetIdx=0) { return hitPatterns_.at(agetIdx); }
	const mfm::DynamicBitset & hitPattern(const uint8_t & agetIdx=0) const { return hitPatterns_.at(agetIdx); }
	void setMultiplicity(const uint8_t & agetIdx, const uint16_t & mult);
	const uint16_t & multiplicity(const uint8_t & agetIdx=0) const;
	const uint32_t & windowOut() const { return windowOut_; }
	void setWindowOut(const uint32_t & windowOut) { windowOut_ = windowOut; }
	void setLastCell(const uint8_t & agetIdx, const uint16_t & lastCell);
	const uint16_t & lastCell(const uint8_t & agetIdx=0) const;
	size_t sampleCount() const;
private:
	void decodeSamples(mfm::Frame & frame);
	void decodeSamples_hardcoded(mfm::Frame & frame);
	struct add_channel_samples
	{
		size_t operator()(size_t n, const std::pair< ChannelKey, Channel > & c)
		{
			return n + c.second.sampleCount();
		}
	};
private:
	static const size_t AGET_PER_ASAD = 4u;
	uint64_t eventTime_; ///< Time offset of the frame as compared to start of run. Expressed in units of experiment clock.
	uint32_t eventIdx_; ///< The event to which this frame belongs within a run.
	uint8_t coboIdx_; ///< Index of CoBo subsystem from which the frame originates.
	uint8_t asadIdx_; ///< Index of AsAd board from which the frame data originates, relative to the CoBo subsystem to which it is connected.
	uint16_t readOffset_; ///< AGET readout offset.
	uint8_t status_; ///< Frame status info.
	std::vector< mfm::DynamicBitset > hitPatterns_; ///< Bit patterns indicating which of the 68 channels of AGET's have been hit.
	std::vector< uint16_t > multiplicities_; ///< AGET multiplicities.
	uint32_t windowOut_; ///< Output of the sliding window at the time of the event trigger.
	std::vector< uint16_t > lastCellIdxs_; ///< Last cell index read by AGET's.
	ChannelMap channels_; ///< Map of channels indexed by the channel index and the AGet chip index.
};
//_________________________________________________________________________________________________
class CoBoEvent::ChannelNotFound : public ::utl::Exception
{
public:
	ChannelNotFound(const uint8_t & chanIdx, const uint8_t & agetIdx);
	virtual void writeMessage(std::ostream& stream) const;
public:
	uint8_t chanIdx;
	uint8_t agetIdx;
};
//_________________________________________________________________________________________________
} /* namespace get */
Q_DECLARE_METATYPE(get::CoBoEvent)

#endif /* get_CoBoEvent_h_INCLUDED */
