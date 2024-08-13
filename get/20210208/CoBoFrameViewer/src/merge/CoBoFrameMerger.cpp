/**
 * @file CoBoFrameMerger.cpp
 * @date Oct 14, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: CoBoFrameMerger.cpp 1508 2015-10-15 16:24:29Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the CoBoFrameViewer software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
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
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "CoBoFrameMerger.h"
#include <mfm/Frame.h>
#include <mfm/Field.h>
#include <mfm/FrameDictionary.h>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include <cerrno>
#include <exception>
#include <fstream>
#include <set>
#include <string>

using namespace mfm;

namespace get
{
namespace data
{
//__________________________________________________________________________________________________
CoBoFrameMerger::CoBoFrameMerger()
{
}
//__________________________________________________________________________________________________
void CoBoFrameMerger::clear()
{
	streams_.clear();
	mergingMap_.clear();
}
//__________________________________________________________________________________________________
void CoBoFrameMerger::print() const
{
	LOG_DEBUG() << "Index contains " << mergingMap_.size() << " frame(s)";
}
//__________________________________________________________________________________________________
/**
 * Adds reference of each frame in given file to map sorted by event ID.
 */
bool CoBoFrameMerger::addFrames(const std::string & fileName)
{
	if (streams_.count(fileName)) return false;

	// Open file
	LOG_DEBUG() << "Indexing file '" << fileName << "'...";
	IStreamPtr streamPtr = IStreamPtr(new std::ifstream);
	streamPtr->exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		streamPtr->open(fileName.c_str(), std::ios::in | std::ios::binary);
	}
	catch (const std::ifstream::failure & e)
	{
		LOG_ERROR() << "Could not open file '" << fileName << "': " << strerror(errno);
		streams_.erase(fileName);
		return false;
	}

	// Loop over frames in input file
	std::auto_ptr<Frame> frame;
	try
	{
		do
		{
			try
			{
				// Read frame
				FrameRef ref;
				ref.frameStream = streamPtr;
				ref.framePos_B = streamPtr->tellg();
				frame = Frame::read(*streamPtr.get());

				// Skip frames with anything other than CoBo data
				const uint16_t frameTypeCode = frame->header().frameType();
				uint32_t eventID = 0;
				if (0x1 == frameTypeCode or 0x2 == frameTypeCode)
				{
					eventID = frame->headerField("eventIdx").value<uint32_t>();
				}
				else if (0x30 == frameTypeCode) // RIBF
				{
					// See http://wiki.ganil.fr/gap/browser/Documents/SPIRAL2_DAQ/DataFormat/GANIL_MFM_RIBF_data_format_v1.0.pdf
					eventID = frame->headerField(14, 4).value<uint32_t>();
				}
				else
				{
					LOG_DEBUG() << "Skipping frame of type " << frameTypeCode << "...";
					continue;
				}

				// Insert into map
				ref.frameSize_B = frame->header().frameSize_B();
				mergingMap_.insert( MergingMap::value_type(eventID, ref));
			}
			catch (const std::ifstream::failure & e)
			{
				if (streamPtr->rdstate() & std::ifstream::eofbit)
				{
					LOG_DEBUG() << "EOF reached.";
					break;
				}
				else
				{
					LOG_ERROR() << "Error reading frame: " << e.what();
				}
				return false;
			}
		}
		while (streamPtr->good());
		streamPtr->clear(); // clear fail and eof bits
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
		return false;
	}

	streams_[fileName] = streamPtr;
	//LOG_DEBUG() << "Map of frames now contains " << mergingMap_.size() << " frame(s)";
	return true;
}
//__________________________________________________________________________________________________
std::string CoBoFrameMerger::buildOutputFilename(const std::string & name, size_t & fileIndex)
{
	std::string actualName = name;

	if (fileIndex > 0)
	{
#if BOOST_FILESYSTEM_VERSION >= 3
		std::string oldext = fs::path(name).extension().string();
#else
		std::string oldext = fs::path(name).extension();
#endif
		std::ostringstream newext;
		newext << oldext << '.' << fileIndex;
		actualName = fs::path(name).replace_extension(newext.str()).string();
	}
	fileIndex++;
	return actualName;
}
//__________________________________________________________________________________________________
/**
 * Creates file with layered frames gathering frames sharing the same event identifier.
 * @param fileName Name of output file
 * @param minEventID Minimum event identifier
 * @param maxMergedEvents Maximum number of merged frames (or zero to merge all)
 */
void CoBoFrameMerger::mergeByEventId(const std::string & fileName, const uint32_t minEventID, const uint32_t maxEventID, const size_t maxMergedEvents, const uint64_t maxFileSize_B)
{
	const size_t numBasicFrames = mergingMap_.size();
	if (numBasicFrames <= 0) return;

	// Create output file
	size_t numFiles=0;
	std::ofstream ofs(buildOutputFilename(fileName, numFiles).c_str(), std::ios::out | std::ios::binary);
	if (not ofs.is_open())
	{
		LOG_FATAL() << "Could not create file '" << fileName << '\'';
		return;
	}
	const FrameFormat & layeredFormat = mfm::FrameDictionary::instance().findLatestFormat(0xFF01);
	const uint64_t layeredHeaderSize_B = layeredFormat.headerSize_blk()*layeredFormat.blkSize_B();

	// As long as map of frames is not empty
	size_t numLayeredFrames = 0;
	size_t numBasicFramesMerged = 0;
	const bool limitMergedEventCount = (maxMergedEvents > 0);
	std::set< size_t > frameCountsPerMergedEvent;
	std::set< size_t > mergedEventIDs;
	while (mergingMap_.size() > 0)
	{
		// Stop if user only wants to generate a few merged events
		if (limitMergedEventCount and numLayeredFrames >= maxMergedEvents)
			break;

		// Get smallest event ID
		const uint32_t  eventID = mergingMap_.lower_bound(0)->first;

		// Skip if event ID is too small or too large
		if (eventID >= minEventID and eventID <= maxEventID)
		{
			// Create new file if old one is too large
			size_t fileSize_B = ofs.tellp();
			if (maxFileSize_B > 0 and fileSize_B > maxFileSize_B)
			{
				ofs.close();
				ofs.open(buildOutputFilename(fileName, numFiles).c_str(), std::ios::out | std::ios::binary);
				if (not ofs.is_open())
				{
					LOG_FATAL() << "Could not create file '" << fileName << '\'';
					return;
				}
			}

			// Get iterators on frames with same event ID
			mergedEventIDs.insert(eventID);
			frameCountsPerMergedEvent.insert(mergingMap_.count(eventID));
			std::pair< MergingMap::iterator, MergingMap::iterator > frameRange = mergingMap_.equal_range(eventID);

			// Compute minimal size of layered frame
			uint64_t layeredFrameSize_B = layeredHeaderSize_B;
			for (MergingMap::iterator frameIter = frameRange.first; frameIter != frameRange.second; ++frameIter)
			{
				FrameRef ref = frameIter->second;
				layeredFrameSize_B += ref.frameSize_B;
			}

			// Create layered frame
			FrameFormat customFormat(layeredFormat);
			customFormat.set_blkSize_B(PrimaryHeader::minimal_blkSize_B(layeredFrameSize_B));
			//LOG_DEBUG() << "Creating layered frame of " << layeredFrameSize_B << " B with " << customFormat.blkSize_B() << " B blocks";
			std::auto_ptr<Frame> layeredFrame = customFormat.createFrame();

			// Loop on frames in range
			for (MergingMap::iterator frameIter = frameRange.first; frameIter != frameRange.second; ++frameIter)
			{
				// Read CoBo frame
				FrameRef ref = frameIter->second;
				ref.frameStream->seekg(ref.framePos_B, std::ios_base::beg);
				std::auto_ptr<Frame> coboFrame = mfm::Frame::read(*ref.frameStream.get());

				// Embed Cobo frame within layered frame
				layeredFrame->addFrame(*coboFrame);
				++numBasicFramesMerged;
			}

			// Add event ID to header of embedding frame
			layeredFrame->headerField("eventIdx").setValue(eventID);

			// Write layered frame to disk
			layeredFrame->write(ofs);

			numLayeredFrames++;
		}

		// Remove frames from list
		mergingMap_.erase(eventID);
	}

	// Close output file
	ofs.close();

	LOG_INFO() << "Gathered " << numBasicFramesMerged << " frame(s) out of " << numBasicFrames << " into " << numLayeredFrames << " layered frame(s) in file '" << fileName << "'";
	if (numBasicFramesMerged <= 0) return;
	LOG_INFO() << "The number of frames per layered frame ranges from " << *frameCountsPerMergedEvent.begin() << " to " << *frameCountsPerMergedEvent.rbegin();
	LOG_INFO() << "The event ID of merged frames varies from " << *mergedEventIDs.begin() << " to " << *mergedEventIDs.rbegin();
}
//__________________________________________________________________________________________________
} /* namespace data */
} /* namespace get */
