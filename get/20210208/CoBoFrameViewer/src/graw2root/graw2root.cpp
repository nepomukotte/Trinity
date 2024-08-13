/*
 * $Id$
 * @file graw2root.cpp
 * @created 10 juil. 2012
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

#include "graw2root.h"

#include <mfm/BitField.h>
#include <mfm/Field.h>
#include <mfm/Frame.h>
#include <mfm/Exception.h>
#include <utl/Logging.h>
#include "GDataFrame.h"
#include "GDataChannel.h"
#include "GDataSample.h"
#include "T2KEvent.h"
#include "get/CoBoEvent.h"
#include "get/Channel.h"

#include <TFile.h>
#include <TTree.h>
#include <TObjString.h>
#include <RVersion.h> // ROOT_VERSION_CODE

//#include <boost/timer.hpp>

#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>

using mfm::Frame;
using std::strerror;

#if ROOT_VERSION_CODE > 393216
using namespace GET; // ROOT 6
#else
using namespace get;
#endif

namespace get {
//_________________________________________________________________________________________________
/**
 * Embed configuration into UserInfo
 */
bool add_config_to_tree(const std::string & filename, TTree* tree)
{
	if (filename.empty()) return true;

	LOG_INFO() << "Embedding file '" << filename << "' into ROOT tree...";
	try
	{
		//  Read file into a string
		std::string contents;
		std::ifstream file(filename.c_str(), std::ios::in);
		if (! file.is_open())
		{
			LOG_ERROR() << "File '" << filename << "' does not exist!";
			return false;
		}
		{
			file.seekg(0, std::ios::end);
			contents.resize(file.tellg());
			file.seekg(0, std::ios::beg);
			file.read(&contents[0], contents.size());
			file.close();
		}
		if (contents.empty()) return true;

		TObjString* cfgObj = new TObjString(contents.c_str());
		tree->GetUserInfo()->Add(cfgObj);
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Failed to embed '" << filename << "' into TTree. " << e.what();
		return false;
	}
	return true;
}
//_________________________________________________________________________________________________
bool graw2root(const std::string & input_filename, const std::string & output_filename, const size_t maxFrames, const size_t frameOffset, const std::string & cfg)
{
	LOG_INFO() << "Converting frame contents from raw items to ROOT tree...";
	LOG_DEBUG() << "Opening file '" << input_filename << "'...";
	std::ifstream in;
	in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		in.open(input_filename.c_str(), std::ios::in | std::ios::binary);
	}
	catch (const std::ifstream::failure & e)
	{
		LOG_ERROR() << "Could not open file '" << input_filename << "': " << strerror(errno);
		return false;
	}

	// Create output file
	LOG_INFO() << "Creating file '" << output_filename << "'...";
	TFile* out = new TFile(output_filename.c_str(), "RECREATE");

	// Create tree
	TTree* tree = new TTree("tree", "A tree of GET data frames");
	GDataFrame dataFrame;
	tree->Branch("GDataFrame", &dataFrame);

	// Embed configuration into UserInfo
	add_config_to_tree(cfg, tree);

	// Loop over frames in input file
	std::auto_ptr<Frame> frame;
	Frame::seekFrame(in, frameOffset);
	size_t frameCount = 0;
	get::CoBoEvent event;
	try {
		do
		{
			if (maxFrames > 0 and frameCount >= maxFrames) break;
			//boost::timer timer;
			// Read frame
			try
			{
				frame = Frame::read(in);
			}
			catch (const std::ifstream::failure & e)
			{
				if (in.rdstate() & std::ifstream::eofbit)
				{
					LOG_WARN() << "EOF reached.";
					break;
				}
				else
				{
					LOG_ERROR() << "Error reading frame: " << e.what();
				}
				return false;
			}

			frameCount++;
			// Skip frames with anything other than CoBo data
			if (0x1 != frame->header().frameType() and 0x2 != frame->header().frameType())
			{
				LOG_DEBUG() << "Skipping frame #" << (frameOffset + frameCount) << "...";
				continue;
			}
			else if (frameCount%50 == 1)
			{
				LOG_DEBUG() << "Adding frame #" << (frameOffset + frameCount) << " to tree...";
			}

			// Decode frame
			event.fromFrame(*frame.get());
			//event.fromFrame_hardcoded(*frame.get());

			// Reset ROOT frame
			dataFrame.Clear();
			// Save current Object count
			Int_t objectCount = TProcessID::GetObjectCount();

			// Get meta-data
			dataFrame.fHeader.fRevision = frame->header().revision();
			dataFrame.fHeader.fDataSource = frame->header().dataSource();
			dataFrame.fHeader.fEventTime = event.eventTime();
			dataFrame.fHeader.fEventIdx = event.eventIdx();
			dataFrame.fHeader.fCoboIdx = event.coboIdx();
			dataFrame.fHeader.fAsadIdx = event.asadIdx();
			dataFrame.fHeader.fReadOffset = event.readOffset();
			dataFrame.fHeader.fStatus = event.status();
			for (size_t agetId=0; agetId < 4u; ++agetId)
			{
				const mfm::DynamicBitset & pattern = event.hitPattern(agetId);
				for (size_t chanId=0; chanId < GFrameHeader::MAX_CHANNELS; ++chanId)
				{
					dataFrame.fHeader.fHitPatterns[agetId].SetBitNumber(GFrameHeader::MAX_CHANNELS - 1u - chanId, pattern[chanId]);
				}
				dataFrame.fHeader.fMult[agetId] = event.multiplicity(agetId);
			}
			dataFrame.fHeader.fWindowOut = event.windowOut();
			for (size_t agetId=0; agetId < 4u; ++agetId)
			{
				dataFrame.fHeader.fLastCellIdx[agetId] = event.lastCell(agetId);
			}

			// Loop over items
			const uint32_t maxAget = 4u;
			const size_t maxChanPerAget = 68u;
			for (size_t chanIdx=0; chanIdx < maxChanPerAget; ++chanIdx)
			{
				for (size_t agetIdx=0; agetIdx < maxAget; ++agetIdx)
				{
					try
					{
						// The const version of CoBoEvent::channel throws an exception for non-existing channels
						// while the non-const version creates the channel
						const get::CoBoEvent & constEvent = event;
						const get::Channel & channel = constEvent.channel(chanIdx, agetIdx);

						GDataChannel* chan = dataFrame.FindChannel(agetIdx, chanIdx);

						const std::vector< uint16_t > & bucketIndexes = channel.buckIndexes();
						const std::vector< uint16_t > & values = channel.sampleValues();
						const size_t numBuckets = channel.sampleCount();
						for (size_t buckIdx=0; buckIdx < numBuckets; ++buckIdx)
						{
							GDataSample* sample = dataFrame.AddSample();
							sample->Set(bucketIndexes[buckIdx], values[buckIdx]);
							chan->AddSample(sample);
						};
					}
					catch (const get::CoBoEvent::ChannelNotFound &)
					{
					}
				}
			}
			//LOG_DEBUG() << "fNchannels=" << dataFrame.GetNchannels();

			//Restore Object count (see JetEvent class in tutorials)
			//To save space in the table keeping track of all referenced objects
			//we assume that our events do not address each other. We reset the
			//object count to what it was at the beginning of the event.
			TProcessID::SetObjectCount(objectCount);

			// Fill tree
			tree->Fill();
			//LOG_DEBUG() << "Frame converted in " << timer.elapsed() << " s";
		}
		while (in.good());
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}

	// Close input file
	in.close();

	// Close output file
	LOG_INFO() << "Saving ROOT file";
	tree->Write();
	LOG_DEBUG() << "ROOT tree contains " << tree->GetEntries() << " event(s)";
	out->Close();
	delete out;

	return true;
}
//_________________________________________________________________________________________________
bool graw2t2k(const std::string & input_filename, const std::string & output_filename, const size_t maxFrames, const size_t frameOffset, const std::string & cfg)
{
	LOG_INFO() << "Converting frame contents from raw items to ROOT tree of T2K events...";
	LOG_DEBUG() << "Opening file '" << input_filename << "'...";
	std::ifstream in;
	in.exceptions(std::ifstream::failbit | std::ifstream::badbit);
	try
	{
		in.open(input_filename.c_str(), std::ios::in | std::ios::binary);
	}
	catch (const std::ifstream::failure & e)
	{
		LOG_ERROR() << "Could not open file '" << input_filename << "': " << strerror(errno);
		return false;
	}

	// Create output file
	LOG_INFO() << "Creating file '" << output_filename << "'...";
	TFile* out = new TFile(output_filename.c_str(), "RECREATE");
	// Create tree
	TTree* tree = new TTree("TT2K", "A tree with T2K events.");
	T2KEvent t2kEvent;
	tree->Branch("T2KEvent", &t2kEvent);

	// Embed configuration into UserInfo
	if (not add_config_to_tree(cfg, tree))
		return false;

	// Loop over frames in input file
	std::auto_ptr<Frame> frame;
	Frame::seekFrame(in, frameOffset);
	size_t frameCount = 0;
	get::CoBoEvent event;
	try {
		do
		{
			if (maxFrames > 0 and frameCount >= maxFrames) break;
			// Read frame
			//boost::timer timer;
			try
			{
				//LOG_INFO() << "Reading frame...";
				frame = Frame::read(in);
			}
			catch (const std::ifstream::failure & e)
			{
				if (in.rdstate() & std::ifstream::eofbit)
				{
					LOG_WARN() << "EOF reached.";
					break;
				}
				else
				{
					LOG_ERROR() << "Error reading frame: " << e.what();
				}
				return false;
			}

			frameCount++;
			// Skip frames with anything other than CoBo data
			if (0x1 != frame->header().frameType() and 0x2 != frame->header().frameType())
			{
				LOG_DEBUG() << "Skipping frame #" << (frameOffset + frameCount) << "...";
				continue;
			}

			//LOG_DEBUG() << "Decoding frame...";
			//event.fromFrame(*frame.get());
			event.fromFrame_hardcoded(*frame.get());

			const uint32_t maxAget = 4u;
			const size_t maxChanPerAget = 68u;
			const size_t maxTimeBuckets = 512;
			t2kEvent.SetfNChannels(maxChanPerAget*maxAget);
			t2kEvent.SetfNBins(maxTimeBuckets);
			t2kEvent.SetfNTBins();

			// Reset event contents (in case not all values are provided)
			for (size_t chanIdx=0; chanIdx < maxChanPerAget; ++chanIdx)
				for (size_t agetIdx=0; agetIdx < maxAget; ++agetIdx)
					for (size_t buckIdx=0; buckIdx < maxTimeBuckets; ++buckIdx)
					{
						t2kEvent.SetPoint(agetIdx*maxChanPerAget+chanIdx, buckIdx, 0);
					}

			for (size_t chanIdx=0; chanIdx < maxChanPerAget; ++chanIdx)
			{
				for (size_t agetIdx=0; agetIdx < maxAget; ++agetIdx)
				{
					try
					{
						const get::Channel & channel = event.channel(chanIdx, agetIdx);
						const std::vector< uint16_t > & bucketIndexes = channel.buckIndexes();
						const std::vector< uint16_t > & values = channel.sampleValues();
						const size_t numBuckets = bucketIndexes.size();
						for (size_t buckIdx=0; buckIdx < numBuckets; ++buckIdx)
						{
							t2kEvent.SetPoint(agetIdx*maxChanPerAget+chanIdx,
									bucketIndexes[buckIdx], values[buckIdx]);

						};
					}
					catch (const get::CoBoEvent::ChannelNotFound &)
					{
					}
				}
			}

			// Fill tree
			frameCount++;
			if (frameCount%50 == 1)
			{
				LOG_DEBUG() << "Adding frame #" << (frameOffset + frameCount) << " to tree...";
			}
			tree->Fill();
			//LOG_DEBUG() << "Frame converted in " << timer.elapsed() << " s";
		}
		while (in.good());
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Error: " << e.what();
	}

	// Close input file
	in.close();

	// Close output file
	LOG_INFO() << "Saving ROOT file";
	tree->Write();
	out->Close();
	delete out;

	return true;
}
//_________________________________________________________________________________________________
}; // namespace get
