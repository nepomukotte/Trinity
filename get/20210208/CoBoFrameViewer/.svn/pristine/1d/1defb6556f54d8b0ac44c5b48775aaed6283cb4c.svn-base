/**
 * @file GDataFrame.cpp
 * @date 26 févr. 2013
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the CoBoFrameViewer software project.
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

#include "GDataFrame.h"
#include "GDataChannel.h"
#include "GDataSample.h"

#include <TClonesArray.h>
#include <TCollection.h>



#if ROOT_VERSION_CODE > 393216
ClassImp(GET::GDataFrame)
namespace GET {
#else
ClassImp(get::GDataFrame)
namespace get {
#endif
//_________________________________________________________________________________________________
TClonesArray* GDataFrame::fgChannels   = 0;
TClonesArray* GDataFrame::fgSamples   = 0;
//_________________________________________________________________________________________________
GDataFrame::GDataFrame() : fNchannels(0), fNsamples(0)
{
	// Create a GDataFrame object.
	// When the constructor is invoked for the first time, the class static
	// variables fgxxx are 0 and the TClonesArray fgxxx are created.
#if ROOT_VERSION_CODE > 393216
	const char* nsName = "GET";
#else
	const char* nsName = "get";
#endif
	if (!fgChannels) fgChannels = new TClonesArray(TString::Format("%s::GDataChannel", nsName), 2048);
	fChannels = fgChannels;
	if (!fgSamples) fgSamples = new TClonesArray(TString::Format("%s::GDataSample", nsName), 139264);
	fSamples = fgSamples;
}
//_________________________________________________________________________________________________
GDataFrame::~GDataFrame()
{
	Reset();
}
//_________________________________________________________________________________________________
GDataChannel* GDataFrame::AddChannel(const UShort_t & agetIdx, const UShort_t & chanIdx)
{
	// Add a new channel to the list of channels for this frame.
	TClonesArray &channels = *fChannels;
	GDataChannel* channel = new(channels[fNchannels++]) GDataChannel(agetIdx, chanIdx);
	return channel;
}
//_________________________________________________________________________________________________
GDataChannel* GDataFrame::SearchChannel(const UShort_t & agetIdx, const UShort_t & chanIdx)
{
	// Searches for a channel in the list of channels for this frame.
	GDataChannel* channel = 0;
	TIter iter(fChannels);
	while ((channel = (GDataChannel*) iter.Next()))
	{
		if (channel->fAgetIdx == agetIdx && channel->fChanIdx == chanIdx) return channel;
	}
	return channel;
}
//_________________________________________________________________________________________________
GDataChannel* GDataFrame::FindChannel(const UShort_t & agetIdx, const UShort_t & chanIdx)
{
	// Searches for a channel or creates a new one
	GDataChannel* channel = SearchChannel(agetIdx, chanIdx);
	if (! channel)
	{
		channel = AddChannel(agetIdx, chanIdx);
	}
	return channel;
}
//_________________________________________________________________________________________________
GDataSample* GDataFrame::AddSample()
{
	// Add a new sample to the list of samples for this frame.
	TClonesArray &samples = *fSamples;
	GDataSample* sample = new(samples[fNsamples++]) GDataSample();
	return sample;
}
//_________________________________________________________________________________________________
GDataSample* GDataFrame::AddSample(const UShort_t & agetIdx, const UShort_t & chanIdx, const UShort_t & buckIdx, const UShort_t & value)
{
	GDataChannel* channel = FindChannel(agetIdx, chanIdx);

	GDataSample* sample = AddSample();
	sample->Set(buckIdx, value);

	channel->AddSample(sample);

	return sample;
}
//_________________________________________________________________________________________________
void GDataFrame::Clear(Option_t *option)
{
	fHeader.Clear(option);
	fChannels->Clear(option);
	fNchannels = 0;
	fSamples->Clear(option);
	fNsamples = 0;
}
//_________________________________________________________________________________________________
void GDataFrame::Reset(Option_t * /* option */)
{
	// Static function to reset all static objects for this event
	delete fgChannels;
	fgChannels = 0;
	delete fgSamples;
	fgSamples = 0;
}
//_________________________________________________________________________________________________
} /* namespace get */
