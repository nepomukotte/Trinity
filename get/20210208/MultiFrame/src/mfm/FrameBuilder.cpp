/*
 * $Id: FrameBuilder.cpp 1568 2014-04-03 11:57:52Z psizun $
 * @file FrameBuilder.cpp
 * @created 13 avr. 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * -----------------------------------------------------------------------------
 */

#include "mfm/FrameBuilder.h"
#include "mfm/PrimaryHeader.h"
#include "mfm/Frame.h"
#include "utl/Logging.h"
#include <iomanip>

namespace mfm {
//______________________________________________________________________
/**
 * Default constructor.
 */
FrameBuilder::FrameBuilder()
	: serializer_(0), frameSize_(0)
{
}
//______________________________________________________________________
FrameBuilder::~FrameBuilder()
{
	reset();
}
//______________________________________________________________________
/**
 * Appends a chunk of data to the frame being built.
 */
void FrameBuilder::addDataChunk(const mfm::Byte* begin, const mfm::Byte* end)
{
	size_t const chunkSize_B = end - begin;
	LOG_DEBUG() << "Adding chunk of " << chunkSize_B << " B to frame builder.";

	// Resize buffer
	size_t const oldCapacity = serializer_.capacity();
	serializer_.setCapacity(oldCapacity + chunkSize_B);
	serializer_.set_size_B(oldCapacity + chunkSize_B);

	// Add new data chunk to buffer
	serializer_.read(chunkSize_B, begin, oldCapacity);

	// Build complete frames
	buildFrames(end);
}
//______________________________________________________________________
void FrameBuilder::buildFrames(const mfm::Byte* end)
{
	while (true)
	{
		// Check if primary header can be reconstructed
		if (serializer_.size_B() < mfm::PrimaryHeader::SPEC_SIZE_B)
			break;

		// Decode primary header and get frame size
		if (frameSize_ <= 0)
		{
			std::auto_ptr<mfm::PrimaryHeader> headerPtr = mfm::PrimaryHeader::decodePrimaryHeader(serializer_.inputStream());
			frameSize_ = headerPtr->frameSize_B();
			processHeader(*headerPtr);
			LOG_DEBUG() << "Expecting frame of " << frameSize_ << " B.";
		}

		// Check if frame is complete
		if (serializer_.size_B() < frameSize_)
		{
			LOG_DEBUG() << "Still missing " << frameSize_ - serializer_.size_B() << " B out of " << frameSize_ << " B";
			break;
		}

		// Reconstruct frame
		mfm::Frame frame(mfm::Serializer(serializer_, frameSize_));

		// Process frame
		LOG_DEBUG() << "Processing frame of " << frameSize_ << " B.";
		processFrame(frame);

		// Remove frame data from buffer
		size_t const newCapacity = serializer_.size_B() - frameSize_;
		serializer_.setCapacity(newCapacity);
		serializer_.set_size_B(newCapacity);
		frameSize_ = 0;
		if (newCapacity > 0)
		{
			LOG_DEBUG() << "Frame builder contains " << newCapacity << " B left.";
			serializer_.read(newCapacity, end - newCapacity);
		}
	}
	return;
}
//______________________________________________________________________
/**
 * @fn void FrameBuilder::processFrame(mfm::Frame & frame)
 * Processes the complete frame.
 * @param frame Complete reconstructed frame.
 * This method should be implemented in derived classes.
 */
//______________________________________________________________________
/**
 * Clears buffer to be ready to build a new frame from scratch.
 */
void FrameBuilder::reset()
{
	LOG_DEBUG() << "Resetting frame builder.";
	serializer_.setCapacity(0u);
	serializer_.set_size_B(0u);
	frameSize_ = 0;
}
//______________________________________________________________________
} /* namespace mfm */
