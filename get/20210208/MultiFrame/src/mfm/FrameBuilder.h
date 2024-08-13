/*
 * $Id: FrameBuilder.h 1568 2014-04-03 11:57:52Z psizun $
 * @file FrameBuilder.h
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

#ifndef mfm_FrameBuilder_h_INCLUDED
#define mfm_FrameBuilder_h_INCLUDED

#include "mfm/Serializer.h"

namespace mfm {
// Forward declarations
class Frame;
class PrimaryHeader;
//______________________________________________________________________
/**
 * A class to reconstruct MFM frames from chunks of data.
 */
class FrameBuilder
{
public:
	FrameBuilder();
	virtual ~FrameBuilder();
	void addDataChunk(const mfm::Byte* begin, const mfm::Byte* end);
	virtual void reset();
protected:
	virtual void processFrame(mfm::Frame & frame) = 0;
	virtual void processHeader(const mfm::PrimaryHeader &) {};
private:
	void buildFrames(const mfm::Byte* end);
	mfm::Serializer serializer_; ///< Buffer storing chunks of frame
	size_t frameSize_; ///< Size [Bytes] of the frame being currently built.
};
//______________________________________________________________________
} /* namespace mfm */
#endif /* mfm_FrameBuilder_h_INCLUDED */
