/**
 * @file GFrameHeader.h
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

#ifndef get_GFrameHeader_h_INCLUDED
#define get_GFrameHeader_h_INCLUDED

#include <TRint.h>
#include <Rtypes.h>
#include <TObject.h>
#include <TBits.h>
#include <RVersion.h> // ROOT_VERSION_CODE

#if ROOT_VERSION_CODE > 393216
namespace GET {
#else
namespace get {
#endif
//__________________________________________________________________________________________________
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GFrameHeader                                                         //
//                                                                      //
// Header of a GET physics data frame.                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class GFrameHeader :  public TObject
{
public:
	GFrameHeader();
	void Clear(Option_t *option = "");

	UShort_t fDataSource; // Data source ID (1 Byte).
	UShort_t fRevision; // GET format revision (1 Byte).
	ULong_t fEventTime; // Time offset of the frame as compared to start of run. Expressed in units of experiment clock. (6 Bytes)
	UInt_t fEventIdx; // The event to which this frame belongs within a run. (4 Bytes)
	UShort_t fCoboIdx; // Index of CoBo subsystem from which the frame originates. (1 Byte)
	UShort_t fAsadIdx; // Index of AsAd board from which the frame data originates, relative to the CoBo subsystem to which it is connected. (1 Byte)
	UShort_t fReadOffset; // AGET readout offset. (2 Bytes)
	UShort_t fStatus; // Frame status info. (1 Byte)
	TBits  fHitPatterns[4]; // Bit patterns indicating which of the 68 channels of AGET's have been hit. (4 x 9 Bytes)
	UShort_t fMult[4]; // AGET multiplicities (4 x 2 Bytes).
	UInt_t fWindowOut; // Output of the sliding window at the time of the event trigger. (4 Bytes)
	UShort_t fLastCellIdx[4]; // Last cell index read by AGET's. (4 x 2 Bytes)
	static const UShort_t MAX_CHANNELS = 68;
private:
	ClassDef(GFrameHeader, 1);
};
//__________________________________________________________________________________________________
} /* namespace get */
#endif /* get_GFrameHeader_h_INCLUDED */
