/**
 * @file GDataSample.h
 * @date 25 févr. 2013
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

#ifndef get_GDataSample_h_INCLUDED
#define get_GDataSample_h_INCLUDED

#include <TRint.h>
#include <Rtypes.h>
#include <TObject.h>
#include <RVersion.h> // ROOT_VERSION_CODE

#if ROOT_VERSION_CODE > 393216
namespace GET {
#else
namespace get {
#endif
//__________________________________________________________________________________________________
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// GDataSample                                                          //
//                                                                      //
// An item within a GET physics data frame.                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
class GDataSample :  public TObject
{
public:
	GDataSample(const UShort_t & buckIdx=0, const UShort_t & value=0);
	virtual ~GDataSample() { ; }
	void Set(const UShort_t & buckIdx=0, const UShort_t & value=0);

	UShort_t fBuckIdx; // Time bucket index attached to the digitized data sample. Range: 0 to 511.
	UShort_t fValue; // Digitized value of the signal on the corresponding channel (12 bits).

	ClassDef(GDataSample, 2);
};
//__________________________________________________________________________________________________
} /* namespace get */
#endif /* get_GDataSample_h_INCLUDED */
