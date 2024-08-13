/**
 * @file DiscardStorage.h
 * @date 24 sept. 2012
 * @author sizun
 * SVN Id: $Id$
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Patrick Sizun (patrick.sizun@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This file is part of the GetController software project.
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

#ifndef get_daq_DiscardStorage_h_INCLUDED
#define get_daq_DiscardStorage_h_INCLUDED

#include <get/daq/FrameStorage.h>

namespace get {
namespace daq {
//_________________________________________________________________________________________________
/**
 * A class to simply count he number of bytes and frames received and discard them
 */
class DiscardStorage : public get::daq::FrameStorage
{
public:
	DiscardStorage();
	virtual void processFrame(mfm::Frame & frame); ///< @override from get::daq::::FrameStorage
private:
	DiscardStorage(const DiscardStorage& r);
	DiscardStorage& operator=(const DiscardStorage& r);
private:
	size_t numBytes, numFrames;
};
//_________________________________________________________________________________________________
} /* namespace daq */
} /* namespace get */
#endif /* get_daq_DiscardStorage_h_INCLUDED */
