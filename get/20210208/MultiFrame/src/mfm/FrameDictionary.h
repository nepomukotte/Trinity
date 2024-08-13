/*
 * $Id: FrameDictionary.h 716 2012-04-10 15:17:58Z psizun $
 * @file FrameDictionary.h
 * @created 21 févr. 2012
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

#ifndef mfm_FrameDictionary_h_INCLUDED
#define mfm_FrameDictionary_h_INCLUDED

#include "mfm/FrameFormat.h"
#include <string>
#include <map>
#include <utility> // pair
#include <boost/cstdint.hpp>

namespace mfm {
//______________________________________________________________________
/**
 * Singleton class managing the MFM format descriptions.
 */
class FrameDictionary
{
private:
	/// Constructor is private (this is a singleton class).
	FrameDictionary();
	virtual ~FrameDictionary() {}
	/// Copy constructor (intentionally undefined).
	FrameDictionary(const FrameDictionary &);
	/// Assignment operator (intentionally undefined).
	FrameDictionary & operator=(const FrameDictionary &);
public:
	/// Accessor to single dictionary instance.
	static FrameDictionary & instance();
	/// Add given frame format description to dictionary.
	bool addFormat(const FrameFormat & format);
	/**
	 * Loads all frame format descriptions contained in file with given path.
	 */
	size_t addFormats(const std::string & path);
	size_t listFormats() const;
	typedef uint16_t FrameType; ///< Integer value used to encode the type of a frame within its header.
	typedef uint8_t FormatRevision; ///< Revision of a frame type.
	FrameFormat const & findFormat(FrameType const & , FormatRevision const &) const;
	FrameFormat const & findLatestFormat(FrameType const & type) const;
private:
	/// Dictionary key (frameType, revision)
	typedef std::pair< FrameType, FormatRevision > FormatKey; ///< Format revision unique identifier.
	typedef std::map< FormatKey, FrameFormat > FormatMap; ///< Container used to map frame formats to their (type, revision) identifier.
	FormatMap formats_; ///< Map of formats.
};
//______________________________________________________________________
} /* namespace mfm */
#endif /* mfm_FrameDictionary_h_INCLUDED */
