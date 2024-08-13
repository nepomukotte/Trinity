/*
 * $Id: FrameDictionary.cpp 716 2012-04-10 15:17:58Z psizun $
 * @file FrameDictionary.cpp
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

#include "FrameDictionary.h"
#include <CCfg/CConfig.h>
#include <CCfg/Io/Document.h>
#include <CCfg/View/ClassIterator.h>
using namespace CCfg;
#include <sstream>
#include <iostream>

namespace mfm {
//______________________________________________________________________
FrameDictionary::FrameDictionary()
{
	;
}
//______________________________________________________________________
FrameDictionary & FrameDictionary::instance()
{
	static FrameDictionary dictionary;
	return dictionary;
}
//______________________________________________________________________
size_t FrameDictionary::addFormats(const std::string & path)
{
	size_t numFormats = 0;
	Io::Document doc;
	CConfig config(&doc.loadFromFile(path));

	// Find all concrete objects with class path MFM[*].Frame[*].Revision[*]
	View::ClassIterator formatIter = config("Frame", "*")("Revision").getClassInstances();
	while (formatIter.hasNext())
	{
		// Do not forget to call next before accessing the first element.
		formatIter.next();
		// Test whether object is concrete
		if (not (formatIter->toObject().isDefaultConfig()
				and formatIter->getParent()->isDefaultConfig()))
		{
			std::string formatName = formatIter->getParent()->getIndex();
			std::istringstream iss(formatIter->getIndex());
			int formatRevision;
			iss >> formatRevision;
			addFormat(FrameFormat(path, formatName, formatRevision));
			numFormats++;
		}
	}
	return numFormats;
}
//______________________________________________________________________
bool FrameDictionary::addFormat(const FrameFormat & format)
{
	// Find 'frameType' and 'revision'
	FormatKey key(format.frameType(), format.revision());
	if (formats_.count(key) > 0)
	{
		return false;
		// FIXME: throw DuplicateFormat
	}
	formats_[key] = format;
	return true;
}
//______________________________________________________________________
/**
 * Prints the list of all known frame formats to STDOUT.
 * @return Number of known formats.
 */
size_t FrameDictionary::listFormats() const
{
	FormatMap::const_iterator formatIter = formats_.begin();
	size_t count = 0;
	while (formatIter != formats_.end())
	{
		std::cout << "Format #" << count << " : name='" << formatIter->second.formatName()
				<< "'\ttype=" << formatIter->second.frameType()
				<< "\trevision=" << formatIter->second.revision() << std::endl;
		++formatIter;
		++count;
	}
	return count;
}
//______________________________________________________________________
/**
 * Finds and returns frame format with given type and revision.
 * @param type Frame type numeric identifier.
 * @param revision Revision of frame type.
 * @return Returns the format found.
 * @throws FormatRevisionNotFound If the dictionary contains no format with such type and revision.
 */
FrameFormat const & FrameDictionary::findFormat(FrameType const & type, FormatRevision const & revision) const
{
	FormatKey searchKey(type, revision);
	if (0 == formats_.count(searchKey))
	{
		throw mfm::FormatRevisionNotFound(type, revision);
	}
	return formats_.find(searchKey)->second;
}
//______________________________________________________________________
/**
 * Finds the latest revision format for the frame with the given type.
 * @param type Frame type numeric identifier.
 * @return Returns the format found.
 * @throws FormatNotFound If the dictionary contains no format with such type.
 */
FrameFormat const & FrameDictionary::findLatestFormat(FrameType const & type) const
{
	bool formatFound = false;
	FormatRevision maxRevision = 0;
	FormatMap::const_iterator frameIter = formats_.begin();
	while (frameIter != formats_.end())
	{
		FormatKey formatKey = frameIter->first;
		if (type == formatKey.first)
		{
			formatFound = true;
			if (formatKey.second > maxRevision)
				maxRevision = formatKey.second;
		}
		frameIter++;
	}
	if (not formatFound)
	{
		throw mfm::FormatNotFound(type);
	}
	return findFormat(type, maxRevision);
}
//______________________________________________________________________
} /* namespace mfm */
