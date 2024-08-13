/**
 * @file SubConfigIdSet.cpp
 * @date 15 févr. 2013
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench-GANIL software project.
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

#include "SubConfigIdSet.h"

#include <tinyxml.h>
#include <cstddef>
#include <sstream>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
/**
 * Constructs an empty set.
 */
SubConfigIdSet::SubConfigIdSet()
{
	;
}
//_________________________________________________________________________________________________
/**
 * Constructs set of configuration identifiers from a set of strings.
 */
SubConfigIdSet::SubConfigIdSet(const std::set< std::string > & strings)
	: std::set< std::string >(strings)
{
	;
}
//_________________________________________________________________________________________________
/**
 * Prints the set of identifiers to a stream as an XML ConfigIdSet element.
 */
std::ostream & operator<< (std::ostream & stream, const ::get::rc::SubConfigIdSet & set)
{
	stream << "<SubConfigIdSet>\n";
	SubConfigIdSet::const_iterator idIter;
	for (idIter = set.begin(); idIter != set.end();  idIter++)
	{
		stream << "\t<SubConfigId>" <<  *idIter << "</SubConfigId>\n";
	}
	stream << "</SubConfigIdSet>";
	return stream;
}
//_________________________________________________________________________________________________
/**
 * Returns an XML string representation of the set of identifiers.
 * @return A ConfigIdSet XML element.
 */
std::string SubConfigIdSet::toXml() const
{
	std::ostringstream oss;
	oss << *this;
	return oss.str();
}
//_________________________________________________________________________________________________
/**
 * Searches for a set of identifiers in an XML string an parses it.
 * @param xmlString XML string to parse.
 * @return Whether a valid set has been found and parsed.
 */
bool SubConfigIdSet::fromXml(const std::string & xmlString)
{
	clear();

	TiXmlDocument doc;
	doc.Parse(xmlString.c_str());
	if(doc.Error()) return false;

	TiXmlNode* setNode = doc.FirstChild("SubConfigIdSet");
	if (NULL == setNode) return false;
	TiXmlNode* idNode = setNode->FirstChild("SubConfigId");
	while (idNode != NULL)
	{
		insert(idNode->ToElement()->GetText());
		idNode = setNode->IterateChildren("SubConfigId", idNode);
	}
	return true;
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
