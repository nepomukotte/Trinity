/**
 * @file ConfigIdSet.cpp
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

#include "ConfigIdSet.h"
#include <tinyxml.h>
#include <sstream>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
/**
 * Constructs an empty set of configuration identifiers.
 */
ConfigIdSet::ConfigIdSet()
{
	;
}
//_________________________________________________________________________________________________
/**
 * Inserts a triplet into the set.
 */
bool ConfigIdSet::insert(const ConfigId & triplet)
{
	return std::set< ConfigId >::insert(triplet).second;
}
//_________________________________________________________________________________________________
/**
 * Creates a triplet from three configuration IDs and inserts it into the set.
 */
bool ConfigIdSet::insert(const std::string & describeId, const std::string & prepareId, const std::string & configureId)
{
	ConfigId triplet(describeId, prepareId, configureId);
	return std::set< ConfigId >::insert(triplet).second;
}
//_________________________________________________________________________________________________
/**
 * Prints the set of configuration IDs to a stream as an XML ConfigIdSet element.
 */
std::ostream & operator<< (std::ostream & stream, const ::get::rc::ConfigIdSet & set)
{
	stream << "<ConfigIdSet>\n";
	ConfigIdSet::const_iterator idIter;
	for (idIter = set.begin(); idIter != set.end();  idIter++)
	{
		stream << '\t' << *idIter << '\n';
	}
	stream << "</ConfigIdSet>";
	return stream;
}
//_________________________________________________________________________________________________
/**
 * Returns an XML string representation of the set of configuration identifiers.
 * @return A ConfigIdSet XML element.
 */
std::string ConfigIdSet::toXml() const
{
	std::ostringstream oss;
	oss << *this;
	return oss.str();
}
//_________________________________________________________________________________________________
/**
 * Searches for a set of configuration identifiers in an XML string an parses it.
 * @param xmlString XML string to parse.
 * @return Whether a valid configuration identifier set has been found and parsed.
 */
bool ConfigIdSet::fromXml(const std::string & xmlString)
{
	clear();

	TiXmlDocument doc;
	doc.Parse(xmlString.c_str());
	if(doc.Error()) return false;

	TiXmlNode* setNode = doc.FirstChild("ConfigIdSet");
	if (NULL == setNode) return false;
	TiXmlNode* tripletNode = setNode->FirstChild("ConfigId");
	while (tripletNode != NULL)
	{
		std::ostringstream oss;
		oss << *tripletNode;
		ConfigId triplet;
		if (triplet.fromXml(oss.str()))
			this->insert(triplet);
		tripletNode = setNode->IterateChildren("ConfigId", tripletNode);
	}
	return true;
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
