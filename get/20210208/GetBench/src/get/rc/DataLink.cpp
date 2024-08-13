/**
 * @file DataLink.cpp
 * @date 11 mars 2013
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

#include "DataLink.h"

#include <tinyxml.h>
#include <sstream>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
DataLink::DataLink()
{
}
//_________________________________________________________________________________________________
DataLink::DataLink(const DataSenderId & sender, const DataRouterId & router)
	: sender(sender), router(router)
{
}
//_________________________________________________________________________________________________
/**
 * Prints data link to a stream as an XML DataLink element.
 */
std::ostream & operator<< (std::ostream & stream, const ::get::rc::DataLink & link)
{
	stream << "<DataLink>" << link.sender << link.router << "</DataLink>";
	return stream;
}
//_________________________________________________________________________________________________
/**
 * Returns an XML string representation of the data link.
 * @return A DataLink XML element.
 */
std::string DataLink::toXml() const
{
	std::ostringstream oss;
	oss << *this;
	return oss.str();
}
//_________________________________________________________________________________________________
/**
 * Searches for a data link in an XML string an parses it.
 * @param xmlString XML string to parse.
 * @return Whether a valid data link has been found and parsed.
 */
bool DataLink::fromXml(const std::string & xmlString)
{
	clear();

	TiXmlDocument doc;
	doc.Parse(xmlString.c_str());
	if(doc.Error()) return false;

	TiXmlNode* link = doc.FirstChild("DataLink");
	if (NULL == link) return false;

	TiXmlNode* node = link->FirstChild("DataSender");
	if (NULL == node) return false;
	TiXmlElement* senderElement = (TiXmlElement*) node->ToElement();
	if (senderElement)
	{
		std::ostringstream oss;
		oss << *senderElement;
		if (not sender.fromXml(oss.str())) return false;
	}

	node = link->FirstChild("DataRouter");
	if (NULL == node) return false;
	TiXmlElement* routerElement = (TiXmlElement*) node->ToElement();
	if (routerElement)
	{
		std::ostringstream oss;
		oss << *routerElement;
		if (not router.fromXml(oss.str())) return false;
	}

	return true;
}
//_________________________________________________________________________________________________
/**
 * Compares two data links.
 */
bool DataLink::operator<(const DataLink & rhs) const
{
	return sender < rhs.sender;
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
