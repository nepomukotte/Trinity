/**
 * @file DataRouterId.cpp
 * @date 15 févr. 2013
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetBench software project.
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

#include "DataRouterId.h"

#include <tinyxml.h>
#include <utl/Logging.h>
#include <boost/lexical_cast.hpp>
#include <sstream>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
DataRouterId::DataRouterId(const std::string & name, const ::utl::net::SocketAddress & sockAddr, const std::string & flowType)
	: name(name), sockAddr(sockAddr), flowType(flowType)
{
}
//_________________________________________________________________________________________________
void DataRouterId::clear()
{
	name.clear();
	sockAddr = ::utl::net::SocketAddress();
}
//_________________________________________________________________________________________________
/**
 * Prints router to a stream as an XML DataRouter element.
 */
std::ostream & operator<< (std::ostream & stream, const ::get::rc::DataRouterId & id)
{
	stream << "<DataRouter "
			<< " name=\"" << id.name << "\""
			<< " ipAddress=\"" << id.sockAddr.ipAddress() << "\""
			<< " port=\"" << id.sockAddr.port() << "\""
			<< " type=\"" << id.flowType << "\" />";
	return stream;
}
//_________________________________________________________________________________________________
/**
 * Returns an XML string representation of the router.
 * @return A DataRouter XML element.
 */
std::string DataRouterId::toXml() const
{
	std::ostringstream oss;
	oss << *this;
	return oss.str();
}
//_________________________________________________________________________________________________
/**
 * Searches for a router in an XML string an parses it.
 * @param xmlString XML string to parse.
 * @return Whether a valid router has been found and parsed.
 */
bool DataRouterId::fromXml(const std::string & xmlString)
{
	clear();

	TiXmlDocument doc;
	doc.Parse(xmlString.c_str());
	if(doc.Error()) return false;

	TiXmlNode* node = doc.FirstChild("DataRouter");
	if (NULL == node) return false;
	TiXmlElement* router = node->ToElement();

	const char* attrValue = router->Attribute("name");
	if (attrValue) name = attrValue;

	attrValue = router->Attribute("ipAddress");
	if (attrValue) sockAddr.ipAddress().fromString(attrValue);

	attrValue = router->Attribute("port");
	if (attrValue)
	{
		try
		{
			sockAddr.port() = boost::lexical_cast< uint16_t >(attrValue);
		}
		catch (const boost::bad_lexical_cast & e)
		{
			LOG_ERROR() << '\'' << attrValue << "' is not a valid port number!";
			return false;
		}
	}

	attrValue = router->Attribute("type");
	if (attrValue) flowType = attrValue;

	return true;
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */



