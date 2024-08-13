/**
 * @file ConfigId.cpp
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

#include "ConfigId.h"

#include <tinyxml.h>

#include <boost/tuple/tuple_comparison.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/case_conv.hpp>
namespace ba = boost::algorithm;
#include <boost/assign/list_of.hpp> // for 'map_list_of()'

template < typename T2 >
class SecondComparator
{
public:
	SecondComparator(T2 const & refValue) : refValue_(refValue) {}
	template< typename T1 > bool operator() ( std::pair< T1, T2 > const & p) const
	{
		return p.second == refValue_;
	}
private:
	T2 const & refValue_;
};

template<class T2>
SecondComparator< T2 > is_second_equal(T2 const & refValue)
{
	return SecondComparator< T2 >(refValue);
}

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
const ConfigId::StringSubTypeMap ConfigId::subTypeMap
	 = boost::assign::map_list_of
	 ("describe", ConfigId::Describe)
	 ("prepare", ConfigId::Prepare)
	 ("configure", ConfigId::Configure)
	 ("hardware", ConfigId::HardwareDescr);
//_________________________________________________________________________________________________
ConfigId::SubType ConfigId::stringToSubType(const std::string & name)
{
	std::string name_(name);
	ba::to_lower(name_);

	StringSubTypeMap::const_iterator it = subTypeMap.find(name);
	if (it == subTypeMap.end())
		throw UnknownSubType(name);

	return it->second;
}
//_________________________________________________________________________________________________
std::string ConfigId::subTypeToString(const ConfigId::SubType & type)
{
	 StringSubTypeMap::const_iterator it = std::find_if(subTypeMap.begin(), subTypeMap.end(),
			 is_second_equal(type));

	return it->first;
}
//_________________________________________________________________________________________________
/**
 * Constructs a triplet of configuration identifiers.
 * @param describeId Identifier if configuration used for Describe
 * @param prepareId Identifier if configuration used for Prepare
 * @param configureId Identifier if configuration used for Configure
 */
ConfigId::ConfigId(const std::string & describeId, const std::string & prepareId, const std::string & configureId)
{
	subConfigIds_.insert(SubConfigIdMap::value_type(Describe, describeId));
	subConfigIds_.insert(SubConfigIdMap::value_type(Prepare, prepareId));
	subConfigIds_.insert(SubConfigIdMap::value_type(Configure, configureId));
}
//_________________________________________________________________________________________________
void ConfigId::insert(const SubType & subType, const std::string & subId)
{
	subConfigIds_.insert(SubConfigIdMap::value_type(subType, subId));
}
//_________________________________________________________________________________________________
size_t ConfigId::erase(const SubType & subType)
{
	return subConfigIds_.erase(subType);
}
//_________________________________________________________________________________________________
/**
 * Clears the three identifiers.
 */
void ConfigId::clear()
{
	subConfigIds_.clear();
}
//_________________________________________________________________________________________________
/**
 * Clears the three identifiers.
 */
bool ConfigId::isComplete() const
{
	return subConfigIds_.count(Describe) and
			subConfigIds_.count(Prepare) and
			subConfigIds_.count(Configure);
}
//_________________________________________________________________________________________________
std::string ConfigId::subConfigId(ConfigId::SubType subType) const
{
	SubConfigIdMap::const_iterator iter = subConfigIds_.find(subType);
	if (iter == subConfigIds_.end())
	{
		throw InvalidConfigId(*this);
	}
	return iter->second;
}
//_________________________________________________________________________________________________
/**
 * Prints configuration ID to a stream as an XML ConfigId element.
 */
std::ostream & operator<< (std::ostream & stream, const ::get::rc::ConfigId & id)
{
	stream << "<ConfigId>";
	ConfigId::SubConfigIdMap::const_iterator subIdIter;
	for (subIdIter = id.subConfigIds_.begin(); subIdIter != id.subConfigIds_.end(); ++subIdIter)
	{
		std::string typeString = ConfigId::subTypeToString(subIdIter->first);
		stream << "<SubConfigId type=\"" << typeString << "\">" << subIdIter->second << "</SubConfigId>";
	}
	stream << "</ConfigId>";
	return stream;
}
//_________________________________________________________________________________________________
/**
 * Returns an XML string representation of the configuration ID.
 * @return A ConfigId XML element.
 */
std::string ConfigId::toXml() const
{
	std::ostringstream oss;
	oss << *this;
	return oss.str();
}
//_________________________________________________________________________________________________
/**
 * Searches for a configuration ID in an XML string an parses it.
 * @param xmlString XML string to parse.
 * @return Whether a valid triplet of sub-configuration IDs has been found and parsed.
 */
bool ConfigId::fromXml(const std::string & xmlString)
{
	clear();

	TiXmlDocument doc;
	doc.Parse(xmlString.c_str());
	if(doc.Error()) return false;

	TiXmlNode* configIdNode = doc.FirstChild("ConfigId");
	if (NULL == configIdNode) return false;

	TiXmlNode* subConfigIdNode = configIdNode->FirstChild("SubConfigId");
	while (subConfigIdNode != NULL)
	{
		TiXmlElement* element = subConfigIdNode->ToElement();
		ConfigId::SubType subConfigType = ConfigId::stringToSubType(element->Attribute("type"));
		subConfigIds_.insert(SubConfigIdMap::value_type(subConfigType, element->GetText()));
		subConfigIdNode = configIdNode->IterateChildren("SubConfigId", subConfigIdNode);
	}
	return isComplete();
}
//_________________________________________________________________________________________________
/**
 * Compares two configuration IDs.
 */
bool ConfigId::operator<(const ConfigId & rhs) const
{
	return subConfigIds_ < rhs.subConfigIds_;
}
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
