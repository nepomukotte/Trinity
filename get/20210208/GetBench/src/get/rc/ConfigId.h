/**
 * @file ConfigId.h
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

#ifndef get_rc_ConfigId_h_INCLUDED
#define get_rc_ConfigId_h_INCLUDED

#include "get/rc/SubConfigIdSet.h"
#include <utl/Exception.h>
#include <string>
#include <ostream>
#include <map>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
/**
 * A class to store a triplet of configuration identifiers.
 */
class ConfigId
{
public:
	class UnknownSubType;
	class InvalidConfigId;
	typedef enum { Describe, Prepare, Configure, HardwareDescr, MaxConfigTypes} SubType;
	typedef std::map< std::string, SubType > StringSubTypeMap;
	static const  StringSubTypeMap subTypeMap;
	static SubType stringToSubType(const std::string & name);
	static std::string subTypeToString(const SubType & type);
	typedef std::multimap< ConfigId::SubType, std::string > SubConfigIdMap;
public:
	ConfigId(const std::string & describeId="", const std::string & prepareId="", const std::string & configureId="");
	void insert(const SubType & subType, const std::string & subId);
	SubConfigIdMap & subConfigIds() { return subConfigIds_; }
	const SubConfigIdMap & subConfigIds() const { return subConfigIds_; }
	size_t erase(const SubType & subType);
	void clear();
	bool isComplete() const;
	std::string subConfigId(SubType type) const;
	bool fromXml(const std::string & xmlString);
	std::string toXml() const;
	friend std::ostream & operator<< (std::ostream & stream, const ConfigId & triplet);
	bool operator<(const ConfigId & rhs) const;
private:
	SubConfigIdMap subConfigIds_;
};
//_________________________________________________________________________________________________
/** An exception thrown when encountering an unknown configuration type. */
class ConfigId::UnknownSubType : public ::utl::Exception
{
public:
	UnknownSubType(const std::string & subType) : subType(subType) {}
	virtual ~UnknownSubType() throw() {}
	void writeMessage(std::ostream& stream) const { stream << "Configuration sub-type '" << subType << "' does not exist."; }
	const std::string subType;
};
//_________________________________________________________________________________________________
/** An exception thrown when encountering an unknown configuration type. */
class ConfigId::InvalidConfigId : public ::utl::Exception
{
public:
	InvalidConfigId(const ConfigId & configId) : configId(configId) {}
	virtual ~InvalidConfigId() throw() {}
	void writeMessage(std::ostream& stream) const { stream << "Configuration identifier '" << configId << "' is not valid."; }
	const ConfigId configId;
};
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_ConfigId_h_INCLUDED */
