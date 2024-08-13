/**
 * @file ConfigManager.h
 * @date 9 oct. 2012
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

#ifndef get_rc_ConfigManager_h_INCLUDED
#define get_rc_ConfigManager_h_INCLUDED

#include "get/rc/ConfigId.h"
#include <utl/Exception.h>
#include <CCfg/Io/Url.h>
#include <CCfg/Object.h>
#include <memory>
#include <string>
#include <vector>

namespace get
{
namespace rc
{
//_________________________________________________________________________________________________
// Forward declarations
class ConfigIdSet;
class SubConfigIdSet;
//_________________________________________________________________________________________________
/**
 * A class to manage the configurations available to the Electronic Control Core.
 */
class ConfigManager
{
public:
	class ConfigIdNotFound;
	class ConfigNotFound;
public:
	ConfigManager();
	virtual ~ConfigManager() {};
	/// @name Getters
	/// @{
	virtual void getSubConfigIds(const ConfigId::SubType type, SubConfigIdSet & subConfigIds) const = 0;
	void getConfigIds(ConfigIdSet &) const;
	virtual std::auto_ptr< CCfg::Object > getConfig(const CCfg::Io::Url & cfgUrl) const = 0;
	virtual void getSubConfig(const ConfigId::SubType subType, const std::string & subConfigId, std::string & cfg) const;
	std::auto_ptr< CCfg::Object > getDescribeConfig(const std::string & subConfigId) const;
	std::auto_ptr< CCfg::Object > getPrepareConfig(const std::string & subConfigId) const;
	std::auto_ptr< CCfg::Object > getConfigureConfig(const std::string & subConfigId) const;
	std::auto_ptr< CCfg::Object > getHardwareDescription(const std::string & subConfigId) const;
	virtual ConfigId getFullConfigId(const ConfigId & id) const;
	virtual std::string getConfigSet(const ConfigId & configId) const;
	/// @}
	/// @name Setters
	/// @{
	std::string addSubConfig(const ConfigId::SubType subType, const std::string & subConfig, const std::string & subConfigId = "") const;
	std::string addDescribeConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId = "") const;
	std::string addPrepareConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId = "") const;
	std::string addConfigureConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId = "") const;
	std::string addHardwareConfig(std::auto_ptr< CCfg::Object > cfg, const std::string & subConfigId = "") const;
	/// @}
protected:
	virtual std::auto_ptr< CCfg::Object > getSubConfig(const ConfigId::SubType subType, const std::string & subConfigId) const = 0;
	virtual std::string addSubConfig(const ConfigId::SubType  subType, std::auto_ptr< CCfg::Object >, const std::string & subConfigId = "") const = 0;
};
//_________________________________________________________________________________________________
/** An exception thrown when trying to access to a configuration with an unknown ID. */
class ConfigManager::ConfigIdNotFound : public ::utl::Exception
{
public:
	ConfigIdNotFound(const std::string & configID) : configID(configID) {}
	virtual ~ConfigIdNotFound() throw() {}
	void writeMessage(std::ostream& stream) const { stream << "Configuration '" << configID << "' does not exist."; }
	const std::string configID;
};
//_________________________________________________________________________________________________
/** An exception thrown when failing to access a configuration URL. */
class ConfigManager::ConfigNotFound : public ::utl::Exception
{
public:
	ConfigNotFound(const std::string & configUrl) : configUrl(configUrl) {}
	virtual ~ConfigNotFound() throw() {}
	void writeMessage(std::ostream& stream) const { stream << "Configuration '" << configUrl << "' does not exist."; }
	const std::string configUrl;
};
//_________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_ConfigManager_h_INCLUDED */
