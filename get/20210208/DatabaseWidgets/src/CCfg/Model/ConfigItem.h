#ifndef CCFG_MODEL_CONFIGITEM_H_INCLUDED
#define CCFG_MODEL_CONFIGITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/Db/Config.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace CCfg
{
namespace Model
{
// forward declarations
class SchemaItem;
class VersionItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a configuration in the database tree. It can enumerate
 * its versions.
 *
 * @author fchateau
 */
class ConfigItem : public AbstractDatabaseItem
{
	SchemaItem* parent;
	Io::Db::ConfigPtr config;
	boost::ptr_vector<VersionItem> versions;
	bool upToDate;

public:
	ConfigItem(SchemaItem* parent, const Io::Db::ConfigPtr& cfg);
	~ConfigItem();

	// @{ Interface implementation
	bool canSetName() const;
	void setName(const std::string& newName);
	const std::string& getName() const;
	size_t getItemCount();
	AbstractDatabaseItem* getItem(size_t index);
	AbstractDatabaseItem* getParent();
	void refresh();
	void removeItem(size_t index, RemoveOptions options = FromList);
	void clear();
	std::string buildUrl() const;
	ItemType getType() const;
	// @}

	SchemaItem* getSchema() const { return parent; }
	Io::Db::Config& getEntity() { return *config; }

protected:
	void ensureUpToDate();

private:
	ConfigItem(const ConfigItem& rvalue);
	ConfigItem& operator=(const ConfigItem& rvalue);
};

}
}

#endif // CCFG_MODEL_CONFIGITEM_H_INCLUDED
