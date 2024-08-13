#ifndef CCFG_MODEL_VERSIONITEM_H_INCLUDED
#define CCFG_MODEL_VERSIONITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/Db/ConfigVersion.h"

namespace CCfg
{
namespace Model
{
// forward declarations
class ConfigItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a configuration version in the database tree.
 *
 * @author fchateau
 */
class VersionItem : public AbstractDatabaseItem
{
	ConfigItem* parent;
	Io::Db::ConfigVersionPtr version;
	std::string name;
	bool structureColor;///< each structure change in the version array toggles the color.

public:
	VersionItem(ConfigItem* parent, const Io::Db::ConfigVersionPtr& version, bool structColor);
	~VersionItem();

	// @{ Interface implementation
	bool isLeaf() const;
	bool canSetName() const;
	void setName(const std::string& newName);
	const std::string& getName() const;
	const std::string& getDisplayName() const;
	size_t getItemCount();
	AbstractDatabaseItem* getItem(size_t index);
	AbstractDatabaseItem* getParent();
	void refresh();
	void removeItem(size_t index, RemoveOptions options = FromList);
	void clear();
	std::string buildUrl() const;
	ItemType getType() const;
	// @}

	ConfigItem* getConfig() const { return parent; }
	Io::Db::ConfigVersion&  getEntity() const { return *version; }
	bool getStructureColor() const { return structureColor; }

protected:
	void updateDisplayName();

private:
	VersionItem(const VersionItem& rvalue);
	VersionItem& operator=(const VersionItem& rvalue);
};

}
}

#endif // CCFG_MODEL_VERSIONITEM_H_INCLUDED
