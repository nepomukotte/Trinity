#ifndef CCFG_MODEL_SCHEMAITEM_H_INCLUDED
#define CCFG_MODEL_SCHEMAITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/Db/ConfigDatabase.h"
#ifndef Q_MOC_RUN   // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/ptr_container/ptr_vector.hpp>
#endif

namespace CCfg
{
namespace Model
{
// forward declarations
class ConnectionItem;
class ConfigItem;
class SiteItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a database schema in the tree. It can enumerate the
 * configurations it contains.
 *
 * @author fchateau
 */
class SchemaItem : public AbstractDatabaseItem
{
	ConnectionItem* parent;
	Io::Db::ConfigDatabasePtr database;
	boost::ptr_vector<ConfigItem> configurations;
	boost::ptr_vector<SiteItem> sites;
	bool upToDate;

public:
	enum DisplayMode {Config = 0x01, Run = 0x02, All = 0x03};

	SchemaItem(ConnectionItem* parent, const Io::Db::ConfigDatabasePtr& db, DisplayMode mode = All);
	~SchemaItem();

	void setDisplayMode(DisplayMode newMode);

	void changeUserPassword(const std::string& newPassword);
	void createConfig(const std::string& name);
	void createSite(const std::string& name);

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

	ConnectionItem* getConnection() const { return parent; }
	Io::Db::ConfigDao& getDao();

protected:
	void ensureUpToDate();

private:
	SchemaItem(const SchemaItem& rvalue);
	SchemaItem& operator=(const SchemaItem& rvalue);
	DisplayMode displayMode;
};

}
}

#endif // CCFG_MODEL_SCHEMAITEM_H_INCLUDED
