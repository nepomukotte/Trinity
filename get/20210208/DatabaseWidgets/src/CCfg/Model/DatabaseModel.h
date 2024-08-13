#ifndef CCFG_MODEL_DATABASEMODEL_H_INCLUDED
#define CCFG_MODEL_DATABASEMODEL_H_INCLUDED

#include "ConnectionRootItem.h"
#include "CCfg/Io/DatabaseStorage.h"
#include <QAbstractItemModel>
#include <QIcon>

namespace CCfg
{
namespace Model
{
// forward declarations
class SiteItem;
class SchemaItem;
class RunItem;
class VersionItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is a tree model representing the content of a list of configuration
 * databases, as a hierarchy: from connections to configuration versions.
 *
 * @author fchateau
 */
class DatabaseModel : public QAbstractItemModel
{
	ConnectionRootItem root;
	QIcon connectedIcon, disconnectedIcon;
	const static int TypeRole = Qt::UserRole;
public:
	DatabaseModel(QObject* parent = 0);
	virtual ~DatabaseModel();

	void setLoginCallback(const Io::LoginCallbackPtr& callback);

	// @{ Item creation
	void createConnection(const std::string& name, Io::DbmsType dbmsType,
			const std::string& host, const std::string& port,
			const std::string& user,
			const std::string& sid, const std::string& schema);
	void createSchema(ConnectionItem* parent, const std::string& name, const std::string& password);
	void createConfig(SchemaItem* schema, const std::string& name);
	void createSite(SchemaItem* schema, const std::string& name);
	void createRun(SiteItem* site);
	void createCondition(RunItem* run, VersionItem* version);
	// @}

	// @{
	void clear(const QModelIndex& parent);
	void refresh(const QModelIndex& parent);
	// @}

	// @{ Overriden methods
	QModelIndex index(int row, int column, const QModelIndex& parent) const;
	QModelIndex parent(const QModelIndex& index) const;

	bool hasChildren(const QModelIndex& parent) const;
	int rowCount(const QModelIndex& parent) const;
	int columnCount(const QModelIndex& parent) const;

	Qt::ItemFlags flags(const QModelIndex& index) const;
	QVariant headerData(int section, Qt::Orientation orientation, int role) const;
	QVariant data(const QModelIndex& index, int role) const;

	bool setData(const QModelIndex& index, const QVariant& value, int role);
	bool removeRows(int row, int count, const QModelIndex& parent);
	// @}

	// @{ index computation
	static AbstractDatabaseItem* indexToItem(const QModelIndex& index);
	QModelIndex itemToIndex(AbstractDatabaseItem* item) const;
	int getRowInParent(AbstractDatabaseItem* item) const;
	// @}

private:
	DatabaseModel(const DatabaseModel& rvalue);
	DatabaseModel& operator=(const DatabaseModel& rvalue);
};

}
}

#endif // CCFG_MODEL_DATABASEMODEL_H_INCLUDED
