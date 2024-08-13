#include "DatabaseModel.h"
#include "ConnectionItem.h"
#include "SchemaItem.h"
#include "VersionItem.h"
#include "ConfigItem.h"
#include "SiteItem.h"
#include "RunItem.h"
#include "ConditionItem.h"
#include "Utils/Error/ExcepHandling.h"
#include <QBrush>
#include <cassert>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::DatabaseModel class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
DatabaseModel::DatabaseModel(QObject* parent) : QAbstractItemModel(parent),
	root(), connectedIcon(":/dbbrowser/connected"),
	disconnectedIcon(":/dbbrowser/disconnected")
{
}

/**
 * Destructor.
 */
DatabaseModel::~DatabaseModel()
{
}

/**
 * Changes the callback associated to login requests on connection items.
 * The modification will only impact new connection items.
 */
void DatabaseModel::setLoginCallback(const Io::LoginCallbackPtr& callback)
{
	root.setLoginCallback(callback);
}

/**
 * Creates a connection item and appends it at the end of the connection list.
 */
void DatabaseModel::createConnection(const std::string& name, Io::DbmsType dbmsType,
		const std::string& host, const std::string& port,
		const std::string& user,
		const std::string& sid, const std::string& schema)
{
	std::auto_ptr<Model::ConnectionItem> conn(
		new ConnectionItem(name, dbmsType, host, port, user, std::string(), sid, schema));

	int n = root.getItemCount();
	beginInsertRows(itemToIndex(&root), n, n);
	root.appendConnection(conn);
	endInsertRows();
}

/**
 * Creates a new schema on the server of the specified connection.
 */
void DatabaseModel::createSchema(ConnectionItem* parent, const std::string& name, const std::string& password)
{
	parent->createSchema(name, password);
	refresh(itemToIndex(parent));
}

/**
 * Creates a new config in the specified schema.
 */
void DatabaseModel::createConfig(SchemaItem* schema, const std::string& name)
{
	schema->createConfig(name);
	refresh(itemToIndex(schema));
}

/**
 * Creates a new site in the specified schema.
 */
void DatabaseModel::createSite(SchemaItem* schema, const std::string& name)
{
	schema->createSite(name);
	refresh(itemToIndex(schema));
}

/**
 * Creates a new run in the specified site.
 */
void DatabaseModel::createRun(SiteItem* site)
{
	site->createRun();
	refresh(itemToIndex(site));
}

/**
 * Creates a new condition in the specified run.
 */
void DatabaseModel::createCondition(RunItem* run, VersionItem* version)
{
	try {
		run->createCondition(version);
		refresh(itemToIndex(run));
	}
	CATCH_AND_REPORT
}

/**
 * Removes all children of the specified model item.
 */
void DatabaseModel::clear(const QModelIndex& parent)
{
	AbstractDatabaseItem* item = indexToItem(parent);
	int n = item->getItemCount();
	if(n > 0)
	{
		beginRemoveRows(parent, 0, n-1);
		item->clear();
		endRemoveRows();
	}
}

/**
 * Rebuilds the list of children of the specified model item.
 */
void DatabaseModel::refresh(const QModelIndex& parent)
{
	clear(parent);
	AbstractDatabaseItem* item = indexToItem(parent);
	//beginInsertRows(parent, 0, 0);
	item->refresh();
	int n = item->getItemCount();
	beginInsertRows(parent, 0, n-1);
	endInsertRows();
}

/*
 * (see Qt documentation)
 */
QModelIndex DatabaseModel::index(int row, int column, const QModelIndex& parent) const
{
	try
	{
		if(parent.isValid())
		{
			AbstractDatabaseItem* item = indexToItem(parent);
			if(static_cast<size_t>(row) < item->getItemCount() and row >= 0)
			{
				return createIndex(row, column, item->getItem(row));
			}
		}
		else	// root requested
		{
			if(row == 0)
			{
				return createIndex(row, column, const_cast<ConnectionRootItem*>(&root));
				// only return a valid index for row 0 because the root is unique.
			}
		}
	}
	CATCH_AND_REPORT
	return QModelIndex();
}

/*
 * (see Qt documentation)
 */
QModelIndex DatabaseModel::parent(const QModelIndex& index) const
{
	try
	{
		if(index.isValid())
		{
			AbstractDatabaseItem* item = indexToItem(index);
			AbstractDatabaseItem* parent = item->getParent();
			if(parent != 0)
			{
				int row = getRowInParent(parent);
				if(row != -1)
				{
					return createIndex(row, 0, parent);
				}
			}
		}
	}
	CATCH_AND_REPORT
	return QModelIndex();
}

/*
 * (see Qt documentation)
 */
bool DatabaseModel::hasChildren(const QModelIndex& parent) const
{
	try
	{
		if(parent.isValid())
		{
			AbstractDatabaseItem* item = indexToItem(parent);
			return not item->isLeaf();
		}
		else
		{
			return true; // there is always a root node.
		}
	}
	CATCH_AND_REPORT
	return false;
}

/*
 * (see Qt documentation)
 */
int DatabaseModel::rowCount(const QModelIndex& parent) const
{
	try
	{
		if(parent.isValid())
		{
			AbstractDatabaseItem* item = indexToItem(parent);
			return item->getItemCount();
		}
		else
		{
			return 1; // one root node.
		}
	}
	CATCH_AND_REPORT
	return 0;
}

/*
 * (see Qt documentation)
 */
int DatabaseModel::columnCount(const QModelIndex& /*parent*/) const
{
	return 1;
}

/*
 * (see Qt documentation)
 */
Qt::ItemFlags DatabaseModel::flags(const QModelIndex& index) const
{
	Qt::ItemFlags flags;
	try
	{
		if(index.isValid())
		{
			AbstractDatabaseItem* item = indexToItem(index);
			flags |= Qt::ItemFlags(Qt::ItemIsEnabled) | Qt::ItemIsSelectable;
			if(item->canSetName())
			{
				flags |= Qt::ItemIsEditable;
			}
		}
	}
	CATCH_AND_REPORT
	return flags;
}

/*
 * (see Qt documentation)
 */
QVariant DatabaseModel::headerData(int /*section*/, Qt::Orientation /*orientation*/, int /*role*/) const
{
	try
	{
	}
	CATCH_AND_REPORT
	return QVariant();
}

/*
 * (see Qt documentation)
 */
QVariant DatabaseModel::data(const QModelIndex& index, int role) const
{
	try
	{
		if(index.isValid())
		{
			AbstractDatabaseItem* item = indexToItem(index);
			switch(role)
			{
			case Qt::DisplayRole:
				return QString::fromStdString(item->getDisplayName());
			case Qt::EditRole:
				return QString::fromStdString(item->getName());
			case Qt::DecorationRole:
				{
					ConnectionItem* conn = dynamic_cast<ConnectionItem*>(item);
					if(conn != 0)
					{
						if(conn->isConnected())
						{
							return connectedIcon;
						}
						else
						{
							return disconnectedIcon;
						}
					}
					break;
				}
			case TypeRole:
				return item->getType();
			case Qt::BackgroundRole:
				{
					VersionItem* version = dynamic_cast<VersionItem*>(item);
					if(version != 0)
					{
						if(version->getStructureColor())
						{
							return QBrush(QColor(0xAF, 0xE9, 0xDD));
						}
						else
						{
							return QBrush(QColor(0xD7, 0xF4, 0xEE));
						}
					}
					SiteItem* site = dynamic_cast<SiteItem*>(item);
					if (site != 0)
					{
						return QBrush(QColor(0xFF, 0xFA, 0xCD)); // LemonChiffon
					}
					RunItem* run = dynamic_cast<RunItem*>(item);
					if (run != 0)
					{
						return QBrush(QColor(0xB0, 0xE2, 0xFF)); // LightSkyBlue1
					}
					ConditionItem* condition = dynamic_cast<ConditionItem*>(item);
					if (condition != 0)
					{
						return QBrush(QColor(0xD7, 0xF4, 0xEE));
					}
					break;
				}
			}
		}
	}
	CATCH_AND_REPORT
	return QVariant();
}

/*
 * (see Qt documentation)
 */
bool DatabaseModel::setData(const QModelIndex& index, const QVariant& value, int role)
{
	try
	{
		if(index.isValid())
		{
			if(role == Qt::EditRole)
			{
				AbstractDatabaseItem* item = indexToItem(index);
				item->setName(value.toByteArray().data());
				emit dataChanged(index, index);

				// If a version was renamed, update display name of associated run conditions
				if (item->getType() == AbstractDatabaseItem::Version)
				{
					QModelIndex schemaIndex = index.parent().parent();
					QModelIndexList conditions = match(schemaIndex, TypeRole, AbstractDatabaseItem::RunCondition, -1, Qt::MatchRecursive);
					for (int i=0; i < conditions.size(); ++i)
					{
						QModelIndex conditionIndex = conditions.at(i);
						ConditionItem* condition = dynamic_cast<ConditionItem*>(indexToItem(conditionIndex));
						condition->refresh();
						emit dataChanged(conditionIndex, conditionIndex);
					}
				}
			}
		}
	}
	CATCH_AND_REPORT
	return false;
}

/*
 * (see Qt documentation)
 */
bool DatabaseModel::removeRows(int row, int count, const QModelIndex& parent)
{
	try
	{
		if(count > 0)
		{
			AbstractDatabaseItem* item = indexToItem(parent);
			for(int i = 0; i < count; ++i)
			{
				item->removeItem(row + i, AbstractDatabaseItem::FromDataStore);
			}
			// removal from the data store can fail, so we have to do it
			// separately, otherwise there would be exception safety problems.
			// As Qt cannot rollback beginRemoveRows, the item would disappear
			// temporarily, even if we execute endRemoveRows during stack unwinding.
			beginRemoveRows(parent, row, row + count - 1);
			for(int i = 0; i < count; ++i)
			{
				item->removeItem(row + i, AbstractDatabaseItem::FromList);
			}
			endRemoveRows();
			// If a version was removed, update all runs
			if (item->getType() == AbstractDatabaseItem::Config)
			{
				QModelIndex schemaIndex = parent.parent();
				QModelIndexList runs = match(schemaIndex, TypeRole, AbstractDatabaseItem::Run, -1, Qt::MatchRecursive);
				for (int i=0; i < runs.size(); ++i)
				{
					refresh(runs.at(i));
				}
			}
		}
	}
	CATCH_AND_REPORT
	return false;
}

///*
// * (see Qt documentation)
// */
//bool DatabaseModel::removeRows(int row, int count, const QModelIndex& parent)
//{
//
//}

/**
 * Retrieves the database item corresponding to the specified index.
 */
AbstractDatabaseItem* DatabaseModel::indexToItem(const QModelIndex& index)
{
	assert(index.isValid());
	AbstractDatabaseItem* item = reinterpret_cast<AbstractDatabaseItem*>(index.internalPointer());
	assert(item != 0);
	return item;
}

/**
 * Creates an index for the specified database item.
 */
QModelIndex DatabaseModel::itemToIndex(AbstractDatabaseItem* item) const
{
	if(item == 0)
	{
		return QModelIndex();
	}
	else
	{
		return createIndex(getRowInParent(item), 0, item);
	}
}

/**
 * Returns the row of the specified item in its parent.
 */
int DatabaseModel::getRowInParent(AbstractDatabaseItem* item) const
{
	int row;
	AbstractDatabaseItem* parent = item->getParent();
	if(parent != 0)
	{
		size_t index = parent->searchItem(item);
		if(index == AbstractDatabaseItem::INVALID_INDEX)
		{
			row = -1;
		}
		else
		{
			row = index;
		}
	}
	else
	{
		row = 0; // the row of the root is 0.
	}
	return row;
}

}
}
