#include "ConfigModel.h"
#include "UndoableConfigModelAccessor.h"
#include "ConfigModelAccessor.h"
#include "ConfigNode.h"
#include "Column/AttributeName.h"
#include "Column/AttributeValue.h"
#include "Column/AttributeUnit.h"
#include "Column/AttributeRange.h"
#include "Column/AttributeBacktrack.h"
#include "Column/AttributeInclude.h"
#include "Column/ConfigTreeRow.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/Url.h"
#include "CCfg/Path.h"
#include "CCfg/CConfig.h"
#include "CCfg/Object.h"
#include "Utils/Error/ExcepHandling.h"
#include <QVariant>
#include <QUndoStack>
#include <QDir>
#include <QFont>

namespace CCfg
{
namespace Model
{
using Column::ConfigTreeRow;

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModel class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent the owner of this new object.
 */
ConfigModel::ConfigModel(QObject* parent)
	: QAbstractItemModel(parent),
	accessor(new UndoableConfigModelAccessor(this))
{
	columns[NameColumn].reset(new Column::AttributeName(accessor));
	columns[ValueColumn].reset(new Column::AttributeValue(accessor));
	columns[UnitColumn].reset(new Column::AttributeUnit(accessor));
	columns[RangeColumn].reset(new Column::AttributeRange(accessor));
	columns[BacktrackColumn].reset(new Column::AttributeBacktrack());
	columns[IncludeColumn].reset(new Column::AttributeInclude(accessor));

	ConfigModelAccessor* m = accessor->getConfigModelAccessor();
	connect(m, SIGNAL(beginInsertion(ConfigNode*)),
		this, SIGNAL(beginInsertion(ConfigNode*)));
	connect(m, SIGNAL(endInsertion(ConfigNode*, IndexList)),
		this, SIGNAL(endInsertion(ConfigNode*, IndexList)));
	connect(m, SIGNAL(beginRemoval(ConfigNode*, IndexList)),
		this, SIGNAL(beginRemoval(ConfigNode*, IndexList)));
	connect(m, SIGNAL(endRemoval(ConfigNode*)),
		this, SIGNAL(endRemoval(ConfigNode*)));

	connect(m, SIGNAL(beginClear(ConfigNode*)),
		this, SLOT(onBeginClear(ConfigNode*)));
	connect(m, SIGNAL(endClear()),
		this, SLOT(onEndClear()));
	connect(m, SIGNAL(beginRecompute(ConfigNode*)),
		this, SLOT(onBeginRecompute(ConfigNode*)));
	connect(m, SIGNAL(endRecompute()),
		this, SLOT(onEndRecompute()));
	connect(m, SIGNAL(beginUpdate()),
		this, SIGNAL(beginUpdate()));
	connect(m, SIGNAL(endUpdate()),
		this, SIGNAL(endUpdate()));

	int type = QMetaType::type("QModelIndex");
	if(type == 0 or not QMetaType::isRegistered(type))
	{
		qRegisterMetaType<QModelIndex>("QModelIndex");
	}
	type = QMetaType::type("QVariant");
	if(type == 0 or not QMetaType::isRegistered(type))
	{
		qRegisterMetaType<QVariant>("QVariant");
	}
}

/**
 * Destructor.
 */
ConfigModel::~ConfigModel()
{
}

/**
 * Setter of the "document" property.
 * Redefines the document that the model represents.
 * @param newDocument the new document for the model.
 */
void ConfigModel::setDocument(std::auto_ptr<Io::Document> newDocument)
{
	accessor->setDocument(newDocument);
}

/**
 * Getter of the "document" property.
 * @return Returns the document that the model represents.
 */
const Io::Document* ConfigModel::getDocument() const
{
	return accessor->getDocument();
}

/**
 * Returns the document name.
 */
const std::string& ConfigModel::getDocumentName() const
{
	return getDocument()->getName();
}

/**
 * Returns the document URL (the place where it is saved).
 */
const std::string& ConfigModel::getDocumentUrl() const
{
	return getDocument()->getUrl();
}

/**
 * Returns the qualified include URL of the specified node.
 * @param node the node from which we want to retrieve the qualified include URL.
 * It must be a configuration object.
 * @throws std::bad_cast if node is not a configuration object.
 */
Io::Url ConfigModel::getQualifiedIncludeUrl(const ConfigNode* node)
{
	return Io::Document::computeQualifiedIncludeUrl(node->getConfigAttribute()->toObject());
}

/**
 * Returns the directory of the configuration where the specified item belongs.
 */
QDir ConfigModel::getConfigDirectory(const QModelIndex& index)
{
	const ConfigNode* node = ConfigModel::indexToNode(index);
	const ConfigNode* includeNode = node->getIncludeRoot();

	Io::Url qualifiedUrl = getQualifiedIncludeUrl(includeNode);

	QDir configDir(qualifiedUrl.toPath().c_str());
	configDir.makeAbsolute();	// ensure the internal path is absolute (otherwise cdUp may fail)
	configDir.cdUp(); // remove the file part.
	return configDir;
}

/**
 * Returns this model's undo stack.
 * @return Returns this model's undo stack.
 */
QUndoStack* ConfigModel::getUndoStack()
{
	return accessor->getUndoStack();
}

/**
 * Saves the model document using its current URL.
 * The URL where the configuration is saved can be obtained with getDocumentUrl()
 * and changed with setDocumentUrl().
 */
void ConfigModel::save()
{
	accessor->getDocument()->save();
	getUndoStack()->setClean();
	emit saved(QString::fromStdString(getDocumentUrl()));
}

/**
 * Saves the model document using a new URL.
 * This is a convenience method: the same result can be achieved with
 * setDocumentUrl() and save().
 */
void ConfigModel::saveAs(const Io::Url& url)
{
	std::string newUrl = url.format();
	std::string oldUrl = getDocumentUrl();
	accessor->getDocument()->saveToUrl(url, false);
	getUndoStack()->setClean();
	if(newUrl != oldUrl)
	{
		emit documentUrlChanged(newUrl, oldUrl);
	}
	emit saved(QString::fromStdString(newUrl));
}

/**
 * Provides the ConfigNode designated by a QModelIndex.
 * @param index the index for which the ConfigNode is requested.
 * @return Returns the ConfigNode designated by a QModelIndex.
 */
ConfigNode* ConfigModel::indexToNode(const QModelIndex& index)
{
	assert(index.isValid());
	ConfigNode* node = static_cast<ConfigNode*>(index.internalPointer());
	assert(node != 0);
	return node;
}

/**
 * Creates a QModelIndex designating a configuration node.
 * @param node the configuration node for which we need the index.
 * @param column the column that the index must target.
 * @return Returns the QModelIndex designating a configuration node.
 */
QModelIndex ConfigModel::nodeToIndex(ConfigNode* node, int column) const
{
	if(node == 0)
	{
		return QModelIndex();
	}
	else
	{
		return createIndex(node->getRowInParent(), column, node);
	}
}

/**
 * Creates a QModelIndex for the node matching a path.
 * @param path the path to search in the model.
 * @param column the column that the index must target.
 * @return Returns the QModelIndex for the node matching the path.
 */
QModelIndex ConfigModel::pathToIndex(const Path& path, int column) const
{
	ConfigNode* node = accessor->pathToNode(path);
	if(node != 0)
	{
		return nodeToIndex(node, column);
	}
	else
	{
		return QModelIndex();
	}
}

/**
 * Determines whether an item is the descendant of another.
 */
bool ConfigModel::isDescendant(const QModelIndex& index, const QModelIndex& ancestor) const
{
	QModelIndex parent = index.parent();
	while(parent.isValid() and parent != ancestor)
	{
		parent = parent.parent();
	}
	return parent == ancestor;
}

/**
 * Getter of the "columns" property.
 * @param columnIndex the index of the column to retrieve.
 * @return Returns a reference to the column of the model.
 */
AbstractColumn& ConfigModel::getColumn(int columnIndex)
{
	assert(columnIndex >= NameColumn and columnIndex < COLUMN_COUNT);
	return *columns[columnIndex].get();
}

/**
 * Provides the index of the children of the specified item.
 * @see QAbstractItemModel::index(int, int, const QModelIndex&) const
 */
QModelIndex ConfigModel::index(int row, int column, const QModelIndex& parent) const
{
	try
	{
		if(not parent.isValid()) // the configuration root is requested
		{
			if(row != 0 or accessor->getRootNode() == 0)
			{
				// configuration root is unique,
				// requesting a different index is invalid.
				return QModelIndex();
			}
			else
			{
				// We create an index for the configuration root node.
				return createIndex(row, column, accessor->getRootNode());
			}
		}
		else
		{
			ConfigNode* node = indexToNode(parent);
			std::vector<ConfigNode>& children = node->getChildren();

			// check whether the request row and column are valid.
			if(row < 0 or row >= static_cast<int>(children.size()) or
				column < 0 or column >= COLUMN_COUNT)
			{
				return QModelIndex();
			}
			else
			{
				// create an index. Pass the corresponding child as internal pointer.
				return createIndex(row, column, &children[row]);
			}
		}
	}
	CATCH_AND_REPORT
	return QModelIndex();
}

/**
 * Provides the index of the parent of the specified item.
 * @see QAbstractItemModel::parent(const QModelIndex&) const
 */
QModelIndex ConfigModel::parent(const QModelIndex& index) const
{
	try
	{
		if(index.isValid())
		{
			ConfigNode* node = indexToNode(index);
			if(node != accessor->getRootNode())
			{
				assert(node->getParent() != 0);
				return createIndex(node->getParent()->getRowInParent(), NameColumn, node->getParent());
			}
			// else the parent index of the configuration root is the invalid index.
		}
	}
	CATCH_AND_REPORT
	return QModelIndex();
}

/**
 * Tells whether an item has children or not.
 * @see QAbstractItemModel::hasChildren(const QModelIndex&) const
 */
bool ConfigModel::hasChildren(const QModelIndex& parent) const
{
	try
	{
		if(not parent.isValid())
		{
			// check if the configuration root exists.
			return accessor->getRootNode() != 0;
		}
		else
		{
			ConfigNode* node = indexToNode(parent);
			return node->isObject();
		}
	}
	CATCH_AND_REPORT
	return false;
}

/**
 * Provides the row count of an item.
 * @see QAbstractItemModel::rowCount(const QModelIndex&) const
 */
int ConfigModel::rowCount(const QModelIndex& parent) const
{
	try
	{
		if(not parent.isValid())
		{
			// if there is no configuration root return 0, return 1 otherwise.
			return accessor->getRootNode() != 0;
		}
		else
		{
			ConfigNode* node = indexToNode(parent);
			return node->getChildren().size();
		}
	}
	CATCH_AND_REPORT
	return 0;
}

/**
 * Provides the column count of an item.
 * @see QAbstractItemModel::columnCount(const QModelIndex&) const
 */
int ConfigModel::columnCount(const QModelIndex& /*parent*/) const
{
	return COLUMN_COUNT;
}

/**
 * Provides the flags for the specified tree index.
 * @see QAbstractItemModel::flag(const QModelIndex&) const
 */
Qt::ItemFlags ConfigModel::flags(const QModelIndex& index) const
{
	try
	{
		if(index.isValid())
		{
			ConfigNode* node = indexToNode(index);
			ConfigTreeRow row(node);
			return columns[index.column()]->getRowFlags(row);
		}
	}
	CATCH_AND_REPORT
	return Qt::ItemFlags();
}

/**
 * Provides the column or row header data.
 * @see QAbstractItemModel::headerData(int, Qt::Orientation, int) const
 */
QVariant ConfigModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	QVariant data;
	try
	{
		if(orientation == Qt::Horizontal)
		{
			switch(role)
			{
			case Qt::DisplayRole:
				data = columns[section]->getHeaderDisplayName();
				break;
			case Qt::ToolTipRole:
				data = columns[section]->getHeaderTooltip();
				break;
			default:
				data = QAbstractItemModel::headerData(section, orientation, role);
				break;
			}
		}
		else
		{
			data = QAbstractItemModel::headerData(section, orientation, role);
		}
	}
	CATCH_AND_REPORT
	return data;
}

/**
 * Provides the data for the specified index and role.
 * @see QAbstractItemModel::data(const QModelIndex&, int) const
 */
QVariant ConfigModel::data(const QModelIndex& index, int role) const
{
	try
	{
		if(index.isValid())
		{
			ConfigNode* node = indexToNode(index);
			ConfigTreeRow row(node);
			AbstractColumn* column = columns[index.column()].get();

			switch(role)
			{
			case Qt::DisplayRole:
				return column->getRowDisplayData(row);
			case Qt::EditRole:
				return column->getRowEditionData(row);
			case Qt::ToolTipRole:
				return column->getRowTooltip(row);
			case Qt::FontRole:
				return column->getRowFont(row);
			case Qt::TextAlignmentRole:
				return static_cast<int>(column->getRowAlignment(row));
			case Qt::ForegroundRole:
				return column->getRowForegroundTextColor(row);
			case Qt::BackgroundRole:
				return column->getRowBackgroundColor(row);
			default:
				break;
			}
		}
	}
	CATCH_AND_REPORT
	return QVariant();
}

/**
 * Modifies the specified item through the undo framework.
 * @param index the index of the item to modify.
 * @param value the new value for the item.
 * @param role the role of the data that should be modified.
 * @return Returns true if changed were accepted, false otherwise.
 */
bool ConfigModel::setData(const QModelIndex& index, const QVariant& value, int role)
{
	// this signal asynchronously trapped by the onDataSet slot.
	// It is necessary to do so because when setData is called, the editor
	// is not closed yet.  As edition can cause the row being edited to be destroyed
	// temporarily (which cause a crash), therefore edition must be delayed.
	QMetaObject::invokeMethod(this, "doSetData", Qt::QueuedConnection,
		Q_ARG(QModelIndex, index), Q_ARG(QVariant, value), Q_ARG(int, role));
	return true;
}

/**
 * Removes the specified rows from the ConfigModel and propagates update.
 * @param row the index of the row in the parent.
 * @param count the number of rows to removes from the specified "row" index.
 * @param parent the parent of the rows to remove.
 */
bool ConfigModel::removeRows(int row, int count, const QModelIndex& parent)
{
	try
	{
		IndexList indexes(count);
		for(int i = 0; i < count; ++i)
		{
			indexes[i] = row + i;
		}

		accessor->remove(indexToNode(parent), indexes);
		return true;
	}
	CATCH_AND_REPORT
	return false;
}

/**
 * Changes the name of a configuration node.
 * @param nodeIndex the index of the node to alter.
 * @param newName the new name to set.
 */
void ConfigModel::setName(const QModelIndex& nodeIndex, const std::string& newName)
{
	accessor->setName(indexToNode(nodeIndex), newName);
}

/**
 * Inserts an attribute after a node.
 * @param nodeIndex the node after which the new attribute is inserted.
 * @param newAttribute the attribute to insert. The configuration tree takes ownership of it.
 */
void ConfigModel::insertAfter(const QModelIndex& nodeIndex, std::auto_ptr<Attribute> newAttribute)
{
	accessor->insertAfter(indexToNode(nodeIndex), newAttribute);
}

/**
 * Inserts an attribute as the first child of a node.
 * @param parentIndex the parent for the new attribute.
 * @param newAttribute the attribute to insert. The configuration tree takes ownership of it.
 */
void ConfigModel::insertInside(const QModelIndex& parentIndex, std::auto_ptr<Attribute> newAttribute)
{
	accessor->insertInside(indexToNode(parentIndex), newAttribute);
}

/**
 * Overrides an inherited attribute.
 * @param nodeIndex the index of the inherited attribute ConfigNode.
 */
void ConfigModel::override(const QModelIndex& nodeIndex)
{
	ConfigNode* node = indexToNode(nodeIndex);
	std::vector<ConfigNode*> created;
	accessor->override(node, created);
}

/**
 * Removes a list of ConfigNodes from a common parent.
 * @param parent the parent of the ConfigNodes to remove.
 * @param indexes the row indexes of the ConfigNodes to remove.
 * They will be sorted by this method (which in why it does not guarantee constness).
 */
void ConfigModel::remove(const QModelIndex& parent, IndexList& indexes)
{
	std::sort(indexes.begin(), indexes.end());
	accessor->remove(indexToNode(parent), indexes);
}

/**
 * Loads a configuration from an URL and include it in a destination node.
 * "Include" have the same semantics than in C: it means importing a configuration
 * by reference. The 2 configuration involved will stay separate, and the
 * destination will be dependant of the source, as long as this link exists.
 * The model will be constructed as the fusion of the 2 configurations, but if a
 * change is saved in the scope of the included file, it will modify the included
 * file, and not the includer.
 * @param destNode the index of the destination node (it must be a configuration
 * object). The content of the node will be erased (but undoing is possible).
 * @param url the URL of the configuration to load.
 */
void ConfigModel::include(const QModelIndex& destNode, const std::string& url)
{
	accessor->include(indexToNode(destNode), url);
}

/**
 * Loads a configuration from an URL and insert it
 */
void ConfigModel::import(const QModelIndex& destNode, const std::string& url)
{
	accessor->import(indexToNode(destNode), url);
}

/**
 * Saves a sub configuration in a location designated by an URL.
 */
void ConfigModel::exportSubConfig(const QModelIndex& rootObject, const std::string& url)
{
	getDocument()->exportSubConfig(
			indexToNode(rootObject)->getConfigAttribute()->toObject(),
			Io::Url(url));
}

/**
 * Slot that asynchronously process column edition.
 * @param index the index the of node column to modify.
 * @param value the new value to set.
 * @param role the role of the value to set.
 */
void ConfigModel::doSetData(const QModelIndex& index, const QVariant& value, int role)
{
	try
	{
		if(index.isValid() and role == Qt::EditRole)
		{
			ConfigNode* node = indexToNode(index);
			ConfigTreeRow row(node);
			AbstractColumn* column = columns[index.column()].get();

			column->setRowEditionData(value.toString(), row);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the beginning of a node reinitialization.
 * @param node the node being cleared.
 */
void ConfigModel::onBeginClear(ConfigNode* node)
{
	try
	{
		QModelIndex index = nodeToIndex(node);
		beginRemoveRows(index, 0, rowCount(index) - 1);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the end of a node reinitialization.
 */
void ConfigModel::onEndClear()
{
	try
	{
		endRemoveRows();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the beginning of a node recomputation.
 * @param node the node being recomputed.
 */
void ConfigModel::onBeginRecompute(ConfigNode* node)
{
	try
	{
		QModelIndex index = nodeToIndex(node);
		if(index.isValid())
		{
			size_t childrenCount = node->getChildren().size();
			beginInsertRows(index, 0, childrenCount-1);
		}
		else
		{
			beginInsertRows(index, 0, 0);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the end of a node recomputation.
 */
void ConfigModel::onEndRecompute()
{
	try
	{
		endInsertRows();
	}
	CATCH_AND_REPORT
}


}
}
