#ifndef CCFG_MODEL_CONFIGMODEL_H_INCLUDED
#define CCFG_MODEL_CONFIGMODEL_H_INCLUDED

#include "AbstractConfigModelAccessor.h"
#include "ConfigNode.h"
#include <QAbstractItemModel>
#include <memory>
#include <vector>

// forward declarations
class QUndoStack;
class QDir;

namespace CCfg
{
// forward declarations
class Attribute;
class Path;
namespace Io {
	class Document;
	class Url;
}

namespace Model
{
// forward declarations
struct AbstractColumn;
class ConfigNode;
class UndoableConfigModelAccessor;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * The model communicating with QTreeView to display a configuration tree.
 *
 * @author fchateau
 */
class ConfigModel : public QAbstractItemModel
{
	Q_OBJECT
public:
	/// Generates an index for each column, and gives it a user friendly name.
	enum Columns
	{
		NameColumn,	///< the index of the column displaying the node name.
		ValueColumn,///< the index of the column displaying the node value
		UnitColumn,	///< the index of the column displaying the node unit.
		RangeColumn,///< the index of the column displaying the node range.
		BacktrackColumn,///< the index of the column displaying the backtrack level.
		IncludeColumn,	///< the index of the column displaying the include path.
		COLUMN_COUNT	///< the number of column in the model.
	};

protected:
	/// the list of column object defining column-specific rendering behaviour.
	std::auto_ptr<AbstractColumn> columns[COLUMN_COUNT];

	/// the object through which modifications on the model are carried out.
	UndoableConfigModelAccessor* accessor;

public:
	ConfigModel(QObject* parent = 0);
	~ConfigModel();

	// @{ Configuration document getter and setter
	void setDocument(std::auto_ptr<Io::Document> newDocument);
	const Io::Document* getDocument() const;

	const std::string& getDocumentName() const;
	const std::string& getDocumentUrl() const;
	static Io::Url getQualifiedIncludeUrl(const ConfigNode* node);
	static QDir getConfigDirectory(const QModelIndex& index);
	// @}

	// @{ Misc
	QUndoStack* getUndoStack();
	UndoableConfigModelAccessor& getAccessor() const { return *accessor; }
	// @}

	// @{ Save related methods
	void save();
	void saveAs(const Io::Url& url);
	// @}

	// @{ Index, path and node conversions
	static ConfigNode* indexToNode(const QModelIndex& index);
	QModelIndex nodeToIndex(ConfigNode* node, int column = 0) const;
	QModelIndex pathToIndex(const Path& path, int column = 0) const;
	bool isDescendant(const QModelIndex& index, const QModelIndex& ancestor) const;
	// @}

	// columns
	AbstractColumn& getColumn(int columnIndex);

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

	// @{ Mutation methods
	void setName(const QModelIndex& nodeIndex, const std::string& newName);
	void insertAfter(const QModelIndex& nodeIndex, std::auto_ptr<Attribute> newAttribute);
	void insertInside(const QModelIndex& parentIndex, std::auto_ptr<Attribute> newAttribute);
	void override(const QModelIndex& nodeIndex);
	void remove(const QModelIndex& parent, IndexList& indexes);
	void include(const QModelIndex& destNode, const std::string& url);
	void import(const QModelIndex& destNode, const std::string& url);
	void exportSubConfig(const QModelIndex& rootObject, const std::string& url);
	// @}

signals:
	// @{ mutation notifications
	void documentUrlChanged(const std::string& newUrl, const std::string& oldUrl);
	void beginInsertion(ConfigNode* parent);
	void endInsertion(ConfigNode* parent, const IndexList& indexes);
	void beginRemoval(ConfigNode* parent, const IndexList& indexes);
	void endRemoval(ConfigNode* parent);
	void saved(const QString& url);
	// @}

	// @{ general update notification
	void beginUpdate();
	void endUpdate();
	// @}

private:
	// non-copyable
	ConfigModel(const ConfigModel& value);
	ConfigModel& operator=(const ConfigModel& value);

private slots:
	void doSetData(const QModelIndex& index, const QVariant& value, int role);

	// @{
	void onBeginClear(ConfigNode* node);
	void onEndClear();
	void onBeginRecompute(ConfigNode* node);
	void onEndRecompute();
	// @}
};

}
}

#endif // CCFG_MODEL_CONFIGMODEL_H_INCLUDED
