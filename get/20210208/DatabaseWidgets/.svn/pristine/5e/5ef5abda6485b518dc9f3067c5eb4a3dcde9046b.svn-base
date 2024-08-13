#ifndef CCFG_MODEL_ABSTRACTDATABASEITEM_H_INCLUDED
#define CCFG_MODEL_ABSTRACTDATABASEITEM_H_INCLUDED

#include <QFlags>
#include <string>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for all kinds of item that should be displayed in the
 * database browser tree.
 *
 * @author fchateau
 */
class AbstractDatabaseItem
{
public:
	static const size_t INVALID_INDEX = 0xFFFFFFFF;

	enum ItemType {Root = 0, Connection, Schema, Config, Version, Site, Run, RunCondition, TypeCount};
	enum RemoveOption {FromList = 0x1, FromDataStore = 0x2};
	Q_DECLARE_FLAGS(RemoveOptions, RemoveOption)

	AbstractDatabaseItem();
	virtual ~AbstractDatabaseItem();

	size_t searchItem(AbstractDatabaseItem* item);

	// @{ interface
	virtual bool isLeaf() const;
	virtual bool canSetName() const = 0;
	virtual void setName(const std::string& newName) = 0;
	virtual const std::string& getName() const = 0;
	virtual const std::string& getDisplayName() const;
	virtual size_t getItemCount() = 0;
	virtual AbstractDatabaseItem* getItem(size_t index) = 0;
	virtual AbstractDatabaseItem* getParent() = 0;
	virtual void refresh() = 0;
	virtual void removeItem(size_t index, RemoveOptions options = FromList) = 0;
	virtual void clear() = 0;
	virtual std::string buildUrl() const = 0;
	virtual ItemType getType() const = 0;
	// @}

private:
	AbstractDatabaseItem(const AbstractDatabaseItem& rvalue);
	AbstractDatabaseItem& operator=(const AbstractDatabaseItem& rvalue);
};

Q_DECLARE_OPERATORS_FOR_FLAGS(AbstractDatabaseItem::RemoveOptions)

}
}

#endif // CCFG_MODEL_ABSTRACTDATABASEITEM_H_INCLUDED
