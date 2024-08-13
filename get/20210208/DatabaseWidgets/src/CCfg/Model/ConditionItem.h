#ifndef CCFG_MODEL_RUNCONDITIONITEM_H_INCLUDED
#define CCFG_MODEL_RUNCONDITIONITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/Db/RunCondition.h"

namespace CCfg
{
namespace Model
{
// forward declarations
class RunItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a configuration version in the database tree.
 *
 * @author fchateau
 */
class ConditionItem : public AbstractDatabaseItem
{
	RunItem* parent;
	Io::Db::RunConditionPtr runCondition;
	std::string name;

public:
	ConditionItem(RunItem* parent, const Io::Db::RunConditionPtr& condition);
	~ConditionItem();

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

	RunItem* getRun() const { return parent; }
	Io::Db::RunCondition& getEntity() { return *runCondition; }

protected:
	void updateDisplayName();

private:
	ConditionItem(const ConditionItem& rvalue);
	ConditionItem& operator=(const ConditionItem& rvalue);
};

}
}

#endif // CCFG_MODEL_RUNCONDITIONITEM_H_INCLUDED
