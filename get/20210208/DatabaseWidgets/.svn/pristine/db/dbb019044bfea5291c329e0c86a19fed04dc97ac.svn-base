#ifndef CCFG_MODEL_RUNITEM_H_INCLUDED
#define CCFG_MODEL_RUNITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/Db/Run.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace CCfg
{
namespace Model
{
// forward declarations
class SiteItem;
class ConditionItem;
class VersionItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a site in the database tree.
 *
 * @author psizun
 */
class RunItem : public AbstractDatabaseItem
{
	SiteItem* parent;
	Io::Db::RunPtr run;
	boost::ptr_vector<ConditionItem> conditions;
	std::string name;
	bool upToDate;

public:
	RunItem(SiteItem* parent, const Io::Db::RunPtr& run);
	~RunItem();

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

	SiteItem* getSite() const { return parent; }
	Io::Db::Run& getEntity() { return *run; }
	void createCondition(const VersionItem*);
protected:
	void ensureUpToDate();
	void updateDisplayName();

private:
	RunItem(const RunItem& rvalue);
	RunItem& operator=(const RunItem& rvalue);
};

}
}

#endif // CCFG_MODEL_RUNITEM_H_INCLUDED
