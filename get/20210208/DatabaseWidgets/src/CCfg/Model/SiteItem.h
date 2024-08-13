#ifndef CCFG_MODEL_SITEITEM_H_INCLUDED
#define CCFG_MODEL_SITEITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/Db/Site.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace CCfg
{
namespace Model
{
// forward declarations
class SchemaItem;
class RunItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a site in the database tree. It can enumerate
 * its runs.
 *
 * @author psizun
 */
class SiteItem : public AbstractDatabaseItem
{
	SchemaItem* parent;
	Io::Db::SitePtr site;
	boost::ptr_vector<RunItem> runs;
	bool upToDate;

public:
	SiteItem(SchemaItem* parent, const Io::Db::SitePtr& site);
	~SiteItem();

	void createRun();

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
	Io::Db::Site& getEntity() { return *site; }
	Io::Db::ConfigDao& getDao();
protected:
	void ensureUpToDate();

private:
	SiteItem(const SiteItem& rvalue);
	SiteItem& operator=(const SiteItem& rvalue);
};

}
}

#endif // CCFG_MODEL_SITEITEM_H_INCLUDED
