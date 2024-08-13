#ifndef CCFG_MODEL_CONNECTIONROOTITEM_H_INCLUDED
#define CCFG_MODEL_CONNECTIONROOTITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/DatabaseStorage.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include <memory>

namespace CCfg
{
namespace Model
{
// forward declarations
class ConnectionItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class is the root of the database model tree. It can enumerate
 * registered connections.
 *
 * @author fchateau
 */
class ConnectionRootItem : public AbstractDatabaseItem
{
	Io::LoginCallbackPtr loginCallback; // the callback responding to login requests when establishing connections.
	boost::ptr_vector<ConnectionItem> connections;
	bool upToDate;

public:
	ConnectionRootItem();
	~ConnectionRootItem();

	// @{ Interface implementation
	bool canSetName() const;
	void setName(const std::string& newName);
	const std::string& getName() const;
	size_t getItemCount();
	AbstractDatabaseItem* getItem(size_t index);
	AbstractDatabaseItem* getParent();
	void removeItem(size_t index, RemoveOptions options = FromList);
	void clear();
	void refresh();
	std::string buildUrl() const;
	ItemType getType() const;
	// @}

	void setLoginCallback(const Io::LoginCallbackPtr& callback) { loginCallback = callback; }
	void appendConnection(std::auto_ptr<ConnectionItem> conn);
	ConnectionItem* search(const std::string& connName);

protected:
	void ensureUpToDate();
	void loadConnections();
	void saveConnections();

private:
	ConnectionRootItem(const ConnectionRootItem& rvalue);
	ConnectionRootItem& operator=(const ConnectionRootItem& rvalue);

	static const char DBMS[];
	static const char HOST[];
	static const char PORT[];
	static const char USER[];
	static const char SID[];
	static const char SCHEMA[];
	static const char DEFAULT_DBMS[];
	static const char DEFAULT_HOST[];
	static const char DEFAULT_PORT[];
	static const char DEFAULT_USER[];
	static const char DEFAULT_SID[];
	static const char DEFAULT_SCHEMA[];
};

}
}

#endif // CCFG_MODEL_CONNECTIONROOTITEM_H_INCLUDED
