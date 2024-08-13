#ifndef CCFG_MODEL_CONNECTIONITEM_H_INCLUDED
#define CCFG_MODEL_CONNECTIONITEM_H_INCLUDED

#include "AbstractDatabaseItem.h"
#include "CCfg/Io/DatabaseStorage.h"
#include <boost/ptr_container/ptr_vector.hpp>

namespace CCfg
{
namespace Model
{
// forward declarations
class ConnectionRootItem;
class SchemaItem;

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a registered connection in the database tree. It can
 * enumerate the schemas of the database.
 *
 * @author fchateau
 */
class ConnectionItem : public AbstractDatabaseItem
{
	enum State {Initial, Connected, Disconnected, Failed, Canceled};

	ConnectionRootItem* parent;
	std::string name;
	Io::DatabaseStorage storage;
	boost::ptr_vector<SchemaItem> schemas;
	State state;
	bool dirty;

public:
	ConnectionItem();
	ConnectionItem(const std::string& name, Io::DbmsType dbmsType,
			const std::string& host, const std::string& port,
			const std::string& user, const std::string& pwd,
			const std::string& sid, const std::string& schema);
	~ConnectionItem();

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

	// @{
	void connect();
	void disconnect();
	bool isConnected() const;
	// @}

	void createSchema(const std::string& name, const std::string& password);

	// @{ accessors
	void setParent(ConnectionRootItem* newParent) { parent = newParent; }
	void setLoginCallback(const Io::LoginCallbackPtr& callback);
	ConnectionRootItem* getRoot() { return parent; }
	Io::DbmsType getDbmsType() const { return storage.getDbmsType(); }
	const std::string& getUser() const { return storage.getUser(); }
	const std::string& getPassword() const { return storage.getPassword(); }
	const std::string& getHost() const { return storage.getHost(); }
	const std::string& getPort() const { return storage.getPort(); }
	const std::string& getSID() const { return storage.getSID(); }
	const std::string& getSchema() const { return storage.getSchema(); }
	void setDbmsType(Io::DbmsType newDbms);
	void setUser(const std::string& newUser);
	void setPassword(const std::string& newPasswd);
	void setHost(const std::string& newHost);
	void setPort(const std::string& newPort);
	void setSID(const std::string& newSID);
	void setSchema(const std::string& newSchema);
	void setCurrentConfigName(const std::string& curConfigName);
	void setCurrentConfigDate(const std::tm& curConfigDate);
	// @}

protected:
	void ensureUpToDate();
	void rebuildSchemaList();

private:
	ConnectionItem(const ConnectionItem& rvalue);
	ConnectionItem& operator=(const ConnectionItem& rvalue);
};

}
}

#endif // CCFG_MODEL_CONNECTIONITEM_H_INCLUDED
