#include "ConnectionItem.h"
#include "SchemaItem.h"
#include "ConnectionRootItem.h"

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConnectionItem class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ConnectionItem::ConnectionItem() : parent(0), storage(), state(Initial), dirty(false)
{
}

/**
 * Full constructor.
 */
ConnectionItem::ConnectionItem(const std::string& name, Io::DbmsType dbmsType,
		const std::string& host, const std::string& port,
		const std::string& user, const std::string& pwd,
		const std::string& sid, const std::string& schema) :
	parent(0), name(name), storage(dbmsType, host, port, user, pwd, sid, schema),
	state(Initial), dirty(false)
{
}

/**
 * Destructor.
 */
ConnectionItem::~ConnectionItem()
{
}

/*
 * (see inherited doc)
 */
bool ConnectionItem::canSetName() const
{
	return true;
}

/*
 * (see inherited doc)
 */
void ConnectionItem::setName(const std::string& newName)
{
	name = newName;
	refresh();
}

/*
 * (see inherited doc)
 */
const std::string& ConnectionItem::getName() const
{
	return name;
}

/*
 * (see inherited doc)
 */
size_t ConnectionItem::getItemCount()
{
	ensureUpToDate();
	return schemas.size();
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConnectionItem::getItem(size_t index)
{
	return &schemas[index];
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConnectionItem::getParent()
{
	return parent;
}

/*
 * (see inherited doc)
 */
void ConnectionItem::refresh()
{
	dirty = true;
}

/*
 * (see inherited doc)
 */
void ConnectionItem::removeItem(size_t index, RemoveOptions options)
{
	if(options & FromDataStore)
	{
		storage.getDao().removeDatabaseAndUser(schemas[index].getName());
	}
	if(options & FromList)
	{
		schemas.erase(schemas.begin() + index);
	}
}

/*
 * (see inherited doc)
 */
void ConnectionItem::clear()
{
	schemas.clear();
}

/*
 * (see inherited doc)
 */
std::string ConnectionItem::buildUrl() const
{
	return Io::DatabaseStorage::formatUrl(getDbmsType(),
			getUser(),
			std::string(),
			getHost(),
			getPort(),
			getSID(),
			getSchema(),
			std::string());
}

/*
 * (see inherited doc)
 */
ConnectionItem::ItemType ConnectionItem::getType() const
{
	return Connection;
}

/**
 * Connects to the server using the parameters of this object.
 * This method has no effect if the connection is already established.
 */
void ConnectionItem::connect()
{
	if(state != Connected)
	{
		try
		{
			storage.open();
		}
		catch(...)
		{
			state = Failed;
			throw;
		}
		if(storage.isOpened())
		{
			state = Connected;
			dirty = true;
		}
		else
		{
			state = Canceled;
		}
	}
}

/**
 * Disconnects from the server.
 * This method has no effect if the connection is not established.
 */
void ConnectionItem::disconnect()
{
	if(state == Connected)
	{
		storage.close();
		state = Disconnected;
	}
}

/**
 * Returns whether this item has established the connection to its server.
 */
bool ConnectionItem::isConnected() const
{
	return storage.isOpened();
}

/**
 * Creates a schema on the server where we are connected.
 */
void ConnectionItem::createSchema(const std::string& name, const std::string& password)
{
	if(state == Connected)
	{
		Io::Db::ConfigDao& dao = storage.getDao();
		dao.createDatabaseAndUser(name, password);
		dao.setWorkingDatabase(name);
		dao.createAllTables();
	}
}

/**
 * Changes the callback used to request login informations when establishing connections.
 */
void ConnectionItem::setLoginCallback(const Io::LoginCallbackPtr& callback)
{
	storage.setLoginCallback(callback);
}

void ConnectionItem::setDbmsType(Io::DbmsType newDbms)
{
	storage.setDbmsType(newDbms);
}

void ConnectionItem::setUser(const std::string& newUser)
{
	storage.setUser(newUser);
}

void ConnectionItem::setPassword(const std::string& newPasswd)
{
	storage.setPassword(newPasswd);
}

void ConnectionItem::setHost(const std::string& newHost)
{
	storage.setHost(newHost);
}

void ConnectionItem::setPort(const std::string& newPort)
{
	storage.setPort(newPort);
}

void ConnectionItem::setSID(const std::string& newSID)
{
	storage.setSID(newSID);
}

void ConnectionItem::setSchema(const std::string& newSchema)
{
	storage.setSchema(newSchema);
}

void ConnectionItem::setCurrentConfigName(const std::string& curConfigName)
{
	storage.setCurrentConfigName(curConfigName);
}

void ConnectionItem::setCurrentConfigDate(const std::tm& curConfigDate)
{
	storage.setCurrentConfigDate(curConfigDate);
}

void ConnectionItem::ensureUpToDate()
{
	if(state == Initial)
	{
		connect();
	}
	if(state == Connected and dirty)
	{
		rebuildSchemaList();
		dirty = false;
	}
}

void ConnectionItem::rebuildSchemaList()
{
	using namespace Io::Db;
	schemas.clear();
	ConfigDao& dao = storage.getDao();
	ConfigDatabasePtrSetPtr databases = dao.getDatabases();
	ConfigDatabaseNameIdx& nameIdx = databases->get<Name>();
	for(ConfigDatabaseNameIdx::iterator curDb = nameIdx.begin(), endDb = nameIdx.end();
		curDb != endDb; ++curDb)
	{
		schemas.push_back(new SchemaItem(this, *curDb));
	}
}

}
}
