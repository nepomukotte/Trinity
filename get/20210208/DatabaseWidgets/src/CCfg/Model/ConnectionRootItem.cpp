#include "ConnectionRootItem.h"
#include "ConnectionItem.h"
#include "CCfg/Model/EditorConfiguration.h"
#include "Utils/Error/ExcepHandling.h"
#include "CCfg/CConfig.h"
#include <stdexcept>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConnectionRootItem class methods body
////////////////////////////////////////////////////////////////////////////////
const char ConnectionRootItem::DBMS[] = "dbms";
const char ConnectionRootItem::HOST[] = "host";
const char ConnectionRootItem::PORT[] = "port";
const char ConnectionRootItem::USER[] = "user";
const char ConnectionRootItem::SID[]  = "sid";
const char ConnectionRootItem::SCHEMA[] = "schema";
const char ConnectionRootItem::DEFAULT_DBMS[] = "invalid";
const char ConnectionRootItem::DEFAULT_HOST[] = "localhost";
const char ConnectionRootItem::DEFAULT_PORT[] = "0";
const char ConnectionRootItem::DEFAULT_USER[] = "";
const char ConnectionRootItem::DEFAULT_SID[]  = "";
const char ConnectionRootItem::DEFAULT_SCHEMA[] = "";

/**
 * Default constructor.
 */
ConnectionRootItem::ConnectionRootItem() : connections(), upToDate(false)
{
}

/**
 * Destructor.
 */
ConnectionRootItem::~ConnectionRootItem()
{
	try
	{
		saveConnections();
	}
	CATCH_AND_REPORT
}

/*
 * (see inherited doc)
 */
bool ConnectionRootItem::canSetName() const
{
	return false;
}

/*
 * (see inherited doc)
 */
void ConnectionRootItem::setName(const std::string& /*newName*/)
{
	throw std::runtime_error("cannot set the name of the connection root");
}

/*
 * (see inherited doc)
 */
const std::string& ConnectionRootItem::getName() const
{
	static std::string name("Connections list");
	return name;
}

/*
 * (see inherited doc)
 */
size_t ConnectionRootItem::getItemCount()
{
	ensureUpToDate();
	// we only call ensureUpToDate here, and not in getItem because we assume
	// that a call to getItemCount always precedes a call to getItem.
	return connections.size();
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConnectionRootItem::getItem(size_t index)
{
	return &connections[index];
}

/*
 * (see inherited doc)
 */
AbstractDatabaseItem* ConnectionRootItem::getParent()
{
	return 0;
}

/*
 * (see inherited doc)
 */
void ConnectionRootItem::removeItem(size_t index, RemoveOptions options)
{
	if(options & FromList)
	{
		connections.erase(connections.begin() + index);
	}
}

/*
 * (see inherited doc)
 */
void ConnectionRootItem::clear()
{
	connections.clear();
}

/*
 * (see inherited doc)
 */
void ConnectionRootItem::refresh()
{
	upToDate = false;
}

/*
 * (see inherited doc)
 */
std::string ConnectionRootItem::buildUrl() const
{
	return std::string();
}

/*
 * (see inherited doc)
 */
ConnectionRootItem::ItemType ConnectionRootItem::getType() const
{
	return Root;
}

/**
 * @fn void ConnectionRootItem::setLoginCallback(const Io::LoginCallbackPtr& callback)
 * Changes the login callback affected to newly created connections items.
 * It does not retroactively affects existing connection items.
 */

/**
 * Inserts the specified connection item at the end of the connection list.
 */
void ConnectionRootItem::appendConnection(std::auto_ptr<ConnectionItem> conn)
{
	conn->setLoginCallback(loginCallback);
	conn->setParent(this);
	connections.push_back(conn);
}

/**
 * Searches a connection of the specified name and return it, if found.
 */
ConnectionItem* ConnectionRootItem::search(const std::string& connName)
{
	for(boost::ptr_vector<ConnectionItem>::iterator curConn = connections.begin(),
			endConn = connections.end(); curConn != endConn; ++curConn)
	{
		if(curConn->getName() == connName)
		{
			return &*curConn;
		}
	}
	return 0;
}

/**
 * Ensures the list of connections is up to date.
 */
void ConnectionRootItem::ensureUpToDate()
{
	if(not upToDate)
	{
		loadConnections();
	}
}

/**
 * Loads the list of connections from the program configuration.
 */
void ConnectionRootItem::loadConnections()
{
	CConfig cfg(&EditorConfiguration::root());

	// calling create allows to ensure the connection list exists
	CConfig connList(&cfg("database")("connections").create(OBJECT_INST));
	CConfig::Iterator iConn = connList.iterate("connection");
	while(iConn.hasNext())
	{
		iConn.next();
		CConfig& connCfg = iConn.get();
		if(connCfg.getIndex() != "*")
		{
			std::auto_ptr<ConnectionItem> connItem(
				new ConnectionItem(connCfg.getIndex(),
					Io::DatabaseStorage::dbmsFromString(connCfg(DBMS)),
					connCfg(HOST),
					connCfg(PORT).getValueAsString(),
					connCfg(USER),
					"",
					connCfg(SID),
					connCfg(SCHEMA)
				)
			);
			connItem->setLoginCallback(loginCallback);
			appendConnection(connItem);
		}
	}
	upToDate = true;
}

/**
 * Saves the list of connection to the program configuration.
 */
void ConnectionRootItem::saveConnections()
{
	Attribute& cfg = EditorConfiguration::root();
	Attribute& connList = cfg("database")("connections");

	// update the default connection
	Attribute& defConn = connList.requireObject("connection", "*");
	defConn.requireScalar(DBMS) = DEFAULT_DBMS;
	defConn.requireScalar(HOST) = DEFAULT_HOST;
	defConn.requireScalar(PORT).setIntValue(DEFAULT_PORT);
	defConn.requireScalar(USER) = DEFAULT_USER;
	defConn.requireScalar(SID)  = DEFAULT_SID;
	defConn.requireScalar(SCHEMA) = DEFAULT_SCHEMA;

	// update connections in the configuration from the list in the tree.
	for(boost::ptr_vector<ConnectionItem>::const_iterator curConn = connections.begin(),
		endConn = connections.end(); curConn != endConn; ++curConn)
	{
		const std::string& dbms = Io::DatabaseStorage::dbmsToString(curConn->getDbmsType());
		const std::string& host = curConn->getHost();
		const std::string& port = curConn->getPort();
		const std::string& user = curConn->getUser();
		const std::string& sid = curConn->getSID();
		const std::string& schema = curConn->getSchema();

		Attribute& connec = connList.requireObject("connection", curConn->getName());
		if(dbms != DEFAULT_DBMS)	connec.requireScalar(DBMS) = dbms;
		else						connec.discard(DBMS);
		if(host != DEFAULT_HOST)	connec.requireScalar(HOST) = host;
		else						connec.discard(HOST);
		if(port != DEFAULT_PORT)	connec.requireScalar(PORT).setIntValue(port);
		else						connec.discard(PORT);
		if(user != DEFAULT_USER)	connec.requireScalar(USER) = user;
		else						connec.discard(USER);
		if(sid != DEFAULT_SID)		connec.requireScalar(SID) = sid;
		else						connec.discard(SID);
		if(schema != DEFAULT_SCHEMA)connec.requireScalar(SCHEMA) = schema;
		else						connec.discard(SCHEMA);
	}
	// remove connections which are not in the tree, from the configuration.
	for(Attribute::iterator curAttr = connList.begin(), endAttr = connList.end(); curAttr != endAttr; )
	{
		if(curAttr->getName() == "connection" and curAttr->getIndex() != "*")
		{
			ConnectionItem* conn = search(curAttr->getIndex());
			if(conn == 0)
			{
				Attribute::iterator rmAttr = curAttr++;
				connList.remove(rmAttr); // assumes the iterator curAttr stays valid (ie: not vector container)
			}
			else
			{
				++curAttr;
			}
		}
		else
		{
			++curAttr;
		}
	}
	EditorConfiguration::instance().save();
}

}
}
