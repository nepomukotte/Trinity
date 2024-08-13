#include "DatabaseStorageFactory.h"
#include "Utils/Date/StdTmHelpers.h"
#include "CCfg/Io/Url.h"
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <vector>
#include <algorithm>

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::DatabaseStorageFactory class methods body
////////////////////////////////////////////////////////////////////////////////

DatabaseStorageFactory::DatabaseStorageFactory()
{
}

DatabaseStorageFactory::~DatabaseStorageFactory()
{
}

/*
 * (see inherited doc)
 */
std::auto_ptr<Storage> DatabaseStorageFactory::createStorage(const Url& url) const
{
	std::auto_ptr<Storage> storage;
	const std::string& scheme = url.getScheme();
	DbmsType dbms = DatabaseStorage::dbmsFromString(scheme);

	std::vector< DbmsType > enabledTypes;
	enabledTypes.push_back(MySql);
#if defined(CCFG_USE_ORACLE)
	enabledTypes.push_back(Oracle);
#endif
#if defined(CCFG_USE_POSTGRESQL)
	enabledTypes.push_back(PostgreSql);
#endif
#if defined(CCFG_USE_SQLITE3)
	enabledTypes.push_back(Sqlite3);
#endif
	if (std::find(enabledTypes.begin(), enabledTypes.end(), dbms) != enabledTypes.end())
	{
		const std::string& path = url.getPath();

		enum {SidPart = 0, SchemaPart, ConfigPart};
		std::string pathComponent[3];
		std::string portStr = "";

		if (dbms == Sqlite3)
		{
			// break the url path in its DB path, schema and config components.
			std::string::size_type pos = path.find_last_of('/');
			pathComponent[ConfigPart] = path.substr(pos+1);
			std::string str = path.substr(0, pos);
			pos = str.find_last_of('/');
			pathComponent[SchemaPart] = str.substr(pos+1);
			pathComponent[SidPart] = str.substr(0, pos);
		}
		else
		{
			// break the url path in its sid, schema and config components.
			typedef boost::char_separator<char> Separator;
			typedef boost::tokenizer<Separator> Tokenizer;

			Separator sep("/");
			Tokenizer tok(path, sep);
			int i = (dbms == Oracle or dbms == PostgreSql ? 0 : 1);
			Tokenizer::iterator cur = tok.begin(), end = tok.end();
			for(; cur != end and i < 3; ++cur, ++i)
			{
				pathComponent[i] = *cur;
			}
			if (dbms == PostgreSql and pathComponent[ConfigPart].empty())
			{
				pathComponent[ConfigPart] = pathComponent[SchemaPart];
				pathComponent[SchemaPart] = "";
			}
			portStr = boost::lexical_cast<std::string>(url.getPort());
		}
		// retrieve the requested version date.
		std::string versionDateStr = url.getQueryItem(DatabaseStorage::VERSION_DATE_KEY);

		// build the storage
		DatabaseStorage* dbStorage = new DatabaseStorage(dbms,
				url.getHost(),
				portStr,
				url.getUser(),
				url.getPassword(),
				pathComponent[SidPart],
				pathComponent[SchemaPart]
		);
		storage.reset(dbStorage);

		// parameterize the storage
		dbStorage->setLoginCallback(loginCallback);
		dbStorage->setCurrentConfigName(pathComponent[ConfigPart]);
		if(not versionDateStr.empty())
		{
			// TODO: bad date error handling ?
			dbStorage->setCurrentConfigDate(boost::lexical_cast<std::tm>(versionDateStr));
		}
	}
	return storage;
}

}
}
