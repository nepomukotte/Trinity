#include "DatabaseFormat.h"
#include "DatabaseStorage.h"
#include "Utils/Date/StdTmHelpers.h"
#include "Db/Config.h"
#include "Db/ConfigVersion.h"
#include "Db/ConfigDao.h"
#include "CCfg/Object.h"
#include <sstream>
#include <ctime>

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::DatabaseFormat class methods body
////////////////////////////////////////////////////////////////////////////////

using namespace Db;

/**
 * Default constructor.
 */
DatabaseFormat::DatabaseFormat()
{
}

/**
 * Destructor.
 */
DatabaseFormat::~DatabaseFormat()
{
}

/*
 * (see inherited doc)
 */
std::auto_ptr<CCfg::Object> DatabaseFormat::parse(Storage& storage, IncludeHandler& /*handler*/)
{
	std::auto_ptr<Object> cfgRoot;
	DatabaseStorage& database = dynamic_cast<DatabaseStorage&>(storage);
	if(not database.isOpened())
	{
		database.open();
		if(not database.isOpened())
		{
			// case where connection is canceled through the login callback
			return cfgRoot;
		}
	}
	ConfigDao& dao = database.getDao();
	const std::string& configName = database.getCurrentConfigName();
	const std::tm& configDate = database.getCurrentConfigDate();

	if(configName.empty())
	{
		throw MappingError(Extraction, "configuration name should be specified");
	}

	ConfigPtr config = dao.searchConfig(configName);
	if(config.get() == 0)
	{
		throw MappingError(Extraction, (std::string("cannot find the configuration \"") += configName) += '\"');
	}

	ConfigVersionPtr version;
	if(configDate.tm_year == 0)
	{
		// assume that the date is not set => retrieve the latest version.
		version = config->getLastVersion();
	}
	else
	{
		version = config->getVersionAt(configDate);
	}

	if(version.get() == 0)
	{
		std::ostringstream errMsg;
		errMsg << "there is no version at: " << configDate;
		throw MappingError(Extraction, errMsg.str());
	}
	cfgRoot.reset(static_cast<Object*>(version->extractConfig().release()));
	return cfgRoot;
}

/*
 * (see inherited doc)
 */
void DatabaseFormat::format(const Object *cfgRoot, Storage& storage, IncludeHandler& /*handler*/)
{
	DatabaseStorage& database = dynamic_cast<DatabaseStorage&>(storage);
	if(not database.isOpened())
	{
		database.open();
		if(not database.isOpened())
		{
			// case where connection is canceled through the login callback
			return;
		}
	}
	ConfigDao& dao = database.getDao();

	const std::string& configName = database.getCurrentConfigName();
	if(configName.empty())
	{
		throw MappingError(Insertion, "configuration name should be specified");
	}
	ConfigPtr config = dao.searchConfig(configName);
	if(config.get() == 0)
	{
		// if the configuration does not exist, we create it.
		config = Config::createConfig(configName, &dao);
		dao.insert(*config);
	}
	if(updateController == 0)
	{
		updateController.reset(new Db::DefaultUpdateController);
	}
	config->insertNewVersion(cfgRoot, *updateController);
}

/*
 * (see inherited doc)
 */
std::string DatabaseFormat::getConfigName(Storage& storage)
{
	DatabaseStorage& database = dynamic_cast<DatabaseStorage&>(storage);
	return database.getCurrentConfigName();
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::DatabaseFormat::MappingError class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
DatabaseFormat::MappingError::MappingError(Operation op, const std::string& message)
	: op(op), msg(message)
{
}

/**
 * Copy constructor.
 */
DatabaseFormat::MappingError::MappingError(const MappingError& rvalue)
	: Exception(rvalue), op(rvalue.op), msg(rvalue.msg)
{
}

/**
 * Assignment operator.
 */
DatabaseFormat::MappingError& DatabaseFormat::MappingError::operator=(const MappingError& rvalue)
{
	op = rvalue.op;
	msg = rvalue.msg;
	Exception::operator=(rvalue);
	return *this;
}

/*
 * (see inherited doc)
 */
std::string DatabaseFormat::MappingError::formatMessage() const
{
	std::ostringstream errorMsg;
	errorMsg <<
		"Configuration " <<
		((op == Insertion) ? "insertion" : "extraction") <<
		" failed: " << msg;
	return errorMsg.str();
}

}
}
