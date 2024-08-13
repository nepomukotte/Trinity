#include "DatabaseFormatFactory.h"
#include "DatabaseFormat.h"
#include "DatabaseStorage.h"
#include "CCfg/Io/Url.h"

namespace CCfg
{
namespace Io
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::DatabaseFormatFactory class methods body
////////////////////////////////////////////////////////////////////////////////

DatabaseFormatFactory::DatabaseFormatFactory()
{
}

DatabaseFormatFactory::~DatabaseFormatFactory()
{
}

/*
 * (see inherited doc)
 */
std::auto_ptr<Format> DatabaseFormatFactory::createFormat(const Url& url) const
{
	std::auto_ptr<Format> format;
	const std::string& scheme = url.getScheme();
	if(scheme == DatabaseStorage::SCHEME[Oracle] or
		scheme == DatabaseStorage::SCHEME[MySql] or
		scheme == DatabaseStorage::SCHEME[PostgreSql] or
		scheme == DatabaseStorage::SCHEME[Sqlite3])
	{
		DatabaseFormat* dbFormat = new DatabaseFormat;
		format.reset(dbFormat);
		dbFormat->setUpdateController(updateController);
	}
	return format;
}

}
}
