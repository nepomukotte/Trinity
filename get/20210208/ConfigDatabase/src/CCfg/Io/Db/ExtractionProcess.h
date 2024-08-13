#ifndef CCFG_IO_DB_EXTRACTIONPROCESS_H_INCLUDED
#define CCFG_IO_DB_EXTRACTIONPROCESS_H_INCLUDED

#include "ConfigVersion.h"
#include "ConfigAttribute.h"
#include "ConfigClassInfo.h"
#include "ConfigClass.h"
#include "ConfigVersion.h"
#include "CCfg/Attribute.h"

namespace CCfg
{
namespace Io
{
namespace Db
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Io::Db::ExtractionProcess class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * This class encapsulates the process of translating a database configuration
 * into a compound config.
 *
 * @author fchateau
 */
class ExtractionProcess
{
	ConfigVersion& version;
	std::auto_ptr<Attribute> result;

	ConfigClassPtrSetPtr classes;
	ConfigClassInfoPtrSetPtr classInfos;
	ConfigAttributePtrSetPtr attributes;

	ConfigClassIdIdx* clsIdIdx;
	ConfigClassInfoIdIdx* infoIdIdx;
	ConfigAttributeIdIdx* attrIdIdx;

public:
	class Error;

	ExtractionProcess(ConfigVersion& version);
	virtual ~ExtractionProcess();

	std::auto_ptr<Attribute>& extract();

private:
	void loadData();
	ConfigClass& getClass(UID id);
	ConfigClassInfo& getClassInfo(UID id);
	ConfigAttribute& getAttribute(UID id);
};

/**
 * Reports an error while extracting a configuration from a database.
 */
class ExtractionProcess::Error : public Exception
{
public:
	Error();
	~Error() throw() {}
	Error(const Error& rvalue);
	Error& operator=(const Error& rvalue);

private:
	std::string formatMessage() const;
};

}
}
}

#endif // CCFG_IO_DB_EXTRACTIONPROCESS_H_INCLUDED
