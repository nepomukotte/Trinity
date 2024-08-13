#ifndef CCFG_MODEL_CONFIGMODEL2_H_INCLUDED
#define CCFG_MODEL_CONFIGMODEL2_H_INCLUDED

#include "CCfg/Model/ConfigModel.h"

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is extension to the configuration model adding methods and attributes to
 * make it work with databases.
 *
 * @author fchateau
 */
class ConfigModel2 : public ConfigModel
{
	std::string dbUrl;

public:
	ConfigModel2(QObject* parent = 0);
	virtual ~ConfigModel2();

	void commit();
	void commitInto(const Io::Url& url);

	const std::string& getUrlInDatabase() const { return dbUrl; }
	void setUrlInDatabase(const std::string& url) { dbUrl = url; }

	void setDocumentUrl(const std::string& url);

private:
	ConfigModel2(const ConfigModel2& rvalue);
	ConfigModel2& operator=(const ConfigModel2& rvalue);
};

}
}

#endif // CCFG_MODEL_CONFIGMODEL2_H_INCLUDED
