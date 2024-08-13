#include "EditorConfiguration.h"
#include "CCfg/Object.h"
#include <QFile>

namespace CCfg
{
namespace Model
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::EditorConfiguration class methods body
////////////////////////////////////////////////////////////////////////////////

std::auto_ptr<EditorConfiguration> EditorConfiguration::programConfig;
const char EditorConfiguration::CONFIG_FILE_NAME[] = "settings.xml";

/**
 * Default constructor.
 */
EditorConfiguration::EditorConfiguration() : doc()
{
	if(not QFile::exists(CONFIG_FILE_NAME))
	{
		doc.setConfig(std::auto_ptr<Object>(new Object("editor")));
		doc.saveToFile(CONFIG_FILE_NAME);
	}
	else
	{
		doc.loadFromFile(CONFIG_FILE_NAME);
	}
}

/**
 * Destructor.
 */
EditorConfiguration::~EditorConfiguration()
{
}

/**
 * Retrieves the singleton instance.
 */
EditorConfiguration& EditorConfiguration::instance()
{
	if(programConfig.get() == 0)
	{
		programConfig.reset(new EditorConfiguration);
	}
	return *programConfig;
}

/**
 * Directly retrieve the configuration root of the singleton instance.
 */
Attribute& EditorConfiguration::root()
{
	return instance().getConfig();
}

/**
 * Returns the configuration root.
 */
Attribute& EditorConfiguration::getConfig()
{
	return *doc.getConfig();
}

/**
 * Saves the configuration.
 */
void EditorConfiguration::save()
{
	doc.save();
}

}
}
