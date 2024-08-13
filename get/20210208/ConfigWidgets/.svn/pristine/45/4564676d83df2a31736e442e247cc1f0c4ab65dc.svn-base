#ifndef CCFG_MODEL_EDITORCONFIGURATION_H_INCLUDED
#define CCFG_MODEL_EDITORCONFIGURATION_H_INCLUDED

#include "CCfg/Io/Document.h"
#include <memory>

namespace CCfg
{
namespace Model
{

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This singleton centralizes the access to the configuration of the program.
 *
 * @author fchateau
 */
class EditorConfiguration
{
	static std::auto_ptr<EditorConfiguration> programConfig;
	static const char CONFIG_FILE_NAME[];
	Io::Document doc;

public:
	~EditorConfiguration();

	static EditorConfiguration& instance();
	static Attribute& root();

	Attribute& getConfig();
	void save();

private:
	EditorConfiguration();
	EditorConfiguration(const EditorConfiguration& rvalue);
	EditorConfiguration& operator=(const EditorConfiguration& rvalue);
};

}
}

#endif // CCFG_MODEL_EDITORCONFIGURATION_H_INCLUDED
