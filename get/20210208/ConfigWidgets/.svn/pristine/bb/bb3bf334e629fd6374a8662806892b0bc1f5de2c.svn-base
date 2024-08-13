#ifndef CCFG_MODEL_COMMAND_IMPORTCONFIG_H_INCLUDED
#define CCFG_MODEL_COMMAND_IMPORTCONFIG_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "../AbstractConfigModelAccessor.h"
#include "CCfg/Path.h"

namespace CCfg
{
namespace Model
{
// forward declarations
class ConfigNode;

namespace Command
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This command imports a configuration resource into a config node.
 *
 * @author fchateau
 */
class ImportConfig : public AbstractConfigModelCommand
{
	AbstractConfigModelCommand* deleteCmd;
	Path destPath;
	std::string url;

public:
	ImportConfig(ConfigNode* destination, const std::string& url,
			ConfigModelAccessor* accessor,
			AbstractConfigModelCommand* parent = 0);
	virtual ~ImportConfig();

	int id() const;

	void redo();
	void undo();

private:
	QString generateText();
	static IndexList listIndexesToRemove(ConfigNode* node);

	// non-default-constructible, non-copyable
	ImportConfig();
	ImportConfig(const ImportConfig& rvalue);
	ImportConfig& operator=(const ImportConfig& rvalue);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_IMPORTCONFIG_H_INCLUDED
