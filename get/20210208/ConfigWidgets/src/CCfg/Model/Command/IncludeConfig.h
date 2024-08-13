#ifndef CCFG_MODEL_COMMAND_INCLUDECONFIG_H_INCLUDED
#define CCFG_MODEL_COMMAND_INCLUDECONFIG_H_INCLUDED

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
 * This command includes a configuration resource into a config node.
 *
 * @author fchateau
 */
class IncludeConfig : public AbstractConfigModelCommand
{
	AbstractConfigModelCommand* deleteCmd;
	AbstractConfigModelCommand* setUrlCmd;
	Path destPath;
	std::string url;

public:
	IncludeConfig(ConfigNode* destination, const std::string& url,
			ConfigModelAccessor* accessor,
			AbstractConfigModelCommand* parent = 0);
	virtual ~IncludeConfig();

	int id() const;

	void redo();
	void undo();

private:
	QString generateText();
	static IndexList listIndexesToRemove(ConfigNode* node);

	// non-default-constructible, non-copyable
	IncludeConfig();
	IncludeConfig(const IncludeConfig& rvalue);
	IncludeConfig& operator=(const IncludeConfig& rvalue);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_INCLUDECONFIG_H_INCLUDED
