#ifndef CCFG_MODEL_COMMAND_OVERRIDECONFIGNODE_H_INCLUDED
#define CCFG_MODEL_COMMAND_OVERRIDECONFIGNODE_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include <vector>
#include <string>

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
 * This command overrides an inherited configuration node.
 *
 * @author fchateau
 */
class OverrideConfigNode : public AbstractConfigModelCommand
{
	typedef std::vector<std::string> PathList;
	const std::string nodePath; ///< the path of the node to override.
	PathList createdPath; ///< the path of the beginning of the created overridden branches

public:
	OverrideConfigNode(ConfigNode* node,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent = 0);
	~OverrideConfigNode();

	int id() const;

	void redo();
	void undo();

protected:
	QString generateText();

private:
	// non-default-constructible, non-copyable
	OverrideConfigNode();
	OverrideConfigNode(const OverrideConfigNode& right);
	OverrideConfigNode& operator=(const OverrideConfigNode& right);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_OVERRIDECONFIGNODE_H_INCLUDED
