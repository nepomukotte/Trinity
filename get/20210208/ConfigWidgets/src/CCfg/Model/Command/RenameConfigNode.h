#ifndef CCFG_MODEL_COMMAND_RENAMECONFIGNODE_H_INCLUDED
#define CCFG_MODEL_COMMAND_RENAMECONFIGNODE_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "CCfg/Path.h"
#include <string>

namespace CCfg
{
namespace Model
{
// forward declarations
class ConfigNode;

namespace Command
{
///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This command changes the name and index of a configuration node.
 *
 * @author fchateau
 */
class RenameConfigNode : public AbstractConfigModelCommand
{
	Path nodePath;	///< the path of the node to modify.
	std::string oldValue;	///< the value before command execution.
	std::string newValue;	///< the value that the command will set in redo().
public:
	RenameConfigNode(ConfigNode* node, const std::string& newValue,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent = 0);
	~RenameConfigNode();

	int id() const;
	bool mergeWith(const QUndoCommand* command);

	void redo();
	void undo();

private:
	QString generateText();

	// non-default-constructible, non-copyable
	RenameConfigNode();
	RenameConfigNode(const RenameConfigNode& cmd);
	RenameConfigNode& operator=(const RenameConfigNode& cmd);
};

}

}

}

#endif // CCFG_MODEL_COMMAND_RENAMECONFIGNODE_H_INCLUDED
