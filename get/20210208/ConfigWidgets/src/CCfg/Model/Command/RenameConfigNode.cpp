#include "RenameConfigNode.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "Commands.h"
#include <QTextStream>
#include <QModelIndex>
#include <cassert>

namespace CCfg
{
namespace Model
{
namespace Command
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::RenameConfigNode class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param node the node to rename.
 * @param newValue the new name and index for the node.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
RenameConfigNode::RenameConfigNode(ConfigNode* node,
		const std::string& newValue,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	nodePath(node->getPath()),
	oldValue(node->getFullName()),
	newValue(newValue)
{
	setText(generateText());
}

RenameConfigNode::~RenameConfigNode()
{
}

/**
 * Returns the ID of this command.
 * This command supports compression, so its ID is not -1.
 * @return Returns the ID of this command.
 */
int RenameConfigNode::id() const
{
	return RenameConfigNodeCommand;
}

/**
 * Attempts to merge this command with the specified command.
 * This command supports compression, so under certain conditions it may be
 * merged with another command, applying the transitivity rule.
 * If command C1 on attribute A changes x -> y
 * and command C2 on attribute A changes y -> z
 * Then the merge C3 = M(C1, C2) gives a command  C3 on attribute A
 * changing x -> z.
 * @note Qt only merges consecutive commands.
 * @return Returns true on success; otherwise returns false.
 */
bool RenameConfigNode::mergeWith(const QUndoCommand* command)
{
	if(command->id() == id())
	{
		const RenameConfigNode* cmd = static_cast<const RenameConfigNode*>(command);
		if(succeeded and cmd->succeeded and // we only merge successful commands;
			nodePath == cmd->nodePath and	// with the same path;
			newValue == cmd->oldValue)	// and matching transitivity rule.
		{
			newValue = cmd->newValue;	// apply transitivity.
			setText(generateText());	// update text.
			return true;
		}
	}
	return false;
}

/**
 * Applies the change to the model.
 */
void RenameConfigNode::redo()
{
	try
	{
		ConfigNode* node = accessor->pathToNode(nodePath);
		assert(node != 0);
		succeeded = accessor->setNameAndIndex(node, newValue);
		// update the path (renaming the node changed it).
		nodePath.getTail().fromString(newValue);
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void RenameConfigNode::undo()
{
	try
	{
		// only undo if redo succeeded
		if(succeeded)
		{
			ConfigNode* node = accessor->pathToNode(nodePath);
			accessor->setNameAndIndex(node, oldValue);
			// update the path (renaming the node changed it).
			nodePath.getTail().fromString(oldValue);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString RenameConfigNode::generateText()
{
	QString text;
	QTextStream stream(&text, QIODevice::WriteOnly);
	stream << "Rename " << oldValue.c_str() << " to " << newValue.c_str();

	return text;
}

}
}
}
