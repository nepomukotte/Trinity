#include "SetConfigNodeUnit.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "Commands.h"
#include <QTextStream>
#include <QModelIndex>

namespace CCfg
{
namespace Model
{
namespace Command
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::SetConfigNodeUnit class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param node the node whose unit must be changed.
 * @param newUnit the new unit of the node.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
SetConfigNodeUnit::SetConfigNodeUnit(ConfigNode* node,
		const std::string& newUnit,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	nodePath(node->getPath()),
	oldUnit(node->getUnit()),
	newUnit(newUnit)
{
	setText(generateText());
}

/**
 * Destructor.
 */
SetConfigNodeUnit::~SetConfigNodeUnit()
{
}

/**
 * Returns the ID of this command.
 * This command supports compression, so its ID is not -1.
 * @return Returns the ID of this command.
 */
int SetConfigNodeUnit::id() const
{
	return SetConfigNodeUnitCommand;
}

/**
 * Attempts to merge this command with the specified command.
 * @see bool RenameConfigNode::mergeWith(const QUndoCommand* command)
 */
bool SetConfigNodeUnit::mergeWith(const QUndoCommand* command)
{
	if(command->id() == id())
	{
		const SetConfigNodeUnit* cmd = static_cast<const SetConfigNodeUnit*>(command);
		if(succeeded and cmd->succeeded and // we only merge successful commands;
			nodePath == cmd->nodePath and	// with the same path;
			newUnit == cmd->oldUnit)	// and matching transitivity rule.
		{
			newUnit = cmd->newUnit;	// apply transitivity.
			setText(generateText());	// update text.
			return true;
		}
	}
	return false;
}

/**
 * Applies the change to the model.
 */
void SetConfigNodeUnit::redo()
{
	try
	{
		ConfigNode* node = accessor->pathToNode(nodePath);
		succeeded = accessor->setUnit(node, newUnit);
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void SetConfigNodeUnit::undo()
{
	try
	{
		// only undo if redo succeeded
		if(succeeded)
		{
			ConfigNode* node = accessor->pathToNode(nodePath);
			accessor->setUnit(node, oldUnit);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString SetConfigNodeUnit::generateText()
{
	QString text;
	QTextStream stream(&text, QIODevice::WriteOnly);
	stream << "Unit changed from " << oldUnit.c_str() << " to " << newUnit.c_str();

	return text;
}

}
}
}
