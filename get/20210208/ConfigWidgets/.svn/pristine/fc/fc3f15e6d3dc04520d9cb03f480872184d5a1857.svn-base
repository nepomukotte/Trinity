#include "SetConfigNodeValue.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "Commands.h"
#include "CCfg/Value.h"
#include <QTextStream>
#include <QModelIndex>

namespace CCfg
{
namespace Model
{
namespace Command
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::SetConfigNodeValue class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param node the node whose value must be changed.
 * @param newValue the new value of the node.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
SetConfigNodeValue::SetConfigNodeValue(ConfigNode* node,
		const Value& newValue,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	nodePath(node->getPath()),
	oldValue(node->getValue()),
	newValue(newValue)
{
	setText(generateText());
}

/**
 * Destructor.
 */
SetConfigNodeValue::~SetConfigNodeValue()
{
}

/**
 * Returns the ID of this command.
 * This command supports compression, so its ID is not -1.
 * @return Returns the ID of this command.
 */
int SetConfigNodeValue::id() const
{
	return SetConfigNodeValueCommand;
}

/**
 * Attempts to merge this command with the specified command.
 * @see bool RenameConfigNode::mergeWith(const QUndoCommand* command)
 */
bool SetConfigNodeValue::mergeWith(const QUndoCommand* command)
{
	if(command->id() == id())
	{
		const SetConfigNodeValue* cmd = static_cast<const SetConfigNodeValue*>(command);
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
void SetConfigNodeValue::redo()
{
	try
	{
		ConfigNode* node = accessor->pathToNode(nodePath);
		succeeded = accessor->setValue(node, newValue);
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void SetConfigNodeValue::undo()
{
	try
	{
		// only undo if redo succeeded
		if(succeeded)
		{
			ConfigNode* node = accessor->pathToNode(nodePath);
			accessor->setValue(node, oldValue);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString SetConfigNodeValue::generateText()
{
	QString text;
	QTextStream stream(&text, QIODevice::WriteOnly);
	stream << "Value changed from " << oldValue.getAsString().c_str() << " to " << newValue.getAsString().c_str();

	return text;
}


}
}
}
