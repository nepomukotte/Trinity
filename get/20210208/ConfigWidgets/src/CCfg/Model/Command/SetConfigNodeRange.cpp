#include "SetConfigNodeRange.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "Commands.h"
#include "CCfg/Range.h"
#include "Utils/Io/PointerValueStreamer.hpp"
#include <QModelIndex>
#include <sstream>

using Utils::Io::writeValue;

namespace CCfg
{
namespace Model
{
namespace Command
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::SetConfigNodeRange class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
SetConfigNodeRange::SetConfigNodeRange(ConfigNode* node,
	RangePtr newRange,
	ConfigModelAccessor* accessor,
	AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	nodePath(node->getPath()),
	oldRange(),
	newRange(newRange)
{
	Range* r = node->getRange();
	if(r != 0)
	{
		oldRange = r->clone();
	}
	setText(generateText());
}

/**
 * Destructor.
 */
SetConfigNodeRange::~SetConfigNodeRange()
{
}

/**
 * Returns the ID of this command.
 */
int SetConfigNodeRange::id() const
{
	return SetConfigNodeRangeCommand;
}

/**
 * Attempts to merge this command with the specified command.
 * @see bool RenameConfigNode::mergeWith(const QUndoCommand* command)
 */
bool SetConfigNodeRange::mergeWith(const QUndoCommand* command)
{
	if(command->id() == id())
	{
		const SetConfigNodeRange* cmd = static_cast<const SetConfigNodeRange*>(command);
		if(succeeded and cmd->succeeded and // we only merge successful commands;
			nodePath == cmd->nodePath and	// with the same path;
			Range::equivalent(newRange.get(), cmd->oldRange.get()))	// and matching transitivity rule.
		{
			// apply transitivity.
			if(cmd->newRange.get() != 0)
			{
				newRange = cmd->newRange->clone();
			}
			else
			{
				newRange.reset();
			}
			setText(generateText());	// update text.
			return true;
		}
	}
	return false;
}

/**
 * Applies the change to the model.
 */
void SetConfigNodeRange::redo()
{
	try
	{
		ConfigNode* node = accessor->pathToNode(nodePath);
		succeeded = accessor->setRange(node, newRange.get() != 0 ? newRange->clone() : newRange);
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void SetConfigNodeRange::undo()
{
	try
	{
		// only undo if redo succeeded
		if(succeeded)
		{
			ConfigNode* node = accessor->pathToNode(nodePath);
			accessor->setRange(node, oldRange.get() != 0 ? oldRange->clone() : oldRange);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Returns a generated title for the command.
 */
QString SetConfigNodeRange::generateText()
{
	std::ostringstream stream;
	stream << "Range of " << nodePath
		<< " changed from " << writeValue(oldRange.get())
		<< " to " << writeValue(newRange.get());

	return QString::fromStdString(stream.str());
}

}
}
}
