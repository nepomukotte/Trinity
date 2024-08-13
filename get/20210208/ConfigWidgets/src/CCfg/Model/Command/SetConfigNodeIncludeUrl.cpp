#include "SetConfigNodeIncludeUrl.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "Commands.h"
#include <sstream>

namespace CCfg
{
namespace Model
{
namespace Command
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::SetConfigNodeIncludeUrl class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
SetConfigNodeIncludeUrl::SetConfigNodeIncludeUrl(ConfigNode* node,
		const std::string& newUrl,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	nodePath(node->getPath()),
	oldUrl(node->getIncludeUrl()),
	newUrl(newUrl)
{
	setText(generateText());
}

/**
 * Destructor.
 */
SetConfigNodeIncludeUrl::~SetConfigNodeIncludeUrl()
{
}

/**
 * Returns the ID of this command.
 */
int SetConfigNodeIncludeUrl::id() const
{
	return SetConfigNodeIncludeUrlCommand;
}

/**
 * Attempts to merge this command with the specified command.
 * @see bool RenameConfigNode::mergeWith(const QUndoCommand* command)
 */
bool SetConfigNodeIncludeUrl::mergeWith(const QUndoCommand* command)
{
	if(command->id() == id())
	{
		const SetConfigNodeIncludeUrl* cmd = static_cast<const SetConfigNodeIncludeUrl*>(command);
		if(succeeded and cmd->succeeded and // we only merge successful commands;
			nodePath == cmd->nodePath and	// with the same path;
			newUrl == cmd->oldUrl)	// and matching transitivity rule.
		{
			newUrl = cmd->newUrl;	// apply transitivity.
			setText(generateText());	// update text.
			return true;
		}
	}
	return false;
}

/**
 * Applies the change to the model.
 */
void SetConfigNodeIncludeUrl::redo()
{
	try
	{
		ConfigNode* node = accessor->pathToNode(nodePath);
		succeeded = accessor->setIncludeUrl(node, newUrl);
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void SetConfigNodeIncludeUrl::undo()
{
	try
	{
		if(succeeded)
		{
			ConfigNode* node = accessor->pathToNode(nodePath);
			accessor->setIncludeUrl(node, oldUrl);
		}
	}
	CATCH_AND_REPORT
}

/*
 *
 */
QString SetConfigNodeIncludeUrl::generateText()
{
	std::ostringstream stream;
	stream << "URL of " << nodePath << " changed from \""
		<< oldUrl << "\" to \"" << newUrl << '"';
	return QString::fromStdString(stream.str());
}


}
}
}
