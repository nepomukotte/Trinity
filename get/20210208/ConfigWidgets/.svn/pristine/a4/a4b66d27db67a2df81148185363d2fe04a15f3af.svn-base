#include "IncludeConfig.h"
#include "DeleteConfigNode.h"
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
// CCfg::Model::Command::IncludeConfig class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
IncludeConfig::IncludeConfig(ConfigNode* destination, const std::string& url,
	ConfigModelAccessor* accessor,
	AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	deleteCmd(new DeleteConfigNode(destination,
			listIndexesToRemove(destination),
			accessor,
			this)),
	setUrlCmd(new SetConfigNodeIncludeUrl(destination, url, accessor, this)),
	destPath(destination->getPath()), url(url)
{
	setText(generateText());
}

/**
 * Destructor.
 */
IncludeConfig::~IncludeConfig()
{
}

/**
 * Returns the ID of this command.
 */
int IncludeConfig::id() const
{
	return IncludeConfigCommand;
}

/**
 * Applies the change to the model.
 */
void IncludeConfig::redo()
{

	try
	{
		AbstractConfigModelCommand::redo();	// calls redo on nested commands.
		if(deleteCmd->hasSucceeded())
		{
			ConfigNode* node = accessor->pathToNode(destPath);
			succeeded = accessor->include(node, url);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void IncludeConfig::undo()
{
	try
	{
		// only undo if redo succeeded
		if(succeeded)
		{
			ConfigNode* node = accessor->pathToNode(destPath);
			accessor->remove(node, listIndexesToRemove(node));
			AbstractConfigModelCommand::undo(); // calls undo on nested commands.
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString IncludeConfig::generateText()
{
	std::ostringstream text;
	text << "Included " << url << " into " << destPath;
	return text.str().c_str();
}

/**
 * Makes the list of the indexes of the nodes to remove before importing the
 * configuration.
 * In fact we need to remove all the nodes.
 */
IndexList IncludeConfig::listIndexesToRemove(ConfigNode* node)
{
	IndexList indexes(node->getChildren().size());
	int i = 0;
	for(ConfigNode::const_iterator child = node->begin(); child != node->end(); ++child, ++i)
	{
		indexes[i] = i;
	}
	return indexes;
}

}
}
}
