#include "ImportConfig.h"
#include "DeleteConfigNode.h"
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
// CCfg::Model::Command::ImportConfig class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
ImportConfig::ImportConfig(ConfigNode* destination, const std::string& url,
	ConfigModelAccessor* accessor,
	AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	deleteCmd(new DeleteConfigNode(destination,
			listIndexesToRemove(destination),
			accessor,
			this)),
	destPath(destination->getPath()), url(url)
{
	setText(generateText());
}

/**
 * Destructor.
 */
ImportConfig::~ImportConfig()
{
}

/**
 * Returns the ID of this command.
 */
int ImportConfig::id() const
{
	return ImportConfigCommand;
}

/**
 * Applies the change to the model.
 */
void ImportConfig::redo()
{

	try
	{
		AbstractConfigModelCommand::redo();	// calls redo on nested commands.
		if(deleteCmd->hasSucceeded())
		{
			ConfigNode* node = accessor->pathToNode(destPath);
			succeeded = accessor->import(node, url);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void ImportConfig::undo()
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
QString ImportConfig::generateText()
{
	std::ostringstream text;
	text << "Imported " << url << " into " << destPath;
	return text.str().c_str();
}

/**
 * Makes the list of the indexes of the nodes to remove before importing the
 * configuration.
 * In fact we need to remove all the nodes.
 */
IndexList ImportConfig::listIndexesToRemove(ConfigNode* node)
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
