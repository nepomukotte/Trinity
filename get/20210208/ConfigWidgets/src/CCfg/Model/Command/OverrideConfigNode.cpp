#include "OverrideConfigNode.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include <CCfg/Attribute.h>
#include <CCfg/Path.h>
#include "Commands.h"
#include <QTextStream>
#include <QModelIndex>

namespace CCfg
{
namespace Model
{
namespace Command
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::OverrideConfigNode class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param node an inherited node to override.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
OverrideConfigNode::OverrideConfigNode(ConfigNode* node,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent), nodePath(node->getPath().toString())
{
	setText(generateText());
}

/**
 * Destructor.
 */
OverrideConfigNode::~OverrideConfigNode()
{
}

/**
 * Returns the ID of this command.
 */
int OverrideConfigNode::id() const
{
	return OverrideConfigNodeCommand;
}

/**
 * Applies the change to the model.
 */
void OverrideConfigNode::redo()
{
	try
	{
		ConfigNode* node = accessor->pathToNode(nodePath);
		std::vector<ConfigNode*> createdNodes;
		createdPath.clear();
		succeeded = accessor->override(node, createdNodes);
		for(std::vector<ConfigNode*>::const_iterator iNode = createdNodes.begin(), endNode = createdNodes.end(); iNode != endNode; ++iNode)
		{
			createdPath.push_back((*iNode)->getPath().toString());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void OverrideConfigNode::undo()
{
	try
	{
		if(succeeded)
		{
			for(PathList::const_iterator iPath = createdPath.begin(), endPath = createdPath.end(); iPath != endPath; ++iPath)
			{
				ConfigNode* node = accessor->pathToNode(*iPath);
				accessor->remove(node->getParent(), node->getRowInParent());
			}
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString OverrideConfigNode::generateText()
{
	QString text;
	QTextStream stream(&text, QIODevice::WriteOnly);
	stream << "Override " << nodePath.c_str();

	return text;
}

}
}
}
