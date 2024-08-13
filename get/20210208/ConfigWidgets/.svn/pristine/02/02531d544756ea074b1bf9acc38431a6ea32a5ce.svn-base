#include "DeleteConfigNode.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "CCfg/Object.h"
#include "CCfg/Attribute.h"
#include "Commands.h"
#include <QTextStream>

namespace CCfg
{
namespace Model
{
namespace Command
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::DeleteConfigNode class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param parentNode the parent of the ConfigNodes to remove.
 * @param indexes the sorted list of row indexes of the ConfigNodes to remove.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
DeleteConfigNode::DeleteConfigNode(ConfigNode* parentNode,
		const IndexList& indexes,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	parentPath(parentNode->getPath()), indexes(indexes),
	deletedContainer(new Object)
{
	setText(generateText());
}

/**
 * Destructor.
 */
DeleteConfigNode::~DeleteConfigNode()
{
}

/**
 *
 */
int DeleteConfigNode::id() const
{
	return DeleteConfigNodeCommand;
}

/**
 *
 */
void DeleteConfigNode::redo()
{
	try
	{
		ConfigNode* parent = accessor->pathToNode(parentPath);
		succeeded = accessor->remove(parent, indexes, deletedContainer.get());
	}
	CATCH_AND_REPORT
}

/**
 *
 */
void DeleteConfigNode::undo()
{
	try
	{
		ConfigNode* parent = accessor->pathToNode(parentPath);
		if(succeeded)
		{
			accessor->insert(parent, deletedContainer.get(), indexes);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString DeleteConfigNode::generateText()
{
	QString text;
	QTextStream stream(&text, QIODevice::WriteOnly);
	stream << "Remove nodes in " << parentPath.toString().c_str();

	return text;
}

}
}
}
