#include "InsertConfigNode.h"
#include "../ConfigModelAccessor.h"
#include "../ConfigNode.h"
#include "Utils/Error/ExcepHandling.h"
#include "CCfg/Attribute.h"
#include "CCfg/PathNode.h"
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
// CCfg::Model::Command::InsertConfigNode class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param referenceNode an existing node that serves as a landmark for insertion.
 * @param position the position of insertion with respect to the reference node.
 * @param attribute the attribute to insert. The command takes its ownership.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
InsertConfigNode::InsertConfigNode(ConfigNode* referenceNode,
		Position position,
		std::auto_ptr<Attribute> attribute,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent)
	: AbstractConfigModelCommand(accessor, parent),
	referencePath(referenceNode->getPath()),
	position(position),
	newAttribute(attribute)
{
	setText(generateText());
}

/**
 * Destructor.
 */
InsertConfigNode::~InsertConfigNode()
{
}

/**
 * Returns the ID of this command.
 */
int InsertConfigNode::id() const
{
	return InsertConfigNodeCommand;
}

/**
 * Applies the change to the model.
 */
void InsertConfigNode::redo()
{
	try
	{
		ConfigNode* reference = accessor->pathToNode(referencePath);

		if(position == Child)
		{
			succeeded = accessor->insertInside(reference, newAttribute->clone());
		}
		else
		{
			succeeded = accessor->insertAfter(reference, newAttribute->clone());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Reverts the change to the model.
 */
void InsertConfigNode::undo()
{
	try
	{
		if(succeeded)
		{
			ConfigNode* reference = accessor->pathToNode(referencePath);
			PathNode nodeName(*newAttribute);

			ConfigNode* parent = (position == Child) ? reference : reference->getParent();
			ConfigNode* child = parent->search(nodeName);
			assert(child != 0);	// this should not happen since succeeded flag is true.
			accessor->remove(parent, child->getRowInParent());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Generates a title for the command.
 * @return Returns the generated title for the command.
 */
QString InsertConfigNode::generateText()
{
	QString text;
	QTextStream stream(&text, QIODevice::WriteOnly);
	PathNode node(*newAttribute);
	stream << "Insert " << node.toString().c_str();

	return text;
}

}
}
}
