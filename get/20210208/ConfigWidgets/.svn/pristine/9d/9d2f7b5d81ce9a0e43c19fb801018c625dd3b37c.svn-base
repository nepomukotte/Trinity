#include "UndoableConfigModelAccessor.h"
#include "ConfigModelAccessor.h"
#include "ConfigNode.h"
#include "Command/InsertConfigNode.h"
#include "Command/OverrideConfigNode.h"
#include "Command/DeleteConfigNode.h"
#include "Command/ImportConfig.h"
#include "Command/IncludeConfig.h"
#include "Command/RenameConfigNode.h"
#include "Command/SetConfigNodeValue.h"
#include "Command/SetConfigNodeUnit.h"
#include "Command/SetConfigNodeRange.h"
#include "Command/SetConfigNodeIncludeUrl.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Attribute.h"
#include "CCfg/Range.h"
#include <QUndoStack>

namespace CCfg
{
namespace Model
{
using namespace Command;

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::CommandBasedConfigModelAccessor class methods body
///////////////////////////////////////////////////////////////////////////////

CommandBasedConfigModelAccessor::CommandBasedConfigModelAccessor(ConfigModelAccessor* accessor)
	: accessor(accessor)
{
}

CommandBasedConfigModelAccessor::~CommandBasedConfigModelAccessor()
{
}

/*
 *
 */
void CommandBasedConfigModelAccessor::setDocument(std::auto_ptr<Io::Document> newDocument)
{
	accessor->setDocument(newDocument);
}

/*
 *
 */
const Io::Document* CommandBasedConfigModelAccessor::getDocument() const
{
	return accessor->getDocument();
}

/*
 *
 */
Io::Document* CommandBasedConfigModelAccessor::getDocument()
{
	return accessor->getDocument();
}

/*
 *
 */
void CommandBasedConfigModelAccessor::setRootNode(std::auto_ptr<ConfigNode> newRoot)
{
	accessor->setRootNode(newRoot);
}

/*
 *
 */
ConfigNode* CommandBasedConfigModelAccessor::getRootNode()
{
	return accessor->getRootNode();
}

/**
 * Searches the ConfigNode designated by a path.
 * @param path the path to search in the model.
 * @return Returns the ConfigNode matching the path if found; and null otherwise.
 */
ConfigNode* CommandBasedConfigModelAccessor::pathToNode(const Path& path) const
{
	return accessor->pathToNode(path);
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setNameAndIndex(ConfigNode* node, const std::string& newNameAndIndex)
{
	return execute(new RenameConfigNode(node, newNameAndIndex, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setName(ConfigNode* node, const std::string& newName)
{
	PathNode nodeName(newName, node->getIndex());
	return execute(new RenameConfigNode(node, nodeName.toString(), accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setIndex(ConfigNode* node, const std::string& newIndex)
{
	PathNode nodeName(node->getName(), newIndex);
	return execute(new RenameConfigNode(node, nodeName.toString(), accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setValue(ConfigNode* node, const Value& newValue)
{
	return execute(new SetConfigNodeValue(node, newValue, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setUnit(ConfigNode* node, const std::string& newUnit)
{
	return execute(new SetConfigNodeUnit(node, newUnit, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setRange(ConfigNode* node, RangePtr range)
{
	return execute(new SetConfigNodeRange(node, range, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::setIncludeUrl(ConfigNode* node, const std::string& url)
{
	return execute(new SetConfigNodeIncludeUrl(node, url, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::insertInside(ConfigNode* parent, std::auto_ptr<Attribute> newAttribute)
{
	return execute(new InsertConfigNode(parent,
		InsertConfigNode::Child,
		newAttribute,
		accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::override(ConfigNode*& node, std::vector<ConfigNode*>& /*createdNodes*/)
{
	return execute(new OverrideConfigNode(node, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::insertAfter(ConfigNode* node, std::auto_ptr<Attribute> newAttribute)
{
	return execute(new InsertConfigNode(node,
		InsertConfigNode::Sibling,
		newAttribute,
		accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::remove(ConfigNode* parent, const IndexList& indexes, Object* /*removed*/)
{
	return execute(new DeleteConfigNode(parent,
		indexes,
		accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::remove(ConfigNode* parent, unsigned index, Object* /*removed*/)
{
	IndexList indexes(1);
	indexes[0] = index;
	return execute(new DeleteConfigNode(parent,
		indexes,
		accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::copy(ConfigNode* /*source*/, ConfigNode* /*dest*/)
{
	// TODO
	return false;
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::move(ConfigNode* /*source*/,
	const IndexList& /*srcIndexes*/,
	ConfigNode* /*dest*/, unsigned /*insertionIndex*/)
{
	// TODO
	return false;
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::include(ConfigNode* node, const std::string& url)
{
	return execute(new IncludeConfig(node, url, accessor));
}

/*
 *
 */
bool CommandBasedConfigModelAccessor::import(ConfigNode* node, const std::string& url)
{
	return execute(new ImportConfig(node, url, accessor));
}


///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::UndoableConfigModelAccessor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 * @param parent the owner of this new object.
 */
UndoableConfigModelAccessor::UndoableConfigModelAccessor(QObject* parent)
	: QObject(parent), CommandBasedConfigModelAccessor(new ConfigModelAccessor(this)),
	undoStack(new QUndoStack(this))
{
}

/**
 * Destructor.
 */
UndoableConfigModelAccessor::~UndoableConfigModelAccessor()
{
}

/*
 *
 */
void UndoableConfigModelAccessor::setDocument(std::auto_ptr<Io::Document> newDocument)
{
	CommandBasedConfigModelAccessor::setDocument(newDocument);
	undoStack->clear();
}

/*
 *
 */
void UndoableConfigModelAccessor::setRootNode(std::auto_ptr<ConfigNode> newRoot)
{
	CommandBasedConfigModelAccessor::setRootNode(newRoot);
	undoStack->clear();
}

/**
 * Executes and stores an undoable command (taking care of a few details).
 * @param cmd the command to execute.
 */
bool UndoableConfigModelAccessor::execute(Command::AbstractConfigModelCommand* cmd)
{
	undoStack->push(cmd);
	bool succeeded = cmd->hasSucceeded();
	if(not succeeded)
	{
		undoStack->undo();
	}
	return succeeded;
}

}
}
