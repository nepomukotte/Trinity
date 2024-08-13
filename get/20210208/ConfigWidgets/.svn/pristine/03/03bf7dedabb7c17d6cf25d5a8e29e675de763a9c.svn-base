#include "ConfigModelAccessor.h"
#include "ConfigNode.h"
#include "Utils/Error/Reporter.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Io/Url.h"
#include "CCfg/View/ClassIterator.h"
#include "CCfg/View/Search.h"
#include "CCfg/Object.h"
#include "CCfg/Range.h"
#include <utl/Logging.h>
#include <sstream>

namespace CCfg
{
namespace Model
{
using namespace Utils;
///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModelAccessor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 * @param parent the owner of this new object.
 */
ConfigModelAccessor::ConfigModelAccessor(QObject* parent)
	: QObject(parent)
{
}

/**
 * Destructor.
 */
ConfigModelAccessor::~ConfigModelAccessor()
{
}

/*
 *
 */
void ConfigModelAccessor::setRootNode(std::auto_ptr<ConfigNode> newRoot)
{
	if(rootNode.get() != 0)
	{
		emit beginClear(0);
		rootNode.reset();
		emit endClear();
	}

	emit beginRecompute(0);
	rootNode = newRoot;
	emit endRecompute();
}

/*
 *
 */
ConfigNode* ConfigModelAccessor::getRootNode()
{
	return rootNode.get();
}

/*
 *
 */
void ConfigModelAccessor::setDocument(std::auto_ptr<Io::Document> newDocument)
{
	doc = newDocument;

	std::auto_ptr<ConfigNode> root(new ConfigNode);
	root->setConfig(CConfig(doc->getConfig()));
	setRootNode(root);
}

/*
 *
 */
const Io::Document* ConfigModelAccessor::getDocument() const
{
	return doc.get();
}

/*
 *
 */
Io::Document* ConfigModelAccessor::getDocument()
{
	return doc.get();
}

/*
 *
 */
ConfigNode* ConfigModelAccessor::pathToNode(const Path& path) const
{
	Path::const_iterator iPath = path.begin();
	if(rootNode->getName() == iPath->name and
		rootNode->getIndex() == iPath->index)
	{
		++iPath;
		if(iPath != path.end())
		{
			return rootNode->search(path, iPath);
		}
		else
		{
			return rootNode.get();
		}
	}
	else
	{
		return 0;
	}
}

/*
 *
 */
bool ConfigModelAccessor::setNameAndIndex(ConfigNode* node, const std::string& newNameAndIndex)
{
	PathNode nodeName;
	// parses the string
	if(nodeName.fromString(newNameAndIndex)) // successful parsing
	{
		if(nodeName.name == node->getName()) // the name has not changed
		{
			if(nodeName.index == node->getIndex())
			{
				return true; // nothing has changed
			}
			else
			{
				// only the index has changed
				return setIndex(node, nodeName.index);
			}
		}
		else // the name has changed
		{
			if(nodeName.index == node->getIndex())
			{
				// only the name has changed
				return setName(node, nodeName.name);
			}
			else
			{
				// both the name and the index have changed
				ConfigNode* parent = node->getParent();
				// check if (name, index) is not already used.
				if(not warnIfAlreadyExist(parent, nodeName))
				{
					ConfigNode backupNode(*node);
					Update update(this, parent);

					backupNode.setName(nodeName.name);
					backupNode.setIndex(nodeName.index);
					if(node == rootNode.get())
					{
						rootNode->updateConfigPath();
					}
					return true;
				}
			}
		}
	}
	else
	{
		InvalidAttributeName error(newNameAndIndex);
		Error::Reporter::report(error);
	}
	return false;
}

/*
 *
 */
bool ConfigModelAccessor::setName(ConfigNode* node, const std::string& newName)
{
	if(node->isRoot())
	{
		// the root requires a special case because it does not has a parent.
		//ConfigNode backupNode(*node);
		Update update(this, 0);
		node->setName(newName);
		node->updateConfigPath(); // Updates the view path of the root node.
		// Note: contrary to regular nodes, Update does not destroy the root node
		// so we can use it directly.
	}
	else
	{
		// Create the list of all instances of the class of the renamed node.
		// Renaming is not immediate, because if a name collision occurs, we can
		// detect it before any change is made.
		std::vector<Attribute*> instances;
		View::ClassIterator iInst(*rootNode->getConfigAttribute(), node->getPath());
		while(iInst.hasNext())
		{
			iInst.next();
			PathNode nodeName(newName, iInst->getIndex());
			// check name collision
			if(not warnIfAlreadyExist(iInst->getParent(), nodeName))
			{
				instances.push_back(iInst.get());
			}
			else
			{
				return false; // fail in case of collision
			}
		}

		for(size_t i = 0; i < instances.size(); ++i)
		{
			// find the ConfigNode associated with the Attribute
			// (because ClassIterator only works at the CompoundConfig level)
			Attribute* instance = instances[i];
			const Path& path = instance->getPath();
			ConfigNode* instanceNode = rootNode->search(path, ++path.begin());
			assert(instanceNode != 0);

			// rename the ConfigNode
			ConfigNode backupNode(*instanceNode);
			Update update(this, instanceNode->getParent());
			backupNode.setName(newName);
		}
	}
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::setIndex(ConfigNode* node, const std::string& newIndex)
{
	PathNode nodeName(node->getName(), newIndex);
	ConfigNode* parent = node->getParent();
	if(not warnIfAlreadyExist(parent, nodeName))
	{
		ConfigNode backupNode(*node);
		Update update(this, parent);
		backupNode.setIndex(nodeName.index);
		if(node == rootNode.get())
		{
			rootNode->updateConfigPath();
		}
		return true;
	}
	return false;
}

/*
 *
 */
bool ConfigModelAccessor::setValue(ConfigNode* node, const Value& value)
{
	Value newValue = value;
	if(newValue.getType() == node->getValue().getType())
	{
		ConfigNode backupNode(*node);
		{
			Update update(this, node->getParent());
			backupNode.setValue(newValue);
		}
		return true;
	}
	else
	{
		ForbiddenTypeChange error(node->getValue().getType(), newValue.getType());
		Error::Reporter::report(error);
		return false;
	}
}

/*
 *
 */
bool ConfigModelAccessor::setUnit(ConfigNode* node, const std::string& unit)
{
	ConfigNode backupNode(*node);
	{
		Update update(this, node->getParent());
		backupNode.setUnit(unit);
	}
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::setRange(ConfigNode* node, RangePtr range)
{
	ConfigNode backupNode(*node);
	{
		Update update(this, node->getParent());
		backupNode.setRange(range);
	}
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::setIncludeUrl(ConfigNode* node, const std::string& url)
{
	ConfigNode backupNode(*node);
	{
		Update update(this, node->getParent());
		backupNode.setIncludeUrl(url);
	}
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::insertAfter(ConfigNode* node, std::auto_ptr<Attribute> newAttribute)
{
	ConfigNode* parent = node->getParent();
	if(not warnIfInherited(parent))
	{
		PathNode nodeName(*newAttribute.get());

		if(not warnIfAlreadyExist(parent, nodeName))
		{
			// the update transaction will destroy *node and all its siblings,
			// therefore we have to make a copy to use it during the transaction.
			ConfigNode backupNode(*node);
			Attribute* newAttr = newAttribute.release();
			IndexList index(1);

			emit beginInsertion(parent);
			{
				Update update(this, parent);
				backupNode.insertAfter(newAttr);
			}
			ConfigNode* insertedNode = parent->search(newAttr);
			assert(insertedNode != 0);
			index[0] = insertedNode->getRowInParent();
			emit endInsertion(parent, index);
			return true;
		}
	}
	return false;
}

/*
 *
 */
bool ConfigModelAccessor::insertInside(ConfigNode* parent, std::auto_ptr<Attribute> newAttribute)
{
	if(not warnIfInherited(parent))
	{
		PathNode nodeName(*newAttribute.get());
		if(not warnIfAlreadyExist(parent, nodeName))
		{
			Attribute* newAttr = newAttribute.release();
			IndexList index(1);

			emit beginInsertion(parent);
			{
				Update update(this, parent);
				parent->insertInside(newAttr);
			}
			ConfigNode* insertedNode = parent->search(newAttr);
			assert(insertedNode != 0);
			index[0] = insertedNode->getRowInParent();
			emit endInsertion(parent, index);
			return true;
		}
	}
	return false;
}

/*
 *
 */
bool ConfigModelAccessor::override(ConfigNode*& node, std::vector<ConfigNode*>& createdNodes)
{
	bool result = false;
	if(node->isInherited())
	{
		const std::string& name = node->getName();
		const std::string& index = node->getIndex();
		ConfigNode* parent = node->getParent();

		bool isTopOverride = !override(parent, createdNodes);
		if(!isTopOverride)
		{
			node = &*parent->search(name, index); // necessary because node have been invalidated
		}

		for(ConfigNode::iterator sibling = parent->begin(), endSibling = parent->end(); sibling != endSibling; ++sibling)
		{
			if((node->getIndex() == CCfg::DEFAULT_OBJ_INDEX or sibling->getIndex() == node->getIndex())
				and sibling->getName() == node->getName())
			{
				if(!sibling->isOverride() or sibling->isInherited())
				{
					std::auto_ptr<Attribute> newAttribute(sibling->getConfigAttribute()->clone(SHALLOW));
					if(newAttribute->isObject())
					{
						newAttribute->setIncludeUrl(std::string()); // clear include url, otherwise the same file would be used at several places.
					}
					if(insertInside(parent, newAttribute))
					{
						LOG_DEBUG() << "created " << sibling->getPath();
						result = true;
						if(isTopOverride)
						{
							createdNodes.push_back(&*sibling);
						}
					}
					else
					{
						LOG_WARN() << "failed creating " << sibling->getPath();
					}
				}
			}
		}
		if(result)
		{
			node = &*parent->search(name, index); // necessary because node have been invalidated
		}
	}
	return result;
}

//bool ConfigModelAccessor::unOverride(ConfigNode* node, Object* removed)
//{
//	struct InstanceSearchVisitor : CCfg::Visitor
//	{
//		InstanceSearchVisitor(const CCfg::Attribute* excluded) : hasInstance(false), excluded(excluded) {}
//		bool onObjectEnter(const CCfg::Object& objectCfg) { return &objectCfg != excluded; }
//		bool onObjectExit(const CCfg::Object& /*objectCfg*/) { return true; }
//		bool onScalar(const CCfg::Scalar& scalar)	{ hasInstance = (&scalar != excluded); return !hasInstance; }
//
//		bool hasInstance;
//		const CCfg::Attribute* excluded;
//	};
//	ConfigNode* parent = node->getParent();
//	if(parent != 0)
//	{
//		CCfg::Attribute* parentObj = parent->getConfigAttribute();
//		InstanceSearchVisitor searcher(node->getConfigAttribute());
//		parentObj->accept(searcher);
//		if(!searcher.hasInstance and parent->isOverride())
//		{
//			unOverride(parent, removed);
//		}
//		else
//		{
//			LOG_INFO() << "removing branch " << node->getPath();
//			IndexList index;
//			index.push_back(node->getRowInParent());
//			remove(parent, index, removed);
//		}
//	}
//}

/*
 *
 */
bool ConfigModelAccessor::remove(ConfigNode* parent, const IndexList& indexes, Object* removed)
{
	if(warnIfInherited(parent))
	{
		return false;
	}

	std::vector<ConfigNode>& children = parent->getChildren();
	// backs up the nodes to remove (because the update mechanism destroy them).
	std::vector<ConfigNode> nodesToRemove(indexes.size());
	for(size_t i = 0; i < indexes.size(); ++i)
	{
		nodesToRemove[i] = children[indexes[i]];
	}

	emit beginRemoval(parent, indexes);
	{
		Update update(this, parent);
		for(size_t i = 0; i < nodesToRemove.size(); ++i)
		{
			std::auto_ptr<Attribute> attr = nodesToRemove[i].remove();
			if(removed != 0)
			{
				removed->append(attr.release());
			}
		}
	}
	emit endRemoval(parent);
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::remove(ConfigNode* parent, unsigned index, Object* removed)
{
	if(warnIfInherited(parent))
	{
		return false;
	}
	ConfigNode target(parent->getChildren()[index]);

	IndexList indexes(1);
	indexes[0] = index;
	emit beginRemoval(parent, indexes);
	{
		Update update(this, parent);
		std::auto_ptr<Attribute> attr = target.remove();
		if(removed != 0)
		{
			removed->append(attr.release());
		}
	}
	emit endRemoval(parent);
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::copy(ConfigNode* /*source*/, ConfigNode* /*dest*/)
{
	//Update update(this, dest);
	// TODO
	return false;
}

/*
 *
 */
bool ConfigModelAccessor::move(ConfigNode* /*source*/,
	const IndexList& /*srcIndexes*/,
	ConfigNode* /*dest*/, unsigned /*insertionIndex*/)
{
	//Update update(this, dest);
	// TODO
	return false;
}

/*
 *
 */
bool ConfigModelAccessor::include(ConfigNode* node, const std::string& url)
{
	ConfigNode backupNode(*node);
	{
		Update update(this, node->getParent());
		doc->include(backupNode.getConfigAttribute()->toObject(), Io::Url(url));
	}
	return true;
}

/*
 *
 */
bool ConfigModelAccessor::import(ConfigNode* node, const std::string& url)
{
	ConfigNode backupNode(*node);
	{
		Update update(this, node->getParent());
		doc->import(backupNode.getConfigAttribute()->toObject(), Io::Url(url));
	}
	return true;
}

/**
 * Inserts a list of attributes in a ConfigNode at arbitrary position.
 * @param parent the ConfigNode where attributes will be inserted.
 * @param attributes this CCfg::Object contains the list of attributes to insert.
 * It will not be inserted itself, and will be cleared during the method
 * execution (the CCfg::Object under the ConfigNode will take ownership of these
 * attributes).
 * @param insertedNodeIndex the post-insertion index of each attribute to insert.
 * Its size must be the same as the number of elements in "attributes".
 * for each i the index where to insert attributes[i] is insertedNodeIndex[i].
 */
void ConfigModelAccessor::insert(ConfigNode* parent, Object* attributes,
	const IndexList& insertedNodeIndex)
{
	if(warnIfInherited(parent))
	{
		return;
	}

	assert(insertedNodeIndex.size() == static_cast<unsigned>(attributes->size()));

	std::vector<ConfigNode>& children = parent->getChildren();

	IndexList existingNodeIndex =
		computePostInsertionIndexes(children.size(), insertedNodeIndex);

	// for each i, insertedNodeIndex[i] has a counterpart in inserterNodeIndex[i]
	// (remember: each node to insert needs an inserter)
	// the value in each inserterNodeIndex[i] is an index in inserterNode vector.
	// this is needed for inserter sharing.
	IndexList inserterNodeIndex(insertedNodeIndex.size());
	std::vector<ConfigNode> inserterNode =
		computeInserterNodes(parent, existingNodeIndex, insertedNodeIndex, inserterNodeIndex);

	beginInsertion(parent);
	{
		Update update(this, parent);

		// for each node to insert, in backward order.
		// this order is necessary because as we cannot use restored attributes
		// as inserters, we have to begin by inserting the latest objects for
		// each inserter.
		// eg: if we have A B C D in a parent and delete B C, we'll have to use
		// A to insert them. As we only have insertAfter method, we first insert
		// C after A, which gives: A C D. We finish by inserting B after A,
		// which gives: A B C D.
		for(int i = insertedNodeIndex.size()-1; i >= 0; --i)
		{
			std::auto_ptr<Attribute> attr = attributes->remove(--attributes->end());
			if(inserterNodeIndex[i] == 0xFFFFFFFF)
			{
				parent->insertInside(attr.release());
			}
			else
			{
				assert(inserterNodeIndex[i] < inserterNode.size());
				ConfigNode& inserter = inserterNode[inserterNodeIndex[i]];
				inserter.insertAfter(attr.release());
			}
		}
	}
	endInsertion(parent, insertedNodeIndex);
}

/**
 * Computes the indexes of an array element after the insertion of several
 * elements at arbitrary position.
 * Example: if we have 5 elements in an array and insert 2 element that should
 * be at index 2 and 4 after insertion, the existing element will have the
 * following indexes after insertion: 0 1 3 5 6
 * @param elementCount the number of elements in the array.
 * @param insertedElementIndex an array containing the indexes of elements to insert.
 * These indexes are the indexes that they should have after insertion.
 * @return returns an array of the post insertion indexes of existing elements.
 */
IndexList ConfigModelAccessor::computePostInsertionIndexes(
	unsigned elementCount, const IndexList& insertedElementIndex) const
{
	// post-insertion indexes of existing elements.
	// i -> existingElementIndex[i] does the mapping between pre insertion indexes
	// to post insertion indexes
	IndexList existingElementIndex(elementCount);
	// first: initialize with pre-insertion index)
	for(size_t i = 0; i < existingElementIndex.size(); ++i)
	{
		existingElementIndex[i] = i;
	}

	// then: compute post insertion indexes.
	for(size_t i = 0; i < insertedElementIndex.size(); ++i) // for each inserted element
	{
		for(size_t j = 0; j < existingElementIndex.size(); ++j) // for each existing element
		{
			// if the current inserted element has a lower or equal index than an existing element
			if(insertedElementIndex[i] <= existingElementIndex[j])
			{
				// the index of the existing element is shifted by one.
				++existingElementIndex[j];
			}
		}
	}
	return existingElementIndex;
}

/**
 * Computes the list of nodes to use for the insertion of several ConfigNodes at
 * arbitrary position.
 * @param parent the parent in which nodes will be inserted.
 * @param existingNodeIndex the post-insertion index of the existing nodes.
 * @param insertedNodeIndex the post-insertion index of the nodes to insert.
 * @param inserterNodeIndex [out] this list will be filled with the index of the
 * inserter (in the returned list) for each node to insert. In other word,
 * for each i, insertedNodeIndex[i] will have to be inserted using the inserter
 * at position inserterNodeIndex[i] in the return list.
 * If the number is 0xFFFFFFFF it means that the node will have to be inserted
 * at the beginning of the parent's nodes. Otherwise it means it will have to be
 * inserted _after_ the specified inserter node.
 * @return Returns the list of ConfigNodes required for a successful insertion.
 * An inserter can be used to insert several nodes. This array should be accessed
 * using the indexes in inserterNodeIndex.
 */
std::vector<ConfigNode> ConfigModelAccessor::computeInserterNodes(
	ConfigNode* parent,
	const IndexList& existingNodeIndex,
	const IndexList& insertedNodeIndex,
	IndexList& inserterNodeIndex) const
{
	std::vector<ConfigNode>& children = parent->getChildren();
	assert(children.size() == existingNodeIndex.size());

	std::vector<ConfigNode> inserterNode;
	// inserterNode.size() is always <= insertedNodeIndex.size() because each node to
	// insert needs only one inserter node. But an inserter can be shared by
	// several nodes to insert.
	inserterNode.reserve(insertedNodeIndex.size());

	// index of the current node to insert (valid in the insertedNodeIndex vector)
	size_t iInsert = 0;
	// index of the current existing node (iExist may be seen as a pre insertion index):
	// - valid in the existingNodeIndex when retrieving a post insertion index;
	// - valid in the children vector to retrieve an existing (pre insertion) node.
	size_t iExist = 0;

	size_t lastInserterIndex = 0xFFFFFFFF;
	size_t candidateInserterIndex = 0xFFFFFFFF;

	do
	{
		// if the current inserted node index overruns the current existing node index.
		if(iExist < existingNodeIndex.size() and
			insertedNodeIndex[iInsert] > existingNodeIndex[iExist])
		{
			// the candidate inserter cannot be an inherited node.
			if(children[iExist].isInherited())
			{
				++iExist; // try the next existing node.
			}
			else
			{
				// we store the index of the existing node, to remember it as a candidate inserter.
				candidateInserterIndex = iExist++;
			}
		}
		else // otherwise we use the current node of insertion to insert the node iInsert.
		{
			if(candidateInserterIndex == 0xFFFFFFFF) // no candidate inserter defined yet.
			{
				inserterNodeIndex[iInsert++] = 0xFFFFFFFF;
			}
			else
			{
				// if the candidate inserter is not the same as the previous
				if(candidateInserterIndex != lastInserterIndex)
				{
					// we need to backup a new inserter
					// (for node iInsert and maybe subsequent nodes).
					inserterNode.push_back(children[candidateInserterIndex]);
					lastInserterIndex = candidateInserterIndex;
				}
				inserterNodeIndex[iInsert++] = inserterNode.size() - 1; // the index of the current inserter.
			}
		}
	}
	while(iInsert < insertedNodeIndex.size());

	return inserterNode;
}

/**
 * Clears the content of a node.
 * @param node the node to clear.
 */
void ConfigModelAccessor::clear(ConfigNode* node)
{
	if(node == 0 or node->getChildren().size() > 0)
	{
		emit beginClear(node);
		if(node != 0)
		{
			node->clear();
		}
		else
		{
			rootNode->clear();
		}
		emit endClear();
	}
}

/**
 * Recompute the content of a node.
 * @param node the node to recompute.
 */
void ConfigModelAccessor::recompute(ConfigNode* node)
{
	ConfigNode* updatedNode = (node != 0) ? node : rootNode.get();

	// we must recompute before sending the signal, because the implementation
	// will need to know the new number of children nodes
	updatedNode->recomputeConfig();
	emit beginRecompute(node);

	emit endRecompute();
}

void ConfigModelAccessor::emitBeginUpdate()
{
	emit beginUpdate();
}

void ConfigModelAccessor::emitEndUpdate()
{
	emit endUpdate();
}

/**
 * Issues a warning message if a node of the specified name already exist in the
 * parent node.
 * @param parent the parent node where the method must perform its search.
 * @param newNodeName the name and index of the node to search.
 * @return Return true if the name already exists; false otherwise.
 */
bool ConfigModelAccessor::warnIfAlreadyExist(ConfigNode* parent, const PathNode& newNodeName)
{
	if(parent != 0)
	{
		return warnIfAlreadyExist(parent->getConfigAttribute(), newNodeName);
	}
	return false;
}

/**
 * Issues a warning message if an attribute of the specified name already exist in the
 * parent object.
 * @param parent the parent object where the method must perform its search.
 * @param newNodeName the name and index of the node to search.
 * @return Return true if the name already exists; false otherwise.
 */
bool ConfigModelAccessor::warnIfAlreadyExist(const Attribute* parent, const PathNode& newNodeName)
{
	if(parent != 0)
	{
		if(parent->search(newNodeName.name, newNodeName.index) != parent->end())
		{
			AttributeAlreadyExists error(newNodeName.toString());
			Error::Reporter::report(error);
			return true;
		}
	}
	return false;
}

/**
 * Issues a warning message if the specified node is inherited (saying that
 * operation is not permitted).
 */
bool ConfigModelAccessor::warnIfInherited(ConfigNode* node)
{
	if(node->isInherited())
	{
		OperationOnInherited error;
		Error::Reporter::report(error);
		return true;
	}
	else
	{
		return false;
	}
}



///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModelAccessor::Update class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor. Begins a model update by clearing the parent node and all
 * nodes inheriting from it.
 * @param accessor the object through which modifications on the model are carried out.
 * @param parent the parent of the updated node.
 */
ConfigModelAccessor::Update::Update(ConfigModelAccessor* accessor, ConfigNode* parent)
	: accessor(accessor), parent(parent), impactedNodes()
{
	accessor->emitBeginUpdate();
	accessor->clear(parent);
	if(parent != 0)
	{
		View::Search<ConfigNode> s(&ConfigNode::populatedFilter);
		ConfigNode* current = parent;
		while((current = s.matchNext(parent->getPath(), *current)))
		{
			if(current->isPopulated())
			{
				impactedNodes.push_back(current);
			}
		}

		for(std::vector<ConfigNode*>::const_iterator iNode = impactedNodes.begin();
			iNode != impactedNodes.end(); ++iNode)
		{
			accessor->clear(*iNode);
		}
	}
}

/**
 * Destructor. Ends a model update by recomputing the content of the parent
 * node and all nodes inheriting from it.
 */
ConfigModelAccessor::Update::~Update()
{
	accessor->recompute(parent);
	if(parent != 0)
	{
		for(std::vector<ConfigNode*>::const_iterator iNode = impactedNodes.begin();
			iNode != impactedNodes.end(); ++iNode)
		{
			accessor->recompute(*iNode);
		}
	}
	accessor->emitEndUpdate();
}


////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModelAccessor::InvalidAttributeName class methods body
////////////////////////////////////////////////////////////////////////////////

ConfigModelAccessor::InvalidAttributeName::InvalidAttributeName(const std::string& fullName)
	: fullName(fullName)
{
}

std::string ConfigModelAccessor::InvalidAttributeName::formatTitle() const
{
	return "Input error";
}

std::string ConfigModelAccessor::InvalidAttributeName::formatBrief() const
{
	std::ostringstream stream;
	stream << "The new name '" << fullName << "' is not a valid attribute name.";
	return stream.str();
}

std::string ConfigModelAccessor::InvalidAttributeName::formatSuggestion() const
{
	return "An attribute name should begin with a letter or an underscore "
		"and may contain letters, digits, dashes or underscores.\n"
		"The optionnal index must be enclosed in brackets, "
		"and be either a number or a name matching "
		"previous requirements.";
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModelAccessor::OperationOnInherited class methods body
////////////////////////////////////////////////////////////////////////////////

std::string ConfigModelAccessor::OperationOnInherited::formatTitle() const
{
	return "Invalid operation";
}

std::string ConfigModelAccessor::OperationOnInherited::formatBrief() const
{
	return "This operation is not allowed on inherited attributes.";
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModelAccessor::AttributeAlreadyExists class methods body
////////////////////////////////////////////////////////////////////////////////

ConfigModelAccessor::AttributeAlreadyExists::AttributeAlreadyExists(const std::string& name)
	: attributeName(name)
{
}

std::string ConfigModelAccessor::AttributeAlreadyExists::formatTitle() const
{
	return "Input error";
}

std::string ConfigModelAccessor::AttributeAlreadyExists::formatBrief() const
{
	std::ostringstream stream;
	stream << "The attribute '" << attributeName << "' already exists in the parent object.";
	return stream.str();
}

std::string ConfigModelAccessor::AttributeAlreadyExists::formatSuggestion() const
{
	return "Please choose another name or another index";
}

////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigModelAccessor::ForbiddenTypeChange class methods body
////////////////////////////////////////////////////////////////////////////////

ConfigModelAccessor::ForbiddenTypeChange::ForbiddenTypeChange(ValueTypeId oldType,
		ValueTypeId newType)
	: oldType(oldType), newType(newType)
{
}

std::string ConfigModelAccessor::ForbiddenTypeChange::formatTitle() const
{
	return "Input error";
}

std::string ConfigModelAccessor::ForbiddenTypeChange::formatBrief() const
{
	return "Type changes are forbidden.";
}

std::string ConfigModelAccessor::ForbiddenTypeChange::formatDetails() const
{
	std::ostringstream stream;
	stream << "The new value interprets as a '" << TypeInfo::getDisplayName(newType) <<
		"' whereas the old one was a '" << TypeInfo::getDisplayName(oldType) << "'";
	return stream.str();
}

}
}
