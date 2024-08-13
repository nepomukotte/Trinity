#include "ConfigNode.h"
#include "CCfg/View/InheritedIterator.h"
#include "CCfg/View/Search.h"
#include "CCfg/Path.h"
#include "CCfg/PathNode.h"
#include "CCfg/Object.h"
#include "CCfg/Scalar.h"
#include "CCfg/Attribute.h"
#include "CCfg/Comment.h"
#include "CCfg/Range.h"
#include "Utils/Error/ExcepHandling.h"
#include <sstream>
#include <iostream>

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigNode class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param rowInParent the row index of the new node in his parent.
 * @param parent the parent of the new node.
 */
ConfigNode::ConfigNode(int rowInParent, const ConfigNode *parent)
	: nodeConfig(), rowInParent(rowInParent), parent(const_cast<ConfigNode*>(parent)),
	children(),	viewResolution(), override(false),
	populated(false), resolved(false)
{
}

/**
 * Copy constructor.
 * Keep the parent pointer but do not modify the parent's list. Children are
 * not copied, they will be recomputed when accessed.
 * @param node the node to copy.
 */
ConfigNode::ConfigNode(const ConfigNode& node)
	: nodeConfig(node.nodeConfig), rowInParent(node.rowInParent), parent(node.parent),
	children(),	viewResolution(), override(false),
	populated(false), resolved(false)
{
}

/**
 * Destructor.
 */
ConfigNode::~ConfigNode()
{
}

/**
 * Assignment operator.
 * Keep the parent pointer but do not modify the parent's list. Children are
 * not copied, they we be recomputed when accessed.
 * @param node the node to copy.
 */
ConfigNode& ConfigNode::operator=(const ConfigNode& node)
{
	nodeConfig = node.nodeConfig;
	rowInParent = node.rowInParent;
	parent = node.parent;

	children.clear();
	viewResolution = View::Result();
	override = false;

	populated = false;
	resolved = false;

	return *this;
}

/**
 * Returns the name of this config node.
 */
const std::string& ConfigNode::getName() const
{
	return getConfigAttribute()->getName();
}

/**
 * Defines the name of this config node.
 */
void ConfigNode::setName(const std::string& newName)
{
	getConfigAttribute()->setName(newName);
}

/**
 * Returns the index of this config node.
 */
const std::string& ConfigNode::getIndex() const
{
	return getConfigAttribute()->getIndex();
}

/**
 * Defines the index of this config node.
 */
void ConfigNode::setIndex(const std::string& newIndex)
{
	getConfigAttribute()->setIndex(newIndex);
}

/**
 * Returns the full name of this node (name ['[' index ']']).
 */
std::string ConfigNode::getFullName() const
{
	PathNode nodeName(getName(), getIndex());
	return nodeName.toString();
}

/**
 * Tells whether this node is an object.
 * @return Returns true if this node is an object.
 */
bool ConfigNode::isObject() const
{
	return getConfigAttribute()->isObject();
}

/**
 * Tells whether this node have a valid value with respect to its range.
 * @returns Returns true if the value is within the range or if there are no range
 * or if it is not a scalar attribute.
 */
bool ConfigNode::hasValidValue() const
{
	const Attribute* attr = getConfigAttribute();
	return attr->isObject() or attr->getRange() == 0 or
		attr->getRange()->contains(attr->getValue());
}

/**
 * Tells whether this node and or subtree contains at least one non-inherited scalar.
 */
bool ConfigNode::hasScalarInstances() const
{
	struct RecursiveSearch
	{
		/**
		 * Nested recursive function to determine if the subtree contains a scalar.
		 */
		static bool containsScalar(const Attribute* parent)
		{
			for(Attribute::const_iterator curAttr = parent->begin(),
				endAttr = parent->end();
				curAttr != endAttr; ++curAttr)
			{
				if(curAttr->isObject())
				{
					if(containsScalar(&*curAttr))
					{
						return true;
					}
				}
				else
				{
					return true;
				}
			}
			return false;
		}
	};

	const Attribute* parent = getConfigAttribute();
	if(not parent->isObject() or isInherited())
	{
		// Scalars cannot contains anything, so there is no need to go further.
		// Inherited objects are not "instances" themselves from the point of
		// view of their parent so whatever scalar they may contain cannot be
		// considered as an "instance".
		return false;
	}
	return RecursiveSearch::containsScalar(parent);
}

/**
 * Returns the value of this configuration scalar.
 * @throws UnsupportedOperation this node does not represents a scalar.
 */
Value ConfigNode::getValue() const
{
	return getConfigAttribute()->getValue();
}

/**
 * Defines the value of this configuration scalar.
 * @throws UnsupportedOperation this node does not represents a scalar.
 */
void ConfigNode::setValue(const Value& value)
{
	getConfigAttribute()->setValue(value);
}

/**
 * Returns the unit of this configuration scalar.
 * @throws UnsupportedOperation this node does not represents a scalar.
 */
const std::string& ConfigNode::getUnit() const
{
	return getConfigAttribute()->getUnit();
}

/**
 * Defines the unit of this configuration scalar.
 * @throws UnsupportedOperation this node does not represents a scalar.
 */
void ConfigNode::setUnit(const std::string& newUnit)
{
	getConfigAttribute()->setUnit(newUnit);
}

/**
 * Returns the range of this configuration scalar.
 * @throws UnsupportedOperation this node does not represents a scalar.
 */
Range* ConfigNode::getRange() const
{
	return getConfigAttribute()->getRange();
}

/**
 * Defines the range of this configuration scalar.
 * @throws UnsupportedOperation this node does not represents a scalar.
 */
void ConfigNode::setRange(RangePtr range)
{
	getConfigAttribute()->setRange(range);
}

/**
 * Returns the include URL of this configuration object.
 * @throws UnsupportedOperation this node does not represents an object.
 */
const std::string& ConfigNode::getIncludeUrl() const
{
	return getConfigAttribute()->getIncludeUrl();
}

/**
 * Returns the first parent node having an include Url.
 */
const ConfigNode* ConfigNode::getIncludeRoot() const
{
	const ConfigNode* current = this;
	if(current->getParent() != 0)
	{
		do
		{
			current = current->getParent();
		}
		while(current->getParent() != 0 and current->getIncludeUrl().empty());
	}
	return current;
}

/**
 * Defines the include URL of this configuration object.
 * @throws UnsupportedOperation this node does not represents an object.
 */
void ConfigNode::setIncludeUrl(const std::string& url)
{
	getConfigAttribute()->setIncludeUrl(url);
}

/**
 * Returns the pre-comment of this node.
 */
Comment& ConfigNode::getPreComment() const
{
	return getConfigAttribute()->getPreComment();
}

/**
 * Returns the post-comment of this node.
 */
Comment& ConfigNode::getPostComment() const
{
	return getConfigAttribute()->getPostComment();
}

/**
 * Returns the path of this node.
 */
const Path& ConfigNode::getPath() const
{
	return nodeConfig.getPath();
}

/**
 * Returns the path of the underlying configuration attribute.
 */
Path ConfigNode::getInternalPath() const
{
	return getConfigAttribute()->getPath();
}

/**
 * Returns the root node.
 */
ConfigNode* ConfigNode::getRoot()
{
	ConfigNode* current = this;
	while(current->parent != 0)
	{
		current = current->parent;
	}
	return current;
}

/**
 * Returns an iterator on this node in its parent.
 * This function is mainly useful to implement the concept of Node in the
 * CCfg::Search class.
 */
ConfigNode::iterator ConfigNode::getIterator()
{
	return parent->getChildren().begin() + rowInParent;
}

/**
 * Returns an iterator on this node in its parent.
 * This function is mainly useful to implement the concept of Node in the
 * CCfg::Search class.
 */
ConfigNode::const_iterator ConfigNode::getIterator() const
{
	return parent->getChildren().begin() + rowInParent;
}

/**
 * Returns the list of children configuration node (all attributes enumerated by the view).
 */
std::vector<ConfigNode>& ConfigNode::getChildren() const
{
	ensurePopulated();
	return children;
}

/**
 * Returns an iterator on the first child.
 */
ConfigNode::iterator ConfigNode::begin()
{
	return getChildren().begin();
}

/**
 * Returns an iterator over the last child.
 */
ConfigNode::iterator ConfigNode::end()
{
	return getChildren().end();
}

/**
 * Returns an iterator on the first child.
 */
ConfigNode::const_iterator ConfigNode::begin() const
{
	return getChildren().begin();
}

/**
 * Returns an iterator over the last child.
 */
ConfigNode::const_iterator ConfigNode::end() const
{
	return getChildren().end();
}


/**
 * Tells whether the underlying configuration attribute overrides another.
 * @return Returns true if the underlying attribute override another one.
 */
bool ConfigNode::isOverride() const
{
	ensureResolved();
	return override;
}

/**
 * Returns the backtrack level that was necessary to retrieve the configuration
 * attribute.
 * This getter uses lazy initialization.
 * @return Returns the backtrack level.
 */
int ConfigNode::getBacktrackLevel() const
{
	ensureResolved();
	return viewResolution.getBacktrackLevel();
}

/**
 * Defines config element associated to this node.
 * @param config the new config for the node.
 */
void ConfigNode::setConfig(const CConfig& config)
{
	nodeConfig = config;
	recomputeConfig();
}

/**
 * Recursively searches the node matching a path.
 * @param path the path to the node to find, either relative of absolute.
 * @param iPath an iterator on the next path node to search in *this.
 * @return Returns the node matching the remaining of the path; returns 0
 * if it cannot be found.
 */
ConfigNode* ConfigNode::search(const Path& path, Path::const_iterator iPath) const
{
	ConfigNode* next = search(*iPath);
	if(next != 0)
	{
		if(++iPath == path.end())
		{
			return next;
		}
		else
		{
			return next->search(path, iPath);
		}
	}
	else
	{
		return 0;
	}
}

/**
 * Searches the child node matching the PathNode.
 * @param node the PathNode (name and index) of the ConfigNode to search.
 * @return Returns the matching node; returns 0 if it cannot be found.
 */
ConfigNode* ConfigNode::search(const PathNode& node) const
{
	std::vector<ConfigNode>& children = getChildren();
	for(size_t i = 0; i < children.size(); ++i)
	{
		Attribute* attr = children[i].getConfigAttribute();
		if(attr->getName() == node.name and
			attr->getIndex() == node.index)
		{
			return &children[i];
		}
	}
	return 0;
}

/**
 * Searches the child node representing the specified Attribute.
 * @param attribute the attribute to search.
 * @return Returns the matching node; returns 0 if it cannot be found.
 */
ConfigNode* ConfigNode::search(const Attribute* attribute) const
{
	std::vector<ConfigNode>& children = getChildren();
	for(size_t i = 0; i < children.size(); ++i)
	{
		if(children[i].getConfigAttribute() == attribute)
		{
			return &children[i];
		}
	}
	return 0;
}

/**
 * Returns an iterator on the child node with the specified name and index.
 */
ConfigNode::iterator ConfigNode::search(const std::string& name, const std::string& index)
{
	for(iterator iChild = begin(); iChild != end(); ++iChild)
	{
		if(iChild->getName() == name and iChild->getIndex() == index)
		{
			return iChild;
		}
	}
	return end();
}

/**
 * Returns an iterator on the child node with the specified name and index.
 */
ConfigNode::const_iterator ConfigNode::search(const std::string& name, const std::string& index) const
{
	for(const_iterator iChild = begin(); iChild != end(); ++iChild)
	{
		if(iChild->getName() == name and iChild->getIndex() == index)
		{
			return iChild;
		}
	}
	return end();
}

/**
 * Returns the list of nodes that this node overrides.
 */
std::vector<ConfigNode*> ConfigNode::getOverriden()
{
	std::vector<ConfigNode*> overriden;

	if(not isInherited()) // pure inherited nodes are not overriden.
	{
		View::Search<ConfigNode> s(&inheritedFilter);
		ConfigNode* current = this;
		while((current = s.searchNext(getPath(), *current)))
		{
			// the overriden node must share the same index
			// (star indexes nodes does not count for overrides unless the index of *this is a star)
			if(current->getIndex() == getIndex())
			{
				overriden.push_back(current);
			}
		}
	}

	return overriden;
}

/**
 * Returns the list of nodes that overrides this one.
 */
std::vector<ConfigNode*> ConfigNode::getOverriding()
{
	std::vector<ConfigNode*> overriding;

	if(not isInherited()) // pure inherited nodes cannot be overriden
	{
		View::Search<ConfigNode> s(&inheritedFilter);
		ConfigNode* current = this;
		while((current = s.matchNext(getPath(), *current)))
		{
			if(current->getIndex() == getIndex())
			{
				overriding.push_back(current);
			}
		}
	}

	return overriding;
}

/**
 * Returns the list of nodes that this one inherits from.
 * Only configuration objects can inherit from something,
 * if *this is not an object the result will be an empty list.
 */
std::vector<ConfigNode*> ConfigNode::getInherited()
{
	std::vector<ConfigNode*> inherited;

	// the notion of inheritance only applies to objects. Although an inherited
	// object has a definition, it does not really inherits from its definition,
	// its just a replication, so we exclude them from this list.
	if(isObject() and not isInherited())
	{
		View::Search<ConfigNode> s(&inheritedFilter);
		ConfigNode* current = this;
		while((current = s.searchNext(getPath(), *current)))
		{
			inherited.push_back(current);
		}
	}

	return inherited;
}

/**
 * Returns the list of nodes that inherits from this one.
 * Only configuration objects can be inherited, if *this is not an object the
 * result will be an empty list.
 */
std::vector<ConfigNode*> ConfigNode::getInheriting()
{
	std::vector<ConfigNode*> inheriting;

	if(isObject() and not isInherited())
	{
		View::Search<ConfigNode> s(&inheritedFilter);
		ConfigNode* current = this;
		while((current = s.matchNext(getPath(), *current)))
		{
			inheriting.push_back(current);
		}
	}

	return inheriting;
}

/**
 * Returns the node where this inherited node is defined.
 * *this must be an inherited node.
 */
ConfigNode* ConfigNode::getDefinition()
{
	assert(isInherited());

	View::Search<ConfigNode> s(&inheritedFilter);
	ConfigNode* current = this;
	while((current = s.searchNext(getPath(), *current)))
	{
		// filter default objects
		if(current->getIndex() == getIndex())
		{
			return current; // return the first non-default object
		}
	}

	assert(false); // nothing found: it means the model is incoherent.
	return current;
}

/**
 * Inserts a new attribute after this one.
 * @param attribute the new attribute to insert.
 */
void ConfigNode::insertAfter(Attribute* attribute)
{
	ConfigNode* parent = getParent();
	Attribute* parentObj = parent->getConfigAttribute();

	Attribute::iterator iInsertion = getConfigAttribute()->getIterator();
	if(iInsertion != parentObj->end()) // if the object is not empty
	{
		parentObj->insert(++iInsertion, attribute);
	}
	else
	{
		parentObj->append(attribute);
	}
}

/**
 * Inserts a new attribute as the first child of this one.
 * @param attribute the new attribute to insert.
 */
void ConfigNode::insertInside(Attribute* attribute)
{
	Attribute* parentObj = getConfigAttribute();
	parentObj->prepend(attribute);
}

/**
 * Removes the underlying attribute of this configuration node.
 * ConfigNodes that stands for an inherited attribute cannot remove it.
 * @return Returns the removed attribute. The called can take ownership of it
 * doing an assignation, or drop the return result: in this case, it will be
 * automatically deleted.
 */
std::auto_ptr<Attribute> ConfigNode::remove()
{
	ConfigNode* parent = getParent();
	try
	{
		if(not isInherited() and parent != 0)
		{
			return parent->getConfigAttribute()->remove(getConfigAttribute());
		}
	}
	CATCH_AND_PRINT
	return std::auto_ptr<Attribute>();
}

/**
 * Returns the underlying configuration attribute.
 * This getter uses lazy initialization.
 * @return Returns the underlying configuration attribute.
 */
Attribute* ConfigNode::getConfigAttribute() const
{
	ensureResolved();
	return viewResolution.getAttribute();
}

/**
 * Populate the children's list by enumerating the view's content.
 */
void ConfigNode::populate() const
{
	if(isObject())
	{
		children.clear();

		int row = 0;
		CConfig::Iterator it = nodeConfig.iterate();
		while(it.hasNext())
		{
			it.next();
			ConfigNode child(row++, this);
			child.setConfig(it.get());
			children.push_back(child);
		}
	}
	populated = true;
}

/**
 * Resolves the view represented by the configuration node.
 */
void ConfigNode::resolve() const
{
	viewResolution = nodeConfig.resolve();

	// if we find this index in one of the default attributes, it means this
	// attribute overrides another attribute.
	std::string index = viewResolution->getIndex();

	bool hasNext = false;
	View::InheritedIterator def = nodeConfig.getInheritedAttributes();

	// browse the list of default object until we find one with the same index.
	while((hasNext = def.hasNext()) and def.next()->getIndex() != index) {}
	//             ^ yes this is an assignment: we need to keep the value of
	// hasNext() before we call next()

	// if hasNext() dit not break the loop, equality of indexes did.
	// so it overrides something.
	override = hasNext;
	resolved = true;
}

/**
 * Ask to recompute children attribute next time they will be accessed.
 */
void ConfigNode::recomputeChildren()
{
	children.clear();
	populated = false;
}

/**
 * Removes all children without removing the corresponding configuration attributes.
 * They won't be recomputed until the node it told to do so (with
 * recomputeChildren() or recomputeConfig()).
 */
void ConfigNode::clear()
{
	children.clear();
	populated = true;
	resolved = true;
}

/**
 * Ask to recompute the view next time its content will be accessed.
 */
void ConfigNode::recomputeConfig()
{
	resolved = false;
	recomputeChildren();
}

void ConfigNode::updateConfigPath()
{
	nodeConfig = CConfig(nodeConfig.getBaseAttribute());
	recomputeConfig();
}

/**
 * Ensures the node is populated, and populate it if necessary.
 */
void ConfigNode::ensurePopulated() const
{
	if(not populated)
	{
		populate();
	}
}

/**
 * Ensure the node view is resolved, and resolve it if necessary.
 */
void ConfigNode::ensureResolved() const
{
	if(not resolved)
	{
		resolve();
	}
}

/*
 * Filter which ignores nodes whose content has not been
 * computed.
 */
bool ConfigNode::populatedFilter(Model::ConfigNode& node)
{
	return not node.isPopulated();
}

/**
 * Filter which ignores inherited nodes.
 */
bool ConfigNode::inheritedFilter(Model::ConfigNode& node)
{
	return node.isInherited();
}

/**
 * Converts a CCfg::Value to a QVariant.
 */
QVariant ConfigNode::convertToVariant(const Value& value)
{
	switch(value.getType())
	{
	case INTEGER:
		return QVariant(value.getInt());
	case HEXADECIMAL:
		return QVariant(value.getHex());
	case REAL:
		return QVariant(value.getReal());
	case BOOLEAN:
		return QVariant(value.getBool());
	case STRING:
		return QVariant(value.getString().c_str());
	default:
		assert(false);
		return QVariant();
	}
}

Value ConfigNode::convertFromVariant(const QVariant& var)
{
	switch(var.type())
	{
	case QVariant::Int:
		return Value(var.toInt());
	case QVariant::UInt:
		return Value(var.toUInt());
	case QVariant::Double:
		return Value(var.toDouble());
	case QVariant::Bool:
		return Value(var.toBool());
	case QVariant::String:
		return Value(var.toString().toStdString());
	default:
		return Value();
	}
}

}
}


// instante template View::Search for ConfigNode.
#include "CCfg/View/Search.hpp"

namespace CCfg
{
namespace View
{
// explicit instanciation
template class Search<Model::ConfigNode>;
}
}
