#include "AbstractConfigModelAccessor.h"
#include "CCfg/Io/Document.h"

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::AbstractConfigModelAccessor class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Destructor.
 */
AbstractConfigModelAccessor::~AbstractConfigModelAccessor()
{
}

/**
 * @fn void AbstractConfigModelAccessor::setRootNode(std::auto_ptr<ConfigNode> newRoot)
 * Redefines the underlying configuration root node.
 * @param newRoot the new configuration root.
 */

/**
 * @fn ConfigNode* AbstractConfigModelAccessor::getRootNode()
 * Returns the root of the underlying configuration.
 */

/**
 * @fn ConfigNode* AbstractConfigModelAccessor::pathToNode(const Path& path) const
 * Searches the ConfigNode designated by a path.
 * @param path the path to search in the model.
 * @return Returns the ConfigNode matching the path if found; and null otherwise.
 */

/**
 * @fn bool AbstractConfigModelAccessor::setNameAndIndex(ConfigNode* node, const std::string& newNameAndIndex)
 * Changes the name and the index of a configuration node.
 * @param node the node to alter.
 * @param newNameAndIndex a string of the form: name ['[' index ']']
 */

/**
 * @fn bool AbstractConfigModelAccessor::setName(ConfigNode* node, const std::string& newName)
 * Changes the name of a configuration node.
 * @param node the node to alter.
 * @param newName the new name to set.
 */

/**
 * @fn bool AbstractConfigModelAccessor::setIndex(ConfigNode* node, const std::string& newIndex)
 * Changes the index of a configuration node.
 * @param node the node to alter.
 * @param newIndex the new index to set.
 */

/**
 * @fn bool AbstractConfigModelAccessor::setValue(ConfigNode* node, const Value& value)
 * Changes the value of a configuration node.
 * @param node the node to alter.
 * @param value the new value to set.
 */

/**
 * @fn bool AbstractConfigModelAccessor::setUnit(ConfigNode* node, const std::string& unit)
 * Changes the unit of a configuration node.
 * @param node the node to alter.
 * @param unit the new unit to set.
 */

/**
 * @fn bool AbstractConfigModelAccessor::setRange(ConfigNode* node, RangePtr range)
 * Changes the range of a configuration node.
 */

/**
 * @fn bool AbstractConfigModelAccessor::setIncludeUrl(ConfigNode* node, const std::string& url)
 * Defines the include URL of a configuration node.
 * @param node the sub-tree below this node will be saved in the specified resource
 * at the next save, and will disappear from the main file. Instead, there will
 * just be a reference.
 * @param url the URL of the resource where to save the sub-tree.
 */

/**
 * @fn bool AbstractConfigModelAccessor::insertInside(ConfigNode* parent, std::auto_ptr<Attribute> newAttribute)
 * Inserts an attribute as the first child of a node.
 * @param parent the parent for the new attribute.
 * @param newAttribute the attribute to insert. The configuration tree takes ownership of it.
 */

/**
 * @fn bool AbstractConfigModelAccessor::insertAfter(ConfigNode* node, std::auto_ptr<Attribute> newAttribute)
 * Inserts an attribute after a node.
 * @param node the node after which the new attribute is inserted.
 * @param newAttribute the attribute to insert. The configuration tree takes ownership of it.
 */

/**
 * @fn bool AbstractConfigModelAccessor::override(ConfigNode* node, std::vector<ConfigNode*>& createdNodes)
 * Overrides an inherited attribute.
 * @param node the inherited ConfigNode.
 */

/**
 * @fn bool AbstractConfigModelAccessor::remove(ConfigNode* parent, const IndexList& indexes, Object* removed)
 * Removes an arbitrary list of nodes, children of the same configuration node.
 * @param parent the parent of the attributes to remove.
 * @param indexes the indexes of the nodes to remove.
 * @param removed [inout] an object (preferably empty) where removed attributes
 * will be put. It can be null if the caller doesn't care about retrieving them.
 */

/**
 * @fn bool AbstractConfigModelAccessor::remove(ConfigNode* parent, unsigned index, Object* removed)
 * Removes a single node of the specified parent node.
 * @param parent the parent of the attributes to remove.
 * @param index the index of the node to remove.
 * @param removed [inout] an object (preferably empty) where removed attributes
 * will be put. It can be null if the caller doesn't care about retrieving them.
 */

/**
 * @fn bool AbstractConfigModelAccessor::copy(ConfigNode* source, ConfigNode* dest)
 * Copies the source node subtree and pastes it as the last child of the destination.
 * @param source the source node to copy (with its subtree)
 * @param dest the destination node (where the source will be appended).
 */

/**
 * @fn bool AbstractConfigModelAccessor::move(ConfigNode* source, const IndexList& srcIndexes, ConfigNode* dest, unsigned insertionIndex)
 * Moves a list of nodes into another parent.
 * @param source the parent of the nodes to move.
 * @param srcIndexes the list of indexes (in the source node) of the nodes to move.
 * The nodes will keep the order of this list.
 * @param dest the new parent of the nodes, after they move.
 * @param insertionIndex the index of the first node in the destination, after they move.
 * The other nodes will have increasing consecutive indexes after this one.
 */

/**
 * @fn bool AbstractConfigModelAccessor::include(ConfigNode* node, const std::string& url)
 * Includes a configuration in a node.
 * @param node the node where the resource will be included. It must be a
 * configuration object. It will be cleared as part of the action.
 * @param url the URL of the resource to load and include.
 */

/**
 * @fn bool AbstractConfigModelAccessor::import(ConfigNode* node, const std::string& url)
 * Imports a configuration in a node.
 * @param node the node where the resource will be imported. It must be a
 * configuration object. It will be cleared as part of the action.
 * @param url the URL of the resource to import.
 */

}
}
