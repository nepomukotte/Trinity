#ifndef CCFG_MODEL_CONFIGNODE_H_INCLUDED
#define CCFG_MODEL_CONFIGNODE_H_INCLUDED

#include "CCfg/CConfig.h"
#include <QVariant>
#include <vector>
#include <memory>

namespace CCfg
{
//forward declarations
struct PathNode;
class Path;
class Comment;
class Range;
typedef std::auto_ptr<Range> RangePtr;

namespace Model
{
//forward declarations
class ConfigModelAccessor;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * A class representing the configuration tree displayed in ConfigEditorTree.
 * This class allows to keep CConfig objects persistent, and organize them in a
 * recursive tree-like structure.
 * It also contains informations to browse the tree, provide direct access to
 * its element, and easily implement operations requested by QAbstractItemModel.
 *
 * @author fchateau
 */
class ConfigNode
{
	// fixed attributes
	CConfig nodeConfig;	///< the config element associated to this node.
	int rowInParent;	///< the row index of this node in his parent.
	ConfigNode* parent;	///< the parent of this node. The root is null.

	// cache attributes
	mutable std::vector<ConfigNode> children;	///< the list of children ("attributes").
	mutable View::Result viewResolution;		///< the view resolution result cache.
	mutable bool override;	///< indicates that the underlying attribute override another one.

	// cache flags
	mutable bool populated;	///< indicates whether the children list has been populated.
	mutable bool resolved;	///< indicates whether the view has been resolved.

public:
	typedef std::vector<ConfigNode>::iterator iterator;
	typedef std::vector<ConfigNode>::const_iterator const_iterator;

	// construction/destruction
	ConfigNode(int rowInParent = 0, const ConfigNode *parent = 0);
	ConfigNode(const ConfigNode& node);
	virtual ~ConfigNode();
	ConfigNode& operator=(const ConfigNode& node);

	// config definition
	void setConfig(const CConfig& config);

	// accessors
	const std::string& getName() const;
	const std::string& getIndex() const;
	std::string getFullName() const;
	Value getValue() const;
	const std::string& getUnit() const;
	Range* getRange() const;
	const std::string& getIncludeUrl() const;
	const ConfigNode* getIncludeRoot() const;
	bool isObject() const;
	bool hasValidValue() const;
	bool hasScalarInstances() const;
	Attribute* getConfigAttribute() const;

	Comment& getPreComment() const;
	Comment& getPostComment() const;

	const Path& getPath() const;
	Path getInternalPath() const;

	int getBacktrackLevel() const;

	// tree structure
	bool isRoot() const;
	ConfigNode* getRoot();
	ConfigNode* getParent();
	const ConfigNode* getParent() const;

	int getRowInParent() const;
	iterator getIterator();
	const_iterator getIterator() const;

	std::vector<ConfigNode>& getChildren() const;
	iterator begin();
	iterator end();
	const_iterator begin() const;
	const_iterator end() const;

	// search
	ConfigNode* search(const Path& path, Path::const_iterator iPath) const;
	ConfigNode* search(const PathNode& node) const;
	ConfigNode* search(const Attribute* attribute) const;
	iterator search(const std::string& name, const std::string& index);
	const_iterator search(const std::string& name, const std::string& index) const;

	// inheritance
	bool isOverride() const;
	bool isInherited() const;
	std::vector<ConfigNode*> getOverriden();
	std::vector<ConfigNode*> getOverriding();
	std::vector<ConfigNode*> getInherited();
	std::vector<ConfigNode*> getInheriting();
	ConfigNode* getDefinition();

	// cache
	bool isPopulated() const;
	bool isResolved() const;

	static bool populatedFilter(Model::ConfigNode& node);
	static bool inheritedFilter(Model::ConfigNode& node);
	static QVariant convertToVariant(const Value& value);
	static Value convertFromVariant(const QVariant& var);

protected:
	friend class ConfigModelAccessor;
	friend class ConfigModel;

	// mutators
	void setName(const std::string& newName);
	void setIndex(const std::string& newIndex);
	void setValue(const Value& value);
	void setUnit(const std::string& newUnit);
	void setRange(RangePtr range);
	void setIncludeUrl(const std::string& url);

	// alteration
	void insertAfter(Attribute* attribute);
	void insertInside(Attribute* attribute);
	std::auto_ptr<Attribute> remove();

	// cache flushing
	void clear();
	void recomputeChildren();
	void recomputeConfig();
	void updateConfigPath();

	PathNode generateNewAttributeName();

private:
	void populate() const;
	void resolve() const;

	void ensurePopulated() const;
	void ensureResolved() const;
};


///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::ConfigNode inline methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Returns the parent of this configuration node.
 * Returns 0 if it is the configuration root.
 */
inline ConfigNode* ConfigNode::getParent()
{
	return parent;
}

/**
 * Returns the parent of this configuration node.
 * Returns 0 if it is the configuration root.
 */
inline const ConfigNode* ConfigNode::getParent() const
{
	return parent;
}

/**
 * Returns the row number of this node in its parent.
 */
inline int ConfigNode::getRowInParent() const
{
	return rowInParent;
}

/**
 * Tells whether this node is the configuration root.
 */
inline bool ConfigNode::isRoot() const
{
	return getParent() == 0;
}

/**
 * Tells whether this node represents an attribute inherited from a default object.
 */
inline bool ConfigNode::isInherited() const
{
	return getBacktrackLevel() > 0;
}

inline bool ConfigNode::isPopulated() const
{
	return populated;
}

inline bool ConfigNode::isResolved() const
{
	return resolved;
}

}
}

#endif // CCFG_MODEL_CONFIGNODE_H_INCLUDED
