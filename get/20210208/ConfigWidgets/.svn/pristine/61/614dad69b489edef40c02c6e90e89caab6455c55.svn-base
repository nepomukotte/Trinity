#ifndef CCFG_MODEL_CONFIGMUTATOR_H_INCLUDED
#define CCFG_MODEL_CONFIGMUTATOR_H_INCLUDED

#include "AbstractConfigModelAccessor.h"
#include "Utils/Error/Exception.h"
#include "CCfg/TypeInfo.h"
#include <QObject>
#include <memory>

namespace CCfg
{
// forward declarations
struct PathNode;
class CConfig;

namespace Model
{
///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a facade for accessing the configuration model.
 * It provide high level methods to modify the model, while ensuring
 * its integrity (e.g: default object semantics). It sends signal to notify
 * higher levels about the nodes to clear or to recompute to keep the
 * model up-to-date and coherent.
 * This class owns a configuration, because independently of which node
 * is currently being processed, it is often necessary to access the whole
 * model.
 *
 * @author fchateau
 */
class ConfigModelAccessor : public QObject, public AbstractConfigModelAccessor
{
	Q_OBJECT

	/// the root of the configuration tree.
	std::auto_ptr<ConfigNode> rootNode;
	/// the document owning the low level configuration tree.
	std::auto_ptr<Io::Document> doc;

public:
	// @{ exception classes
	class InvalidAttributeName;
	class OperationOnInherited;
	class AttributeAlreadyExists;
	class ForbiddenTypeChange;
	// @}

	ConfigModelAccessor(QObject* parent = 0);
	virtual ~ConfigModelAccessor();

	// configuration root access
	void setDocument(std::auto_ptr<Io::Document> newDocument);
	const Io::Document* getDocument() const;
	Io::Document* getDocument();
	void setRootNode(std::auto_ptr<ConfigNode> newRoot);
	ConfigNode* getRootNode();

	// model accessors
	ConfigNode* pathToNode(const Path& path) const;


	// model mutators
	bool setNameAndIndex(ConfigNode* node, const std::string& newNameAndIndex);
	bool setName(ConfigNode* node, const std::string& newName);
	bool setIndex(ConfigNode* node, const std::string& newIndex);
	bool setValue(ConfigNode* node, const Value& value);
	bool setUnit(ConfigNode* node, const std::string& unit);
	bool setRange(ConfigNode* node, RangePtr range);
	bool setIncludeUrl(ConfigNode* node, const std::string& url);
	bool insertInside(ConfigNode* parent, std::auto_ptr<Attribute> newAttribute);
	bool insertAfter(ConfigNode* node, std::auto_ptr<Attribute> newAttribute);
	bool override(ConfigNode*& node, std::vector<ConfigNode*>& createdNodes);
	bool remove(ConfigNode* parent, const IndexList& indexes, Object* removed = 0);
	bool remove(ConfigNode* parent, unsigned index, Object* removed = 0);
	bool copy(ConfigNode* source, ConfigNode* dest);
	bool move(ConfigNode* source, const IndexList& srcIndexes, ConfigNode* dest, unsigned insertionIndex);
	bool include(ConfigNode* node, const std::string& url);
	bool import(ConfigNode* node, const std::string& url);

	void insert(ConfigNode* parent, Object* attributes, const IndexList& insertedNodeIndex);

signals:
	// mutation notifications
	void beginInsertion(ConfigNode* parent);
	void endInsertion(ConfigNode* parent, const IndexList& indexes);
	void beginRemoval(ConfigNode* parent, const IndexList& indexes);
	void endRemoval(ConfigNode* parent);

	// general update notification
	void beginClear(ConfigNode* node);
	void endClear();
	void beginRecompute(ConfigNode* node);
	void endRecompute();
	void beginUpdate();
	void endUpdate();

protected:
	IndexList computePostInsertionIndexes(unsigned elementCount, const IndexList& insertedElementIndex) const;
	std::vector<ConfigNode> computeInserterNodes(ConfigNode* parent, const IndexList& existingNodeIndex, const IndexList& insertedNodeIndex, IndexList& inserterNodeIndex) const;

	class Update;
	friend class Update;

	void clear(ConfigNode* node);
	void recompute(ConfigNode* node);
	void emitBeginUpdate();
	void emitEndUpdate();

	bool warnIfAlreadyExist(ConfigNode* parent, const PathNode& newNodeName);
	bool warnIfAlreadyExist(const Attribute* parent, const PathNode& newNodeName);
	bool warnIfInherited(ConfigNode* parent);
};


/**
 * This class is used to instanciate objects whose lifetime delimitate a model
 * update.
 * When updating a configuration node, there are side effects at several
 * locations of the model which requires to recompute things.
 * Instances of this class carry out necessary operations before and
 * after the alteration on the low level model, in their constructor and
 * destructor, respectively.
 */
class ConfigModelAccessor::Update
{
	/// the object through which modifications on the model are carried out.
	ConfigModelAccessor* accessor;
	/// the parent of the updated node.
	ConfigNode* parent;
	/// the list of impacted nodes (those who inherits from the parent).
	std::vector<ConfigNode*> impactedNodes;

public:
	Update(ConfigModelAccessor* accessor, ConfigNode* parent);
	~Update();
};

/**
 *
 */
class ConfigModelAccessor::InvalidAttributeName : public Utils::Error::Exception
{
	std::string fullName; ///< the invalid name and index
public:
	InvalidAttributeName(const std::string& fullName);
	~InvalidAttributeName() throw() {}

	const std::string& getNameIndex() const { return fullName; }

protected:
	std::string formatTitle() const;
	std::string formatBrief() const;
	std::string formatSuggestion() const;
};

/**
 *
 */
class ConfigModelAccessor::OperationOnInherited : public Utils::Error::Exception
{
public:
	OperationOnInherited() {}
	~OperationOnInherited() throw() {}

protected:
	std::string formatTitle() const;
	std::string formatBrief() const;
};

/**
 *
 */
class ConfigModelAccessor::AttributeAlreadyExists : public Utils::Error::Exception
{
	std::string attributeName;
public:
	AttributeAlreadyExists(const std::string& name);
	~AttributeAlreadyExists() throw() {}

	const std::string& getAttributeName() const { return attributeName; }

protected:
	std::string formatTitle() const;
	std::string formatBrief() const;
	std::string formatSuggestion() const;
};

/**
 *
 */
class ConfigModelAccessor::ForbiddenTypeChange : public Utils::Error::Exception
{
	ValueTypeId oldType;
	ValueTypeId newType;
public:
	ForbiddenTypeChange(ValueTypeId oldType, ValueTypeId newType);
	~ForbiddenTypeChange() throw() {}

protected:
	std::string formatTitle() const;
	std::string formatBrief() const;
	std::string formatDetails() const;
};

}
}

#endif // CCFG_MODEL_CONFIGMUTATOR_H_INCLUDED
