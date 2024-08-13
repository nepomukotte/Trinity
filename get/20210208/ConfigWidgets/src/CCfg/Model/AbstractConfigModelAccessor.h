#ifndef CCFG_MODEL_ABSTRACTCONFIGMODELACCESSOR_H_INCLUDED
#define CCFG_MODEL_ABSTRACTCONFIGMODELACCESSOR_H_INCLUDED

#include <string>
#include <vector>
#include <memory>


namespace CCfg
{
// forward declarations
class Object;
class Attribute;
class Value;
class Path;
class Range;
typedef std::auto_ptr<Range> RangePtr;
namespace Io { class Document; }

namespace Model
{
// forward declarations
class ConfigNode;

// typedefs
typedef std::vector<unsigned> IndexList;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This interface defines high level methods to read and modify a configuration
 * model.
 * @author fchateau
 */
class AbstractConfigModelAccessor
{
public:
	virtual ~AbstractConfigModelAccessor();

	// configuration root access
	virtual void setDocument(std::auto_ptr<Io::Document> newDocument) = 0;
	virtual const Io::Document* getDocument() const = 0;
	virtual Io::Document* getDocument() = 0;
	virtual void setRootNode(std::auto_ptr<ConfigNode> newRoot) = 0;
	virtual ConfigNode* getRootNode() = 0;

	// model accessors
	virtual ConfigNode* pathToNode(const Path& path) const = 0;

	// model mutators
	virtual bool setNameAndIndex(ConfigNode* node, const std::string& newNameAndIndex) = 0;
	virtual bool setName(ConfigNode* node, const std::string& newName) = 0;
	virtual bool setIndex(ConfigNode* node, const std::string& newIndex) = 0;
	virtual bool setValue(ConfigNode* node, const Value& value) = 0;
	virtual bool setUnit(ConfigNode* node, const std::string& unit) = 0;
	virtual bool setRange(ConfigNode* node, RangePtr range) = 0;
	virtual bool setIncludeUrl(ConfigNode* node, const std::string& url) = 0;
	virtual bool insertInside(ConfigNode* parent, std::auto_ptr<Attribute> newAttribute) = 0;
	virtual bool insertAfter(ConfigNode* node, std::auto_ptr<Attribute> newAttribute) = 0;
	virtual bool override(ConfigNode*& node, std::vector<ConfigNode*>& createdNodes) = 0;
	virtual bool remove(ConfigNode* parent, const IndexList& indexes, Object* removed = 0) = 0;
	virtual bool remove(ConfigNode* parent, unsigned index, Object* removed = 0) = 0;
	virtual bool copy(ConfigNode* source, ConfigNode* dest) = 0;
	virtual bool move(ConfigNode* source, const IndexList& srcIndexes, ConfigNode* dest, unsigned insertionIndex) = 0;
	virtual bool include(ConfigNode* node, const std::string& url) = 0;
	virtual bool import(ConfigNode* node, const std::string& url) = 0;
};

}
}

#endif // CCFG_MODEL_ABSTRACTCONFIGMODELACCESSOR_H_INCLUDED
