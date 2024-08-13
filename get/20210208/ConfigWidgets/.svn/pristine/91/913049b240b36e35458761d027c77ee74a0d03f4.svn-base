#ifndef CCFG_MODEL_UNDOABLECONFIGMUTATOR_H_INCLUDED
#define CCFG_MODEL_UNDOABLECONFIGMUTATOR_H_INCLUDED

#include "AbstractConfigModelAccessor.h"
#include <QObject>

// forward declarations
class QUndoStack;

namespace CCfg
{
namespace Model
{
// forward declarations
class ConfigModelAccessor;
namespace Command { class AbstractConfigModelCommand; }

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

class CommandBasedConfigModelAccessor : public AbstractConfigModelAccessor
{
public:
	explicit CommandBasedConfigModelAccessor(ConfigModelAccessor* accessor = 0);
	virtual ~CommandBasedConfigModelAccessor();

	void setAccessor(ConfigModelAccessor* newAccessor) { accessor = newAccessor; }

	// configuration root access
	virtual void setDocument(std::auto_ptr<Io::Document> newDocument);
	virtual const Io::Document* getDocument() const;
	virtual Io::Document* getDocument();
	virtual void setRootNode(std::auto_ptr<ConfigNode> newRoot);
	virtual ConfigNode* getRootNode();

	// model accessors
	virtual ConfigNode* pathToNode(const Path& path) const;

	// model mutators
	virtual bool setNameAndIndex(ConfigNode* node, const std::string& newNameAndIndex);
	virtual bool setName(ConfigNode* node, const std::string& newName);
	virtual bool setIndex(ConfigNode* node, const std::string& newIndex);
	virtual bool setValue(ConfigNode* node, const Value& value);
	virtual bool setUnit(ConfigNode* node, const std::string& unit);
	virtual bool setRange(ConfigNode* node, RangePtr range);
	virtual bool setIncludeUrl(ConfigNode* node, const std::string& url);
	virtual bool insertInside(ConfigNode* parent, std::auto_ptr<Attribute> newAttribute);
	virtual bool insertAfter(ConfigNode* node, std::auto_ptr<Attribute> newAttribute);
	virtual bool override(ConfigNode*& node, std::vector<ConfigNode*>& createdNodes);
	virtual bool remove(ConfigNode* parent, const IndexList& indexes, Object* removed = 0);
	virtual bool remove(ConfigNode* parent, unsigned index, Object* removed = 0);
	virtual bool copy(ConfigNode* source, ConfigNode* dest);
	virtual bool move(ConfigNode* source, const IndexList& srcIndexes, ConfigNode* dest, unsigned insertionIndex);
	virtual bool include(ConfigNode* node, const std::string& url);
	virtual bool import(ConfigNode* node, const std::string& url);

protected:
	virtual bool execute(Command::AbstractConfigModelCommand* cmd) = 0;

	ConfigModelAccessor* accessor;
};

/**
 * This class is a facade for accessing the configuration model.
 * It provides undo features to make high level commands reversible.
 *
 * @see ConfigModelAccessor
 * @author fchateau
 */
class UndoableConfigModelAccessor : public QObject, public CommandBasedConfigModelAccessor
{
	Q_OBJECT

	/// the stack retaining every command that altered the model, making them undoable.
	QUndoStack* undoStack;

public:
	UndoableConfigModelAccessor(QObject* parent = 0);
	virtual ~UndoableConfigModelAccessor();

	QUndoStack* getUndoStack();
	ConfigModelAccessor* getConfigModelAccessor();

	void setDocument(std::auto_ptr<Io::Document> newDocument);
	void setRootNode(std::auto_ptr<ConfigNode> newRoot);

protected:
	bool execute(Command::AbstractConfigModelCommand* cmd);
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::UndoableConfigModelAccessor inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 *
 */
inline QUndoStack* UndoableConfigModelAccessor::getUndoStack()
{
	return undoStack;
}

/**
 *
 */
inline ConfigModelAccessor* UndoableConfigModelAccessor::getConfigModelAccessor()
{
	return accessor;
}

}
}

#endif // CCFG_MODEL_UNDOABLECONFIGMUTATOR_H_INCLUDED
