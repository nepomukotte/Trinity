#ifndef CCFG_MODEL_COMMAND_INSERTCONFIGNODE_H_INCLUDED
#define CCFG_MODEL_COMMAND_INSERTCONFIGNODE_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "CCfg/Path.h"
#include <memory>

namespace CCfg
{
// forward declarations
struct PathNode;
class Attribute;

namespace Model
{
// forward declarations
class ConfigNode;

namespace Command
{
///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This command inserts a new configuration node.
 *
 * @author fchateau
 */
class InsertConfigNode : public AbstractConfigModelCommand
{
public:
	enum Position {Child, Sibling};

private:
	Path referencePath;	///< the path of an existing node that serves as a landmark for insertion.
	Position position;	///< the position of insertion with respect to the reference node.
	std::auto_ptr<Attribute> newAttribute;	///< the attribute to insert.

public:
	InsertConfigNode(ConfigNode* referenceNode,
		Position position,
		std::auto_ptr<Attribute> attribute,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent = 0);
	~InsertConfigNode();

	int id() const;

	void redo();
	void undo();

protected:
	QString generateText();

private:
	// non-default-constructible, non-copyable
	InsertConfigNode();
	InsertConfigNode(const InsertConfigNode& cmd);
	InsertConfigNode& operator=(const InsertConfigNode& cmd);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_INSERTCONFIGNODE_H_INCLUDED
