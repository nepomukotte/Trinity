#ifndef CCFG_MODEL_COMMAND_SETCONFIGNODEUNIT_H_INCLUDED
#define CCFG_MODEL_COMMAND_SETCONFIGNODEUNIT_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "CCfg/Path.h"
#include <string>

namespace CCfg
{
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
 * This command changes the unit of a configuration node.
 *
 * @author fchateau.
 */
class SetConfigNodeUnit : public AbstractConfigModelCommand
{
	Path nodePath;	///< the path of the node to modify.
	std::string oldUnit;	///< the unit before command execution.
	std::string newUnit;	///< the unit that the command will set in redo().
public:
	SetConfigNodeUnit(ConfigNode* node, const std::string& newUnit,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent = 0);
	~SetConfigNodeUnit();

	int id() const;
	bool mergeWith(const QUndoCommand* command);

	void redo();
	void undo();

private:
	QString generateText();

	// non-default-constructible, non-copyable
	SetConfigNodeUnit();
	SetConfigNodeUnit(const SetConfigNodeUnit& cmd);
	SetConfigNodeUnit& operator=(const SetConfigNodeUnit& cmd);
};

}

}

}

#endif // CCFG_MODEL_COMMAND_SETCONFIGNODEUNIT_H_INCLUDED
