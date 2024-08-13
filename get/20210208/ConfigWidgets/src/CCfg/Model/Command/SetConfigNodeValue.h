#ifndef CCFG_MODEL_COMMAND_SETCONFIGNODEVALUE_H_INCLUDED
#define CCFG_MODEL_COMMAND_SETCONFIGNODEVALUE_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "CCfg/Path.h"
#include "CCfg/Value.h"

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
 * This command changes the value of a configuration node.
 *
 * @author fchateau
 */
class SetConfigNodeValue : public AbstractConfigModelCommand
{
	Path nodePath;	///< the path of the node to modify.
	Value oldValue;	///< the value before command execution.
	Value newValue;	///< the value that the command will set in redo().
public:
	SetConfigNodeValue(ConfigNode* node, const Value& newValue,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent = 0);
	~SetConfigNodeValue();

	int id() const;
	bool mergeWith(const QUndoCommand* command);

	void redo();
	void undo();

private:
	QString generateText();

	// non-default-constructible, non-copyable
	SetConfigNodeValue();
	SetConfigNodeValue(const SetConfigNodeValue& cmd);
	SetConfigNodeValue& operator=(const SetConfigNodeValue& cmd);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_SETCONFIGNODEVALUE_H_INCLUDED
