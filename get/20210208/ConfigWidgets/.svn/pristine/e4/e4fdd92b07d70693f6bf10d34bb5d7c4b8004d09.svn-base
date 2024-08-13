#ifndef CCFG_MODEL_COMMAND_SETCONFIGNODERANGE_H_INCLUDED
#define CCFG_MODEL_COMMAND_SETCONFIGNODERANGE_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "CCfg/Path.h"
#include <memory>

namespace CCfg
{
// forward declarations
class Range;
typedef std::auto_ptr<Range> RangePtr;

namespace Model
{
// forward declarations
class ConfigNode;

namespace Command
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This command changes the range of a configuration node.
 *
 * @author fchateau
 */
class SetConfigNodeRange : public AbstractConfigModelCommand
{
	Path nodePath;	///< the path of the node to modify.
	RangePtr oldRange;	///< the range before command execution.
	RangePtr newRange;	///< the unit that the command will set in redo().

public:
	SetConfigNodeRange(ConfigNode* node, RangePtr newRange,
			ConfigModelAccessor* accessor,
			AbstractConfigModelCommand* parent = 0);
	~SetConfigNodeRange();

	int id() const;
	bool mergeWith(const QUndoCommand* command);

	void redo();
	void undo();

private:
	QString generateText();

	// non-default-constructible, non-copyable
	SetConfigNodeRange();
	SetConfigNodeRange(const SetConfigNodeRange& rvalue);
	SetConfigNodeRange& operator=(const SetConfigNodeRange& rvalue);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_SETCONFIGNODERANGE_H_INCLUDED
