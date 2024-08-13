#ifndef CCFG_MODEL_COMMAND_ABSTRACTCONFIGMODELCOMMAND_H_INCLUDED
#define CCFG_MODEL_COMMAND_ABSTRACTCONFIGMODELCOMMAND_H_INCLUDED

#include "../AbstractCommand.h"

namespace CCfg
{
namespace Model
{
// forward declarations
class ConfigModelAccessor;
class ConfigModel;

namespace Command
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This interface must be realized by classes modifying a configuration model.
 * See AbstractCommand to find the list of methods to implement.
 *
 * @author fchateau
 */
class AbstractConfigModelCommand : public AbstractCommand
{
protected:
	ConfigModelAccessor* accessor; ///< the object responsible for model access.
	bool succeeded;	///< indicates whether the last call to redo() succeeded.

public:
	AbstractConfigModelCommand(ConfigModelAccessor* accessor, AbstractCommand* parent = 0);
	virtual ~AbstractConfigModelCommand();

	ConfigModelAccessor* getAccessor();
	bool hasSucceeded() const;
};

}
}
}

#endif // CCFG_MODEL_COMMAND_ABSTRACTCONFIGMODELCOMMAND_H_INCLUDED
