#include "AbstractConfigModelCommand.h"

namespace CCfg
{
namespace Model
{
namespace Command
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Command::AbstractConfigModelCommand class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 * @param accessor the object responsible for model access.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
AbstractConfigModelCommand::AbstractConfigModelCommand(ConfigModelAccessor* accessor,
	AbstractCommand* parent)
	: AbstractCommand(parent), accessor(accessor), succeeded(false)
{
}

/**
 * Destructor.
 */
AbstractConfigModelCommand::~AbstractConfigModelCommand()
{
}

/**
 * Returns the object responsible for model access.
 */
ConfigModelAccessor* AbstractConfigModelCommand::getAccessor()
{
	return accessor;
}

/**
 * Returns true if the command execution succeeded.
 */
bool AbstractConfigModelCommand::hasSucceeded() const
{
	return succeeded;
}

}
}
}
