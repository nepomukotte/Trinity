#include "AbstractCommand.h"

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::AbstractCommand class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent if not null, this is a command in which this one
 * will be nested, thus becoming part of a composite command.
 */
AbstractCommand::AbstractCommand(AbstractCommand* parent)
	: QUndoCommand(parent)
{
}

/**
 * Destructor.
 */
AbstractCommand::~AbstractCommand()
{
}

}
}
