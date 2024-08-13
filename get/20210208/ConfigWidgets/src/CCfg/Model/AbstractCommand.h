#ifndef CCFG_MODEL_ABSTRACTCOMMAND_H_INCLUDED
#define CCFG_MODEL_ABSTRACTCOMMAND_H_INCLUDED

#include <QUndoCommand>

namespace CCfg
{
namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This interface must be realized by classes modifying a model.
 * Every model alteration should be implemented by creating a derived class
 * carrying out the operation. This requirement is mandatory to achieve correct
 * undo/redo semantic.
 * See QUndoCommand to find the list of methods to implement.
 *
 * @author fchateau
 */
class AbstractCommand : public QUndoCommand
{
public:
	AbstractCommand(AbstractCommand* parent = 0);
	virtual ~AbstractCommand();
};

}

}

#endif // CCFG_MODEL_ABSTRACTCOMMAND_H_INCLUDED
