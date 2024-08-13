#ifndef CCFG_MODEL_COMMAND_DELETECONFIGNODE_H_INCLUDED
#define CCFG_MODEL_COMMAND_DELETECONFIGNODE_H_INCLUDED

#include "AbstractConfigModelCommand.h"
#include "../AbstractConfigModelAccessor.h"
#include "CCfg/Path.h"
#include <vector>

namespace CCfg
{
// forward declarations
class Object;

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
 *
 * @author fchateau
 */
class DeleteConfigNode : public AbstractConfigModelCommand
{
	Path parentPath;	///< the parent of the ConfigNodes to remove.
	IndexList indexes;	///< the sorted list of row indexes of the ConfigNodes to remove.
	/// the Object that contains and owns removed attributes, after redo. It is
	/// empty after undo.
	std::auto_ptr<Object> deletedContainer;

public:
	DeleteConfigNode(ConfigNode* parentNode,
		const IndexList& indexes,
		ConfigModelAccessor* accessor,
		AbstractConfigModelCommand* parent = 0);
	~DeleteConfigNode();

	int id() const;

	void redo();
	void undo();

private:
	QString generateText();

	DeleteConfigNode();
	DeleteConfigNode(const DeleteConfigNode& right);
	DeleteConfigNode& operator=(const DeleteConfigNode& right);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_DELETECONFIGNODE_H_INCLUDED
