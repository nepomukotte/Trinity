#ifndef CCFG_MODEL_COMMAND_SETCONFIGNODEINCLUDEURL_H_INCLUDED
#define CCFG_MODEL_COMMAND_SETCONFIGNODEINCLUDEURL_H_INCLUDED

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
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This command changes the include url of a configuration node.
 *
 * @author fchateau
 */
class SetConfigNodeIncludeUrl : public AbstractConfigModelCommand
{
	Path nodePath;	///< the path of the node to modify.
	std::string oldUrl;	///< the include URL before command execution.
	std::string newUrl;	///< the include URL that the command will set in redo().

public:
	SetConfigNodeIncludeUrl(ConfigNode* node, const std::string& newUrl,
			ConfigModelAccessor* accessor,
			AbstractConfigModelCommand* parent = 0);
	~SetConfigNodeIncludeUrl();

	int id() const;
	bool mergeWith(const QUndoCommand* command);

	void redo();
	void undo();

private:
	QString generateText();

	// non-default-constructible, non-copyable
	SetConfigNodeIncludeUrl();
	SetConfigNodeIncludeUrl(const SetConfigNodeIncludeUrl& rvalue);
	SetConfigNodeIncludeUrl& operator=(const SetConfigNodeIncludeUrl& rvalue);
};

}
}
}

#endif // CCFG_MODEL_COMMAND_SETCONFIGNODEINCLUDEURL_H_INCLUDED
