#include "ConfigEditorDelegate.h"
#include "../../Model/ConfigModel.h"
#include "../../Model/ConfigNode.h"
#include "CCfg/Ui/TypeDisplayInfo.h"
#include "CCfg/Value.h"
#include "Utils/Ui/FileChooser.h"
#include "Utils/Error/ExcepHandling.h"


namespace CCfg
{
using Model::ConfigModel;
using Model::ConfigNode;
using Ui::TypeDisplayInfo;
using Utils::Ui::FileChooser;

namespace Control
{
namespace Delegate
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Control::Delegate::ConfigEditorDelegate class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent the parent QObject.
 */
ConfigEditorDelegate::ConfigEditorDelegate(QObject* parent)
	: AbstractDelegate(parent), fileDialogFilter()
{
}

/**
 * Destructor.
 */
ConfigEditorDelegate::~ConfigEditorDelegate()
{
}

/**
 * Returns a suitable widget for editing a model item.
 * @param parent the parent for the editor.
 * @param option unused.
 * @param index the index of the item to edit.
 * @return Returns the
 */
QWidget* ConfigEditorDelegate::createEditor(QWidget* parent,
		const QStyleOptionViewItem& option,
		const QModelIndex& index) const
{
	try
	{
		if(index.column() == ConfigModel::ValueColumn)
		{
			ConfigNode* const node = ConfigModel::indexToNode(index);
			assert(not node->isObject());
			return TypeDisplayInfo::createEditor(
					node->getValue().getType(),
					node->getRange(),
					parent);
		}
		else if(index.column() == ConfigModel::IncludeColumn)
		{
			// Note: we should spawn another editor if the URL is not a path.
			FileChooser* fc = new FileChooser(parent);
			fc->setSelectionMode(FileChooser::FileToSave);
			fc->setPathMode(FileChooser::Relative);
			fc->setFilter(fileDialogFilter);
			fc->setBaseDirectory(ConfigModel::getConfigDirectory(index));
			return fc;
		}
		else
		{
			return QItemDelegate::createEditor(parent, option, index);
		}
	}
	CATCH_AND_REPORT
	return 0;
}

}
}
}
