#include "AbstractDelegate.h"

namespace CCfg
{
namespace Control
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Control::AbstractDelegate class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent the parent QObject.
 */
AbstractDelegate::AbstractDelegate(QObject* parent)
	: QItemDelegate(parent)
{
}

/**
 * Destructor.
 */
AbstractDelegate::~AbstractDelegate()
{
}

/**
 * Updates the editor for the item specified by index according to the style option given.
 * @param editor the editor to update.
 * @param option the style options.
 * @param index the index of the item being edited.
 */
void AbstractDelegate::updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& /*index*/) const
{
	editor->setGeometry(option.rect);
}

}
}
