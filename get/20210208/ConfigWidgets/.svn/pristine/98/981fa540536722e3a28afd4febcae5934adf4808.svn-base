#ifndef CCFG_UI_ABSTRACTDELEGATE_H_INCLUDED
#define CCFG_UI_ABSTRACTDELEGATE_H_INCLUDED

#include <QItemDelegate>

namespace CCfg
{
namespace Control
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for model/view delegates of the application.
 *
 * @author fchateau
 */
class AbstractDelegate : public QItemDelegate
{
public:
	AbstractDelegate(QObject* parent = 0);
	virtual ~AbstractDelegate();

	virtual void updateEditorGeometry(QWidget* editor, const QStyleOptionViewItem& option, const QModelIndex& index) const;
};

}
}

#endif // CCFG_UI_ABSTRACTDELEGATE_H_INCLUDED
