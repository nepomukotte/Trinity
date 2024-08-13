#ifndef CCFG_CONTROL_DELEGATE_CONFIGEDITORDELEGATE_H_INCLUDED
#define CCFG_CONTROL_DELEGATE_CONFIGEDITORDELEGATE_H_INCLUDED

#include "../AbstractDelegate.h"
#include "CCfg/TypeInfo.h"

// forward declarations
class QDir;


namespace CCfg
{
namespace Control
{
namespace Delegate
{
///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This delegate enables edition of configuration properties, whatever their type.
 *
 * @author fchateau
 */
class ConfigEditorDelegate : public AbstractDelegate
{
public:
	ConfigEditorDelegate(QObject* parent = 0);
	~ConfigEditorDelegate();

	QWidget* createEditor(QWidget* parent, const QStyleOptionViewItem& option, const QModelIndex& index) const;

	const QString& getFileDialogFilter() const;
	void setFileDialogFilter(const QString& filter);

private:
	QString fileDialogFilter; ///< the extension filter for file dialogs created by this delegate.
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Control::Delegate::ConfigEditorDelegate inline methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Returns the extension filter for file dialogs created by this delegate.
 */
inline const QString& ConfigEditorDelegate::getFileDialogFilter() const
{
	return fileDialogFilter;
}

/**
 * Defines the extension filter for file dialogs created by this delegate.
 */
inline void ConfigEditorDelegate::setFileDialogFilter(const QString& filter)
{
	fileDialogFilter = filter;
}

}
}
}

#endif // CCFG_CONTROL_DELEGATE_CONFIGEDITORDELEGATE_H_INCLUDED
