#ifndef CCFG_UI_VALUEEDITIONWIDGET_H_INCLUDED
#define CCFG_UI_VALUEEDITIONWIDGET_H_INCLUDED

#include "CCfg/TypeInfo.h"
#include <QWidget>
#include <QString>

// forward declarations
class QStackedLayout;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This widget allows to edit a configuration value whose type can change
 * dynamically while providing a specialized editor for each type.
 *
 * @author fchateau
 */
class ValueEditionWidget : public QWidget
{
	Q_OBJECT

	QStackedLayout* layout;	///< the layout where the different editors are stacked.
	ValueTypeId type;	///< the current type of the configuration value.

public:
	ValueEditionWidget(QWidget* parent = 0);
	virtual ~ValueEditionWidget();

	ValueTypeId getEditedType() const;
	void setEditedType(ValueTypeId newType);

	QString getValue() const;
	void setValue(const QString& newValue);

signals:
	void valueChanged(const QString& newValue);

protected:
	void disconnectEditor(QWidget* ed);
	void connectEditor(QWidget* ed);

protected slots:
	void onTextChanged(const QString& text);
	void onValueChanged(int intValue);
	void onCurrentIndexChanged(const QString& text);

private:
	ValueEditionWidget(const ValueEditionWidget& right);
	ValueEditionWidget& operator=(const ValueEditionWidget& right);
};

}
}

#endif // CCFG_UI_VALUEEDITIONWIDGET_H_INCLUDED
