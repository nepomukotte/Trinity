#include "ValueEditionWidget.h"
#include "TypeDisplayInfo.h"
#include <QStackedLayout>
#include <QLineEdit>
#include <QSpinBox>
#include <QComboBox>
#include <QLabel>
#include <cassert>

namespace CCfg
{
namespace Ui
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::ValueEditionWidget class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 * @param parent the parent widget.
 */
ValueEditionWidget::ValueEditionWidget(QWidget* parent)
	: QWidget(parent), layout(new QStackedLayout(this)),
	type(INTEGER)
{
    layout->setSpacing(0);
    layout->setMargin(0);

    for(int i = 0; i < VALUETYPE_COUNT; ++i)
    {
    	layout->addWidget(TypeDisplayInfo::createEditor(static_cast<ValueTypeId>(i), 0, this));
    }
    setFocusProxy(layout->widget(0));
    setEditedType(STRING);
}

/**
 * Destructor.
 */
ValueEditionWidget::~ValueEditionWidget()
{
}

/**
 * Returns the current type of the edited configuration value.
 */
ValueTypeId ValueEditionWidget::getEditedType() const
{
	return type;
}

/**
 * Defines the type of the edited configuration value.
 * It immediately causes the edition widget to switch to the fitted editor.
 */
void ValueEditionWidget::setEditedType(ValueTypeId newType)
{
	if(newType != type)
	{
		disconnectEditor(layout->currentWidget());

		type = newType;
		layout->setCurrentIndex(newType);
		emit valueChanged(getValue());

		connectEditor(layout->currentWidget());
	}
}

/**
 * Disconnect all slots from the current editor.
 */
void ValueEditionWidget::disconnectEditor(QWidget* ed)
{
	if(ed != 0)
	{
		disconnect(ed, 0, this, 0);
	}
}

/**
 * Connect necessary slots to the current editor.
 */
void ValueEditionWidget::connectEditor(QWidget* ed)
{
	if(qobject_cast<QLineEdit*>(ed))
	{
		connect(ed, SIGNAL(textChanged(QString)),
			this, SLOT(onTextChanged(QString)));
	}
	else if(qobject_cast<QSpinBox*>(ed))
	{
		connect(ed, SIGNAL(valueChanged(int)),
			this, SLOT(onValueChanged(int)));
	}
	else if(qobject_cast<QComboBox*>(ed))
	{
		connect(ed, SIGNAL(currentIndexChanged(const QString &)),
			this, SLOT(onCurrentIndexChanged(const QString &)));
	}
	else
	{
		assert(false);
	}
}

/**
 * Returns the current edited value.
 */
QString ValueEditionWidget::getValue() const
{
	return layout->currentWidget()->property(
		TypeDisplayInfo::getEditorPropertyName(type)).toString();
}

/**
 * Defines the current edited value.
 */
void ValueEditionWidget::setValue(const QString& newValue)
{
	if(getValue() != newValue)
	{
		layout->currentWidget()->setProperty(
			TypeDisplayInfo::getEditorPropertyName(type), newValue);
		emit valueChanged(newValue);
	}
}

/**
 * Executes necessary operations when edition in a "line edit" is finished.
 */
void ValueEditionWidget::onTextChanged(const QString& text)
{
	assert(qobject_cast<QLineEdit*>(layout->currentWidget()) != 0);
	emit valueChanged(text);
}

/**
 * Executes necessary operations when the value of a "spin box" changed.
 */
void ValueEditionWidget::onValueChanged(int intValue)
{
	assert(qobject_cast<QSpinBox*>(layout->currentWidget()) != 0);
	emit valueChanged(QString::number(intValue));
}

/**
 * Executes necessary operations when the current item of a "combo box" changed.
 */
void ValueEditionWidget::onCurrentIndexChanged(const QString& text)
{
	assert(qobject_cast<QComboBox*>(layout->currentWidget()) != 0);
	emit valueChanged(text);
}

}
}
