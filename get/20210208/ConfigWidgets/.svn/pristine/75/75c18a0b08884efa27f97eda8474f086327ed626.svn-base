#include "AttributeBuilderWidget.h"
#include "CCfg/Ui/AttributeBuilder.ui.h"
#include "../Model/AttributeBuilderModel.h"
#include "TypeDisplayInfo.h"
#include "Utils/Error/ExcepHandling.h"
#include <QMessageBox>
#include <QIntValidator>
#include <QRegExp>
#include <limits>
#include <cassert>

namespace CCfg
{
namespace Ui
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::AttributeBuilderWidget class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
AttributeBuilderWidget::AttributeBuilderWidget(QWidget* parent)
	: Base(parent),
	model(0), validTarget(false)
{
	setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Fixed);
	connect(ui_->lineEditName, SIGNAL(textChanged(QString)),
		this, SLOT(onNameChanged(QString)));
	connect(ui_->lineEditIndex, SIGNAL(textChanged(QString)),
		this, SLOT(onIndexChanged(QString)));
	connect(ui_->pushButtonDefault, SIGNAL(toggled(bool)),
		this, SLOT(onDefaultToggled(bool)));
	connect(ui_->comboBoxType, SIGNAL(currentIndexChanged(int)),
		this, SLOT(onTypeChanged(int)));
	connect(ui_->valueEditor, SIGNAL(valueChanged(QString)),
		this, SLOT(onValueChanged(QString)));
	connect(ui_->buttonInsert, SIGNAL(clicked()),
		this, SLOT(onInsertClicked()));
}

/**
 * Destructor.
 */
AttributeBuilderWidget::~AttributeBuilderWidget()
{
}

/**
 * Changes the model of this attribute builder.
 * The content of every widget is updated.
 * @param newModel the new data used to initialize inner widgets.
 */
void AttributeBuilderWidget::setModel(Model::AttributeBuilderModel* newModel)
{
	model = newModel;
	loadModel();
}

/**
 * Returns the model of this attribute builder.
 * Its values reflect those in the widgets.
 */
Model::AttributeBuilderModel* AttributeBuilderWidget::getModel()
{
	if(model == 0)
	{
		model = new Model::AttributeBuilderModel(this);
	}
	return model;
}

/**
 * Loads the model data into inner widgets.
 */
void AttributeBuilderWidget::loadModel()
{
	// we use the accessor even inside the class, because it creates a default
	// model if "model" is null.
	Model::AttributeBuilderModel* mdl = getModel();
	ui_->comboBoxType->setCurrentIndex(mdl->getType());
	ui_->lineEditName->setText(mdl->getName());
	if(mdl->getIndex() == "*")
	{
		ui_->pushButtonDefault->setChecked(true);
	}
	else
	{
		ui_->lineEditIndex->setText(mdl->getIndex());
		ui_->pushButtonDefault->setChecked(false);
	}
	ui_->valueEditor->setValue(mdl->getValue());
	ui_->buttonInsert->setEnabled(isInsertionReady());
}

/**
 * Updates the model data using inner widgets data.
 */
void AttributeBuilderWidget::refreshModel()
{
	Model::AttributeBuilderModel* mdl = getModel();
	mdl->setType((AttributeTypeId)ui_->comboBoxType->currentIndex());
	mdl->setName(ui_->lineEditName->text());
	if(ui_->pushButtonDefault->isChecked())
	{
		mdl->setIndex("*");
	}
	else
	{
		mdl->setIndex(ui_->lineEditIndex->text());
	}
	mdl->setValue(ui_->valueEditor->getValue());
}

/**
 * Enables or disables the "default" button.
 */
void AttributeBuilderWidget::enableDefault(bool isEnabled)
{
	bool wasEnabled = ui_->pushButtonDefault->isEnabled();
	ui_->pushButtonDefault->setEnabled(isEnabled);
	if(wasEnabled xor isEnabled)
	{
		commitIndex();
	}
}

/**
 * Enables or disables index edition.
 */
void AttributeBuilderWidget::enableIndex(bool isEnabled)
{
	bool wasEnabled = ui_->lineEditIndex->isEnabled();
	ui_->lineEditIndex->setEnabled(isEnabled);
	if(wasEnabled xor isEnabled)
	{
		commitIndex();
	}
}

/**
 * Enables or disables value edition.
 */
void AttributeBuilderWidget::enableValue(bool isEnabled)
{
	bool wasEnabled = ui_->valueEditor->isEnabled();
	ui_->valueEditor->setEnabled(isEnabled);
	if(wasEnabled xor isEnabled)
	{
		commitValue();
	}
}

/**
 * Commits the name from the editor to the model.
 */
void AttributeBuilderWidget::commitName()
{
	Model::AttributeBuilderModel* mdl = getModel();
	mdl->setName(ui_->lineEditName->text());	// update the model.
	ui_->buttonInsert->setEnabled(isInsertionReady());
}

/**
 * Commits the index from the editor to the model.
 */
void AttributeBuilderWidget::commitIndex()
{
	Model::AttributeBuilderModel* mdl = getModel();
	if(ui_->pushButtonDefault->isEnabled() and ui_->pushButtonDefault->isChecked())
	{
		mdl->setIndex("*");
	}
	else
	{
		mdl->setIndex(ui_->lineEditIndex->text());
	}
}

/**
 * Commits the type from the editor to the model.
 */
void AttributeBuilderWidget::commitType(AttributeTypeId newType)
{
	Model::AttributeBuilderModel* mdl = getModel();
	mdl->setType(newType);
	commitValue();
}

/**
 * Commits the value from the editor to the model.
 */
void AttributeBuilderWidget::commitValue()
{
	Model::AttributeBuilderModel* mdl = getModel();
	mdl->setValue(ui_->valueEditor->isEnabled() ?
		ui_->valueEditor->getValue() : "");
	ui_->buttonInsert->setEnabled(isInsertionReady());
}

/**
 * Determines if there are enough data for insertion.
 */
bool AttributeBuilderWidget::isInsertionReady()
{
	Model::AttributeBuilderModel* mdl = getModel();
	return validTarget and
		not mdl->getName().isEmpty() and
		(mdl->getType() == OBJECT or
				mdl->getType() == STRING_SCALAR or
				not mdl->getValue().isEmpty());
}

/**
 * Checks whether the form data are correct and issue a warning dialog otherwise.
 * @return Return true if the data are correct; false otherwise.
 */
bool AttributeBuilderWidget::checkData()
{
	if(model->getType() == OBJECT)
	{
		return true;
	}

	ValueTypeId type = TypeInfo::lookup(model->getType());
	const char* const regexpStr = TypeDisplayInfo::getValidationRegexp(type);
	if(regexpStr != 0 and regexpStr[0] != '\0')
	{
		QRegExp regexp(regexpStr);
		if(not regexp.exactMatch(model->getValue()))
		{
			QMessageBox::warning(this, tr("Bad input data"),
				tr("The value \"%1\" is not a valid %2")
					.arg(model->getValue())
					.arg(TypeInfo::getDisplayName(type)));
			return false;
		}
	}
	return true;
}

/**
 * Increments the index if it is a numeric (integer) index.
 */
void AttributeBuilderWidget::incrementNumericIndex()
{
	QIntValidator validator(std::numeric_limits<Int>::min(), std::numeric_limits<Int>::max(), 0);
	QString index(model->getIndex());
	int pos = 0;

	if(validator.validate(index, pos) == QValidator::Acceptable)
	{
		index.setNum(index.toInt() + 1);
		ui_->lineEditIndex->setText(index);
	}
}

/**
 * Reacts to a switch of insertion target.
 * The insertion target is the attribute where the new attribute will be
 * inserted (either as a sibling or as a child).
 * @param path the path of the insertion target.
 * @param isObject indicates whether the insertion target is an object or not.
 * (a scalar cannot have children, so we have to disable the radio button if false).
 * @param isRoot indicates whether the insertion target is the root or not.
 * (a root cannot have siblings, so the corresponding radio button must be disabled too).
 * @param isInherited indicates whether the insertion target is inherited or not.
 * (an inherited node cannot have children, so we have to disable the radio button if true).
 */
void AttributeBuilderWidget::onTargetChanged(const QString& path, bool isObject,
		bool isRoot, bool isInherited)
{
	try
	{
		validTarget = not isInherited;
		ui_->lineEditPath->setText(path);
		ui_->lineEditPath->end(false);
		ui_->radioButtonSibling->setEnabled(not isRoot);
		ui_->radioButtonChild->setEnabled(isObject);
		ui_->radioButtonSibling->setChecked(not isRoot);
		ui_->radioButtonChild->setChecked(isRoot);
		ui_->buttonInsert->setEnabled(isInsertionReady());
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a user modification on the name.
 */
void AttributeBuilderWidget::onNameChanged(const QString& /*newName*/)
{
	try
	{
		commitName();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a user modification on the index.
 */
void AttributeBuilderWidget::onIndexChanged(const QString& /*newIndex*/)
{
	try
	{
		commitIndex();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to the toggling of the "default object" button.
 */
void AttributeBuilderWidget::onDefaultToggled(bool checked)
{
	try
	{
		ui_->lineEditIndex->setEnabled(not checked);
		commitIndex();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a user modification on the type.
 */
void AttributeBuilderWidget::onTypeChanged(int newType)
{
	try
	{
		assert(newType >= OBJECT and newType <= STRING_SCALAR);
		AttributeTypeId type = (AttributeTypeId)newType;

		if(type != OBJECT)
		{
			ui_->valueEditor->setEditedType(TypeInfo::lookup(type));
		}

		enableDefault(type == OBJECT);
		enableIndex(type != OBJECT or !ui_->pushButtonDefault->isChecked());
		enableValue(type != OBJECT);

		commitType(type);
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a user modification on the value.
 */
void AttributeBuilderWidget::onValueChanged(const QString& /*newValue*/)
{
	try
	{
		commitValue();
	}
	CATCH_AND_REPORT
}

/**
 * Reacts to a click on the "Insert" button.
 */
void AttributeBuilderWidget::onInsertClicked()
{
	try
	{
		if(checkData())
		{
			bool insertAsChild = ui_->radioButtonChild->isChecked();
			emit insertClicked(model, ui_->lineEditPath->text(), insertAsChild);

			if(insertAsChild)
			{
				// switch to sibling insertion
				ui_->radioButtonSibling->setChecked(true);
			}
			incrementNumericIndex();
		}
	}
	CATCH_AND_REPORT
}


}
}
