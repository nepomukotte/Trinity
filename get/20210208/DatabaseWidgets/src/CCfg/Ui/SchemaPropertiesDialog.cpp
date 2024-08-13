#include "SchemaPropertiesDialog.h"
#include "CCfg/Ui/SchemaPropertiesDialog.ui.h"
#include <QMessageBox>

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::SchemaPropertiesDialog class methods body
////////////////////////////////////////////////////////////////////////////////

SchemaPropertiesDialog::SchemaPropertiesDialog(QWidget* parent) : QDialog(parent),
	widgets(new Ui_SchemaPropertiesDialog)
{
	widgets->setupUi(this);
}

SchemaPropertiesDialog::~SchemaPropertiesDialog()
{
}

std::string SchemaPropertiesDialog::getSchemaName() const
{
	return widgets->lineEditName->text().toStdString();
}

std::string SchemaPropertiesDialog::getPassword() const
{
	return widgets->lineEditPassword2->text().toStdString();
}

void SchemaPropertiesDialog::setSchemaName(const std::string& newName)
{
	widgets->lineEditName->setText(newName.c_str());
}

void SchemaPropertiesDialog::setNameEditable(bool isEditable)
{
	widgets->lineEditName->setEnabled(isEditable);
}

void SchemaPropertiesDialog::accept()
{
	if(widgets->lineEditPassword1->text() == widgets->lineEditPassword2->text())
	{
		QDialog::accept();
	}
	else
	{
		QMessageBox::information(this, "Input mismatch", "Passwords do not match.", QMessageBox::Close);
	}
}

}
}
