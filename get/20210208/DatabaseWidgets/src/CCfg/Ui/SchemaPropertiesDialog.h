#ifndef CCFG_UI_SCHEMAPROPERTIESDIALOG_H_INCLUDED
#define CCFG_UI_SCHEMAPROPERTIESDIALOG_H_INCLUDED

#include <QDialog>
#include <memory>

// forward declarations
class Ui_SchemaPropertiesDialog;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 *
 * @author fchateau
 */
class SchemaPropertiesDialog : public QDialog
{
	std::auto_ptr<Ui_SchemaPropertiesDialog> widgets;

public:
	SchemaPropertiesDialog(QWidget* parent = 0);
	~SchemaPropertiesDialog();

	std::string getSchemaName() const;
	std::string getPassword() const;

	void setSchemaName(const std::string& newName);

	void setNameEditable(bool isEditable);

public slots:
	void accept();
};

}
}

#endif // CCFG_UI_SCHEMAPROPERTIESDIALOG_H_INCLUDED
