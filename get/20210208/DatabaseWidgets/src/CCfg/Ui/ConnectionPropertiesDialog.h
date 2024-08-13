#ifndef CCFG_UI_CONNECTIONPROPERTIESDIALOG_H_INCLUDED
#define CCFG_UI_CONNECTIONPROPERTIESDIALOG_H_INCLUDED

#include "CCfg/Io/DatabaseStorage.h"
#include <QDialog>
#include <memory>

// forward declarations
class Ui_ConnectionPropertiesDialog;

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This is a dialog box allowing to edit the properties of a connection.
 *
 * @author fchateau
 */
class ConnectionPropertiesDialog : public QDialog
{
	Q_OBJECT

	std::auto_ptr<Ui_ConnectionPropertiesDialog> widgets;

public:
	ConnectionPropertiesDialog(QWidget* parent = 0);
	virtual ~ConnectionPropertiesDialog();

	// @{ getters
	std::string getName() const;
	Io::DbmsType getDbmsType() const;
	std::string getUser() const;
	std::string getHost() const;
	std::string getPort() const;
	std::string getSID() const;
	std::string getSchema();
	// @}
	// @{ setters
	void setName(const std::string& newName);
	void setDbmsType(Io::DbmsType newDbms);
	void setUser(const std::string& newUser);
	void setHost(const std::string& newHost);
	void setPort(const std::string& newPort);
	void setSID(const std::string& newSID);
	void setSchema(const std::string& newSchema);
	// @}

protected:
	void setPortToDefault();

protected slots:
	void onDbmsChanged(const QString & text);
	void onUseDefaultPortToggled(bool checked);
	void onSelectSIDToggled();
};

}
}

#endif // CCFG_UI_CONNECTIONPROPERTIESDIALOG_H_INCLUDED
