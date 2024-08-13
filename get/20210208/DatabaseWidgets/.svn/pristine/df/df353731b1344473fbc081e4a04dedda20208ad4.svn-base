#include "ConnectionPropertiesDialog.h"
#include "CCfg/Ui/ConnectionPropertiesDialog.ui.h"
#include <QRegExpValidator>
#include <QFileDialog>
#include <QDir>
#include <boost/lexical_cast.hpp>
#include <string>

namespace CCfg
{
using namespace Io;
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::ConnectionPropertiesDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Constructor.
 */
ConnectionPropertiesDialog::ConnectionPropertiesDialog(QWidget* parent)
	: QDialog(parent), widgets(new Ui_ConnectionPropertiesDialog)
{
	widgets->setupUi(this);
	connect(widgets->comboBoxDbms, SIGNAL(currentIndexChanged(QString)),
			this, SLOT(onDbmsChanged(QString)));
	connect(widgets->checkBoxDefaultPort, SIGNAL(toggled(bool)),
			this, SLOT(onUseDefaultPortToggled(bool)));
	connect(widgets->buttonSelectSID, SIGNAL(released()),
			this, SLOT(onSelectSIDToggled()));
	widgets->lineEditName->setValidator(new QRegExpValidator(QRegExp("[-_a-zA-Z0-9]*"), this));
	widgets->spinBoxPort->setValue(DatabaseStorage::DEFAULT_PORT[getDbmsType()]);
#ifndef CCFG_USE_ORACLE
	widgets->comboBoxDbms->removeItem(widgets->comboBoxDbms->findText("Oracle"));
#endif
#ifndef CCFG_USE_POSTGRESQL
	widgets->comboBoxDbms->removeItem(widgets->comboBoxDbms->findText("PostgreSql"));
#endif
#ifndef CCFG_USE_SQLITE3
	widgets->comboBoxDbms->removeItem(widgets->comboBoxDbms->findText("SQLite3"));
#endif
}

/**
 * Destructor.
 */
ConnectionPropertiesDialog::~ConnectionPropertiesDialog()
{
}

std::string ConnectionPropertiesDialog::getName() const
{
	return widgets->lineEditName->text().toStdString();
}

Io::DbmsType ConnectionPropertiesDialog::getDbmsType() const
{
	const std::string dmbsString = widgets->comboBoxDbms->currentText().toLower().toStdString();
	return DatabaseStorage::dbmsFromString(dmbsString);
}

std::string ConnectionPropertiesDialog::getUser() const
{
	return widgets->lineEditUser->text().toStdString();
}

std::string ConnectionPropertiesDialog::getHost() const
{
	return widgets->lineEditHost->text().toStdString();
}

std::string ConnectionPropertiesDialog::getPort() const
{
	if (Io::Sqlite3 == getDbmsType())
		return QString("%1").arg(DatabaseStorage::DEFAULT_PORT[Sqlite3]).toStdString();
	return widgets->spinBoxPort->cleanText().toStdString();
}

std::string ConnectionPropertiesDialog::getSID() const
{
	return widgets->lineEditSID->text().toStdString();
}

std::string ConnectionPropertiesDialog::getSchema()
{
	return widgets->lineEditSchema->text().toStdString();
}

void ConnectionPropertiesDialog::setName(const std::string& newName)
{
	widgets->lineEditName->setText(newName.c_str());
}

void ConnectionPropertiesDialog::setDbmsType(Io::DbmsType newDbms)
{
	widgets->comboBoxDbms->setCurrentIndex(
			widgets->comboBoxDbms->findText(DatabaseStorage::dbmsToString(newDbms), static_cast<Qt::MatchFlags> (~Qt::MatchCaseSensitive)));
}

void ConnectionPropertiesDialog::setUser(const std::string& newUser)
{
	widgets->lineEditUser->setText(newUser.c_str());
}

void ConnectionPropertiesDialog::setHost(const std::string& newHost)
{
	widgets->lineEditHost->setText(newHost.c_str());
}

void ConnectionPropertiesDialog::setPort(const std::string& newPort)
{
	widgets->spinBoxPort->setValue(boost::lexical_cast<int>(newPort));
}

void ConnectionPropertiesDialog::setSID(const std::string& newSID)
{
	widgets->lineEditSID->setText(newSID.c_str());
}

void ConnectionPropertiesDialog::setSchema(const std::string& newSchema)
{
	widgets->lineEditSchema->setText(newSchema.c_str());
}

/**
 * Changes the value of the port to its default value for the selected database.
 */
void ConnectionPropertiesDialog::setPortToDefault()
{
	widgets->spinBoxPort->setValue(DatabaseStorage::DEFAULT_PORT[getDbmsType()]);
}

/**
 * Enables and disables some widgets and change the default port each time the
 * selected DBMS changes.
 */
void ConnectionPropertiesDialog::onDbmsChanged(const QString & text)
{
	if (text.isEmpty()) return;
	DbmsType type = DatabaseStorage::dbmsFromString(text.toLower().toStdString());
	widgets->lineEditHost->setEnabled(type != Sqlite3);
	if (type == Sqlite3)
		setHost(std::string());
	else
	{
		if (getHost().empty()) setHost("localhost");
	}
	widgets->lineEditUser->setEnabled(type != Sqlite3);
	widgets->lineEditSID->setEnabled(type == Oracle or type == PostgreSql or type == Sqlite3);
	widgets->lineEditSchema->setEnabled(type == MySql);
	if(widgets->checkBoxDefaultPort->isChecked())
	{
		setPortToDefault();
	}
	widgets->checkBoxDefaultPort->setEnabled(type != Sqlite3);
	widgets->buttonSelectSID->setEnabled(type == Sqlite3);
}

/**
 * Enables or disables the port editor when the "Default" check box state changes.
 */
void ConnectionPropertiesDialog::onUseDefaultPortToggled(bool checked)
{
	widgets->spinBoxPort->setEnabled(not checked);
	if(checked)
	{
		setPortToDefault();
	}
}

/**
 * Lets the user select an SQLite3 database file.
 */
void ConnectionPropertiesDialog::onSelectSIDToggled()
{
	QString fileName = QFileDialog::getOpenFileName(this, tr("Select SQLite3 database"), QDir::homePath());
	if (not fileName.isEmpty())
		widgets->lineEditSID->setText(fileName);
}

}
}
