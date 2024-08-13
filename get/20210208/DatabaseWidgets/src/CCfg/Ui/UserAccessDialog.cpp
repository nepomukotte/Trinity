#include "UserAccessDialog.h"
#include "CCfg/Ui/UserAccessDialog.ui.h"
#include "Utils/Error/ExcepHandling.h"
#include <QPushButton>

namespace CCfg
{
namespace Ui
{
using namespace Io::Db;

////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::UserAccessDialog class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Full constructor.
 */
UserAccessDialog::UserAccessDialog(Io::Db::ConfigDao& dao, QWidget* parent)
	: QDialog(parent), widgets(new Ui_UserAccessDialog), dao(dao)
{
	widgets->setupUi(this);
	initAccessList();
	loadUserList();
	loadSchemaList();
	connect(widgets->comboBoxSchema, SIGNAL(activated(int)),
			this, SLOT(onSchemaActivated(int)));
	connect(widgets->comboBoxUser, SIGNAL(activated(int)),
			this, SLOT(onUserActivated(int)));
	connect(widgets->comboBoxAccess, SIGNAL(activated(int)),
			this, SLOT(onAccessActivated(int)));
	connect(widgets->buttonBox->button(QDialogButtonBox::Apply), SIGNAL(clicked()),
			this, SLOT(apply()));
	connect(widgets->buttonBox->button(QDialogButtonBox::Discard), SIGNAL(clicked()),
			this, SLOT(discard()));
	onSchemaActivated(0);
	setDirty(false);
}

/**
 * Destructor.
 */
UserAccessDialog::~UserAccessDialog()
{
}

/**
 * Set the current schema in the combo box.
 */
void UserAccessDialog::setSchema(const std::string& schema)
{
	widgets->comboBoxSchema->setCurrentIndex(widgets->comboBoxSchema->findText(schema.c_str()));
}

/**
 * Returns whether the access of the current user on the current schema has been changed.
 */
bool UserAccessDialog::isDirty() const
{
	return widgets->buttonBox->button(QDialogButtonBox::Apply)->isEnabled();
}

/**
 * Applies the selected access to the user on the schema.
 */
void UserAccessDialog::apply()
{
	try
	{
		dao.setRole(
			widgets->comboBoxUser->currentText().toStdString(),
			static_cast<ConfigUserRole>(widgets->comboBoxAccess->currentIndex()));
		setDirty(false);
	}
	CATCH_AND_REPORT
}

/**
 * Discards the access modifications and restores its current value.
 */
void UserAccessDialog::discard()
{
	try
	{
		setDirty(false);
		onUserActivated(widgets->comboBoxUser->currentIndex()); // reload the access.
	}
	CATCH_AND_REPORT
}

/**
 * Loads the list of users from the database to the combo box.
 */
void UserAccessDialog::loadUserList()
{
	std::vector<std::string> users = dao.getUserNames();
	for(size_t i = 0; i < users.size(); ++i)
	{
		widgets->comboBoxUser->addItem(users[i].c_str());
	}
}

/**
 * Loads the list of schema from the database to the combo box.
 * Only schema containing all tables of the configuration system are listed.
 */
void UserAccessDialog::loadSchemaList()
{
	ConfigDatabasePtrSetPtr schemaList = dao.getDatabases();
	ConfigDatabaseNameIdx& nameIdx = schemaList->get<Name>();
	for(ConfigDatabaseNameIdx::iterator curDb = nameIdx.begin(), endDb = nameIdx.end();
		curDb != endDb; ++curDb)
	{
		const std::string& schemaName = (*curDb)->name;
		dao.setWorkingDatabase(schemaName);
		if(dao.checkHasAllTables())
		{
			widgets->comboBoxSchema->addItem(schemaName.c_str());
		}
	}
}

/**
 * Fills the access combo box with its items.
 */
void UserAccessDialog::initAccessList()
{
	widgets->comboBoxAccess->addItem("Expert");
	widgets->comboBoxAccess->addItem("Shifter");
	widgets->comboBoxAccess->addItem("Analyst");
	widgets->comboBoxAccess->addItem("No access");
}

/**
 * Updates the access and enables/disables widgets when a different schema is selected.
 */
void UserAccessDialog::onSchemaActivated(int index)
{
	try
	{
		widgets->comboBoxUser->setEnabled(index != -1);
		widgets->comboBoxAccess->setEnabled(index != -1);
		if(index != -1)
		{
			dao.setWorkingDatabase(widgets->comboBoxSchema->itemText(index).toStdString());
			onUserActivated(widgets->comboBoxUser->currentIndex());
		}
	}
	CATCH_AND_REPORT
}

/**
 * Updates the access when a different user is selected.
 */
void UserAccessDialog::onUserActivated(int index)
{
	try
	{
		if(index != -1)
		{
			std::string user = widgets->comboBoxUser->itemText(index).toStdString();
			bool isAdmin = dao.isSchemaAdmin(user);

			ConfigUserRole access = isAdmin ? Expert : dao.getRole(user);
			widgets->comboBoxAccess->setCurrentIndex(access);

			widgets->checkBoxAdmin->setChecked(isAdmin);
			widgets->comboBoxAccess->setEnabled(not isAdmin);
		}
	}
	CATCH_AND_REPORT
}

/**
 * Updates buttons when a different access is selected.
 */
void UserAccessDialog::onAccessActivated(int /*index*/)
{
	try
	{
		setDirty(true);
	}
	CATCH_AND_REPORT
}

/**
 * Set the access as modified or not.
 */
void UserAccessDialog::setDirty(bool dirty)
{
	try
	{
		widgets->buttonBox->button(QDialogButtonBox::Apply)->setEnabled(dirty);
		widgets->buttonBox->button(QDialogButtonBox::Discard)->setEnabled(dirty);
		widgets->comboBoxUser->setEnabled(not dirty);
		widgets->comboBoxSchema->setEnabled(not dirty);
	}
	CATCH_AND_REPORT
}

}
}
