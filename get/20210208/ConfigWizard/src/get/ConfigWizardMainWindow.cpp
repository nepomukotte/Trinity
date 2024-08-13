/**
 * @file ConfigWizardMainwindow.cpp
 * @date Nov 19, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: ConfigWizardMainWindow.cpp 1133 2014-02-11 09:34:39Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the GetConfigWizard software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "ConfigWizardMainWindow.h"
#include "get/ConfigWizardMainWindow.ui.h"
#include "get/cfg/ChannelConfigModel.h"
#include "get/geo/EditMappingDialog.h"
#include "get/geo/GeometryModel.h"
#include "get/root/HistogramDialog.h"
#include "CCfg/CConfig.h"
#include "CCfg/Control/Workspace.h"
#include "CCfg/Control/BaseConfigEditorController.h"
#include "CCfg/Io/Document.h"
#include "CCfg/Model/ConfigModelAccessor.h"
#include "CCfg/Ui/EditorMainWindow.h"
#include "CCfg/Ui/ConfigEditorTree.h"
#include "CCfg/Enum.h"
#include "Utils/Error/ExcepHandling.h"
#include "utl/Logging.h"
using namespace CCfg;

#include <TTree.h>

#include <QCloseEvent>
#include <QCheckBox>
#include <QComboBox>
#include <QDialog>
#include <QDialogButtonBox>
#include <QDockWidget>
#include <QFileDialog>
#include <QFileInfo>
#include <QLabel>
#include <QHBoxLayout>
#include <QInputDialog>
#include <QLabel>
#include <QLineEdit>
#include <QMessageBox>
#include <QMenuBar>
#include <QStatusBar>
#include <QString>
#include <QTableView>
#include <QToolBar>
#include <QVBoxLayout>
#include <QWizard>

#include <vector>
#include <stdexcept>

using CCfg::Model::ConfigNode;


namespace get
{
//_________________________________________________________________________________________________
ConfigWizardMainWindow::ConfigWizardMainWindow(QWidget* parentWidget)
	: Utils::Ui::DesignedWidget<Ui_ConfigWizardMainWindow, QMainWindow>(parentWidget), warnIfNotUndoable_(true)
{
	// Hide central widget
	// Note: Creating a main window without a central widget is not supported. You must have a central widget even if it is just a placeholder.
	QWidget* fakeCentralWidget = new QWidget(this);
	setCentralWidget(fakeCentralWidget);
	fakeCentralWidget->hide();

	// Create geometry mapping
	geoModel_ = new geo::GeometryModel(this);

	// Dock a config editor
	QDockWidget *cfgDockWidget = new QDockWidget(tr("Configuration Editor"), this);
	cfgDockWidget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	cfgEditor_ = new CCfg::Ui::EditorMainWindow(cfgDockWidget);
	cfgEditor_->setWindowFlags(Qt::Widget);
	//cfgEditor->menuBar()->hide();
	cfgDockWidget->setWidget(cfgEditor_);
	cfgDockWidget->setFeatures(cfgDockWidget->features() ^ QDockWidget::DockWidgetClosable);
	addDockWidget(Qt::LeftDockWidgetArea, cfgDockWidget);

	QDockWidget *channelsDockWidget = new QDockWidget(tr("Channel Editor"), this);
	channelsDockWidget->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);
	channelView_ = new QTableView(channelsDockWidget);
	channelView_->horizontalHeader()->setResizeMode(QHeaderView::ResizeToContents);
	channelView_->setSelectionBehavior(QAbstractItemView::SelectRows);
	channelModel_ = new cfg::ChannelConfigModel(channelsDockWidget);
	channelView_->setModel(channelModel_);
	channelsDockWidget->setWidget(channelView_);
	channelsDockWidget->setFeatures(channelsDockWidget->features() ^ QDockWidget::DockWidgetClosable);
	addDockWidget(Qt::RightDockWidgetArea, channelsDockWidget);

	// Create status bar
	/*QLabel* label = new QLabel(tr("Configuration: "));
	statusBar()->addPermanentWidget(label);
	currentConfigLabel_ = new QLabel(tr("-"));
	currentConfigLabel_->setMinimumWidth(10);
	statusBar()->addPermanentWidget(currentConfigLabel_);*/

	// Create actions, menu and toolbar
	createActions();
	createChannelActions();

	setMinimumSize(600, 400);
	setWindowTitle(tr("GET Configuration Wizard"));
}
//_________________________________________________________________________________________________
/**
 * Creates actions, menu and toolbar.
 */
void ConfigWizardMainWindow::createActions()
{
	cfgEditor_->getActionNew()->disconnect();
	connect(cfgEditor_->getActionNew(), SIGNAL(triggered()), SLOT(createNewConfig()));
	//cfgEditor_->getActionSaveAll()->disconnect();
	//cfgEditor_->getActionSaveAll()->setVisible(false);
	//cfgEditor_->getActionSaveAs()->disconnect();
	//cfgEditor_->getActionSaveAs()->setVisible(false);
	//cfgEditor_->getActionSave()->disconnect();
	//cfgEditor_->getActionSave()->setVisible(false);
	cfgEditor_->getActionAbout()->disconnect();
	connect(cfgEditor_->getActionAbout(), SIGNAL(triggered()), this, SLOT(about()));

	QMenu* menuGeometry = new QMenu(tr("Geometry"), this);
	cfgEditor_->getMenubar()->insertMenu(cfgEditor_->getMenuHelp()->menuAction(), menuGeometry);
	QAction* actionEditGeometry = menuGeometry->addAction(tr("&Edit geometry"));
	connect(actionEditGeometry, SIGNAL(triggered()), this, SLOT(editGeometry()));
	QMenu* menuDraw = menuGeometry->addMenu(tr("&Draw..."));
	QAction* actionDrawThreshold = menuDraw->addAction(tr("Draw &threshold"));
	connect(actionDrawThreshold, SIGNAL(triggered()), this, SLOT(drawThreshold()));
	QAction* actionDrawLsbThreshold = menuDraw->addAction(tr("Draw &LSB threshold"));
	connect(actionDrawLsbThreshold, SIGNAL(triggered()), this, SLOT(drawLsbThreshold()));

	// Connect signals to slots
	connect(cfgEditor_->getWorkspace(), SIGNAL(modelActivated(Model::ConfigModel*, Model::ConfigModel*)),
			this, SLOT(onConfigActivated(Model::ConfigModel*, Model::ConfigModel*)));

	// If you want all windows in a Mac application to share one menu bar, you must create a menu bar that does not have a parent.
#ifdef Q_OS_MAC
	QMenuBar* menuBar = cfgEditor_->menuBar();
	menuBar->setParent(0);
	setMenuBar(menuBar);
#endif
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::createChannelActions()
{
	QAction* actionGoToChannel = new QAction(tr("&Go To Channel"), channelView_);
	connect(actionGoToChannel, SIGNAL(triggered()), this, SLOT(goToCurrentChannel()));
	channelView_->addAction(actionGoToChannel);

	QMenu* menuEditThreshold = new QMenu(tr("Edit &DAC threshold..."), channelView_);
	QAction* actionEditLsbThreshold = menuEditThreshold->addAction(tr("Edit LSB threshold"));
	connect(actionEditLsbThreshold, SIGNAL(triggered()), this, SLOT(editLsbThreshold()));
	QAction* actionEditMsbThreshold = menuEditThreshold->addAction(tr("Edit MSB threshold"));
	connect(actionEditMsbThreshold, SIGNAL(triggered()), this, SLOT(editMsbThreshold()));
	channelView_->addAction(menuEditThreshold->menuAction());

	QAction* actionEditGain = new QAction(tr("Edit gain"), channelView_);
	connect(actionEditGain, SIGNAL(triggered()), this, SLOT(editGain()));
	channelView_->addAction(actionEditGain);
	QAction* actionSelectChannel = new QAction(tr("Edit test mode selection"), channelView_);
	connect(actionSelectChannel, SIGNAL(triggered()), this, SLOT(editSelect()));
	channelView_->addAction(actionSelectChannel);
	QAction* actionEditZeroSuppThreshold = new QAction(tr("Edit zero suppression threshold"), channelView_);
	connect(actionEditZeroSuppThreshold, SIGNAL(triggered()), this, SLOT(editZeroSuppressionThreshold()));
	channelView_->addAction(actionEditZeroSuppThreshold);
	QAction* actionEditReadoutCondition = new QAction(tr("Edit readout condition"), channelView_);
	connect(actionEditReadoutCondition, SIGNAL(triggered()), this, SLOT(editReadoutCondition()));
	channelView_->addAction(actionEditReadoutCondition);
	QAction* actionEditInhibition = new QAction(tr("Edit trigger inhibition"), channelView_);
	connect(actionEditInhibition, SIGNAL(triggered()), this, SLOT(editTriggerInhibition()));
	channelView_->addAction(actionEditInhibition);
	channelView_->setContextMenuPolicy(Qt::ActionsContextMenu);
}
//_________________________________________________________________________________________________
/**
 * Reacts to the event closing this main window.
 * It asks the user to save its documents. The user can interrupt the operation
 * if its unsaved documents requires additional processing before saving.
 */
void ConfigWizardMainWindow::closeEvent(QCloseEvent* event)
{
	try
	{
		if(cfgEditor_->getWorkspace()->closeAll())
		{
			event->accept();
			return;
		}
	}
	CATCH_AND_REPORT
	event->ignore();
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editGeometry()
{
	geo::EditMappingDialog dialog(this);
	dialog.setGeometryModel(geoModel_);
	dialog.exec();
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::about()
{
	QMessageBox::about(this, tr("About get-config-wizard"),
		tr("<h2> get-config-wizard </h2>"
		"<br><br>"
		"<b>This tool is an assistant for creating or editing configuration files for the GET run control software."
		"<br><br>"
		"&copy; 2014 <b>CEA, <a href=%1>IRFU</a></b><br><a href=%2>CeCILL</a> license"
		"<br><br>"
		"For further details, refer to the wiki of the GET collaboration at <a href=%3>%3</a>"
		" or send an email to <a href=mailto:%4>%4</a>."
		).arg("http://irfu.cea.fr").arg("http://www.cecill.info")
		.arg("https://dsm-trac.cea.fr/get").arg("get@cea.fr"));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::createNewConfig()
{
	// Get output file
	QStringList filters = QStringList()
			<< QString("%1 (*.xcfg").arg(tr("Configuration files"))
			<< QString("%1 (*.*)").arg("All files");
	QString	newConfigPath = QFileDialog::getSaveFileName(this, tr("Save new configuration as..."), QString("configure-pedestals.xcfg"), filters.join(";;"));
	if (newConfigPath.isEmpty()) return;

	// Create new config
	QString templatePath = QString("%1/../share/doc/GetController/examples/pedestals.xcfg").arg(QCoreApplication::applicationDirPath());
	if (not createNewConfig(templatePath.toStdString(), newConfigPath.toStdString())) return;

	// Close already opened config
	//cfgEditor_->getWorkspace()->closeAll();
	// Display new configuration
	cfgEditor_->openFiles(QStringList() << newConfigPath);
}
//_________________________________________________________________________________________________
/**
 * Creates a new configuration file,
 *  based on a template.
 */
bool ConfigWizardMainWindow::createNewConfig(const std::string & templatePath, const std::string & configPath)
{
	try
	{
		// Load template
		CCfg::Io::Document templateDoc;
		LOG_DEBUG() << "Loading template configuration from '" << templatePath << '\'';
		try
		{
			templateDoc.loadFromFile(templatePath);
		}
		catch (const CCfg::Io::Document::IoError & e)
		{
			LOG_ERROR() << e.what();
			return false;
		}
		Object* templateConfig = templateDoc.getConfig();

		// Create config from template
		std::auto_ptr< Object > newCfgObj = templateConfig->clone();

		// Save document
		CCfg::Io::Document doc;
		doc.setConfig(newCfgObj);
		doc.saveToFile(configPath, false);
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
		return false;
	}

	return true;
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::onConfigActivated(Model::ConfigModel* current, Model::ConfigModel*)
{
	channelModel_->setSourceModel(current);
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::goToCurrentChannel()
{
	QModelIndex channelIdx = channelModel_->mapToSource(channelView_->currentIndex());
	cfgEditor_->getWorkspace()->getCurrentController()->getView()->setCurrentIndex(channelIdx);
}
//_________________________________________________________________________________________________
bool ConfigWizardMainWindow::warnOpIsNotUndoable()
{
	if (not warnIfNotUndoable_) return true;

	QDialog dialog(this);
	QBoxLayout* layout = new QVBoxLayout(&dialog);
	layout->addWidget(new QLabel(tr("This operation is irreversible. Do you want to proceed?")));
	QCheckBox* alwaysWarnBox = new QCheckBox(tr("Always warn me"));
	alwaysWarnBox->setChecked(true);
	layout->addWidget(alwaysWarnBox);
	QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Yes | QDialogButtonBox::No);
	layout->addWidget(buttonBox);
	dialog.setLayout(layout);
	dialog.setWindowTitle(tr("Warning"));
	connect(buttonBox, SIGNAL(accepted()), &dialog, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), &dialog, SLOT(reject()));

	if (dialog.exec())
	{
		if (not alwaysWarnBox->isChecked()) warnIfNotUndoable_ = false;
		return true;
	}
	return false;
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editChannelParameter(const std::string parName, const CCfg::Value & value)
{
	// Store path of selected channels
	QModelIndexList rows = channelView_->selectionModel()->selectedRows();
	std::vector< Path > paths;
	for (int i=0; i < rows.size(); ++i)
	{
		paths.push_back(channelModel_->channelNode(rows[i])->getPath());
	}

	if (not warnOpIsNotUndoable()) return;

	// Save and close file
	cfgEditor_->getWorkspace()->save();
	std::string url = cfgEditor_->getWorkspace()->getCurrentModel()->getDocumentUrl();
	cfgEditor_->getWorkspace()->close();

	// Edit file
	std::auto_ptr< CCfg::Io::Document> doc(new CCfg::Io::Document);
	CCfg::CConfig cfg(&doc->loadFromFile(url));

	for (size_t i=0; i < paths.size(); ++i)
	{
		CCfg::CConfig param(cfg.appendPath(paths[i].subPath(1))(parName));
		param.setValue(value);
	}
	doc->saveToFile(url, false);

	// Reload file
	cfgEditor_->openFiles(QStringList() << QString::fromStdString(url));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editLsbThreshold()
{
	// Get new value
	bool ok;
	int value = QInputDialog::getInt(this, tr("Edit LSB DAC threshold"), tr("Value:"), 15, 0, 15, 1, &ok);
	if (not ok) return;

	// Set parameter
	editChannelParameter("LSBThresholdValue", CCfg::Value(value));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editGain()
{
	// Get allowed values
	QModelIndexList rows = channelView_->selectionModel()->selectedRows();
	ConfigNode* channelNode = channelModel_->channelNode(rows.first());
	const Enum* const enumeration = dynamic_cast<const Enum*>(channelNode->search("Gain", "")->getRange());
	const Enum::Set& values = enumeration->getValues();
	QStringList items;
	for (Enum::Set::const_iterator iVal = values.begin(); iVal != values.end(); ++iVal)
	{
		items << QString::fromStdString(iVal->getAsString());
	}

	// Get new value
	bool ok;
	QString value = QInputDialog::getItem(this, tr("Edit gain"), tr("Value:"), items, 0, false, &ok);
	if (not ok) return;

	// Set parameter
	editChannelParameter("Gain", CCfg::Value(value.toStdString()));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editSelect()
{
	// Get new value
	bool ok;
	QStringList items = QStringList() << "true" << "false";
	QString stringValue = QInputDialog::getItem(this, tr("Select channel(s) for test mode"), tr("Value:"), items, 0, false, &ok);
	if (not ok) return;
	bool value = (QString("true") == stringValue);

	// Set parameter
	editChannelParameter("isSelectedforTestMode", CCfg::Value(value));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editZeroSuppressionThreshold()
{
	// Get new value
	bool ok;
	int value = QInputDialog::getInt(this, tr("Edit zero suppression threshold"), tr("Threshold:"), 600, 0, 4095, 5, &ok);
	if (not ok) return;

	// Set parameter
	editChannelParameter("zeroSuppressionThreshold", CCfg::Value(value));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editReadoutCondition()
{
	// Get allowed values
	QModelIndexList rows = channelView_->selectionModel()->selectedRows();
	ConfigNode* channelNode = channelModel_->channelNode(rows.first());
	const Enum* const enumeration = dynamic_cast<const Enum*>(channelNode->search("Reading", "")->getRange());
	const Enum::Set& values = enumeration->getValues();
	QStringList items;
	for (Enum::Set::const_iterator iVal = values.begin(); iVal != values.end(); ++iVal)
	{
		items << QString::fromStdString(iVal->getAsString());
	}

	// Get new value
	bool ok;
	QString value = QInputDialog::getItem(this, tr("Edit readout condition"), tr("Value:"), items, 0, false, &ok);
	if (not ok) return;

	// Set parameter
	editChannelParameter("Reading", CCfg::Value(value.toStdString()));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editTriggerInhibition()
{
	// Get allowed values
	QModelIndexList rows = channelView_->selectionModel()->selectedRows();
	ConfigNode* channelNode = channelModel_->channelNode(rows.first());
	const Enum* const enumeration = dynamic_cast<const Enum*>(channelNode->search("TriggerInhibition", "")->getRange());
	const Enum::Set& values = enumeration->getValues();
	QStringList items;
	for (Enum::Set::const_iterator iVal = values.begin(); iVal != values.end(); ++iVal)
	{
		items << QString::fromStdString(iVal->getAsString());
	}

	// Get new value
	bool ok;
	QString value = QInputDialog::getItem(this, tr("Edit trigger inhibition"), tr("Value:"), items, 0, false, &ok);
	if (not ok) return;

	// Set parameter
	editChannelParameter("TriggerInhibition", CCfg::Value(value.toStdString()));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editChipParameter(const std::string parameterPath, const CCfg::Value & value)
{
	// Store path of selected chips
	QModelIndexList rows = channelView_->selectionModel()->selectedRows();
	std::vector< Path > paths;
	for (int i=0; i < rows.size(); ++i)
	{
		paths.push_back(channelModel_->agetNode(rows[i])->getPath());
	}

	if (not warnOpIsNotUndoable())
		return;

	if ( QMessageBox::Yes != QMessageBox::warning(this, tr("Warning"),
			tr("You are editing an AGET level parameter.\n"
			"It might affect more than just the selected channels.\nDo you want to proceed?"),
			QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Yes))
		return;

	// Save and close file
	cfgEditor_->getWorkspace()->save();
	std::string url = cfgEditor_->getWorkspace()->getCurrentModel()->getDocumentUrl();
	cfgEditor_->getWorkspace()->close();

	// Edit file
	std::auto_ptr< CCfg::Io::Document> doc(new CCfg::Io::Document);
	CCfg::CConfig cfg(&doc->loadFromFile(url));

	for (size_t i=0; i < paths.size(); ++i)
	{
		CCfg::Path fullPath = paths[i].append(CCfg::Path(parameterPath));
		CCfg::CConfig param(cfg.appendPath(fullPath.subPath(1)));
		//LOG_DEBUG() << param.getPath();
		param.setValue(value);
	}
	doc->saveToFile(url, false);

	// Reload file
	cfgEditor_->openFiles(QStringList() << QString::fromStdString(url));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::editMsbThreshold()
{
	// Get new value
	bool ok;
	int value = QInputDialog::getInt(this, tr("Edit MSB DAC threshold"), tr("Value:"), 7, 0, 7, 1, &ok);
	if (not ok) return;

	// Set parameter
	editChipParameter("Global.Reg1.GlobalThresholdValue", CCfg::Value(value));
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::drawThreshold()
{
#ifdef HAVE_QTROOT
	TH2F* hist = geoModel_->createThresholdHistogram(channelModel_);
	root::HistogramDialog dialog(this);
	dialog.addHistogram(hist, "Threshold", "cThreshold");
	dialog.exec();
#endif
}
//_________________________________________________________________________________________________
void ConfigWizardMainWindow::drawLsbThreshold()
{
#ifdef HAVE_QTROOT
	TH2F* hist = geoModel_->createLsbThresholdHistogram(channelModel_);
	root::HistogramDialog dialog(this);
	dialog.addHistogram(hist, "LSB threshold", "cLsbThreshold");
	dialog.exec();
#endif
}
//_________________________________________________________________________________________________
} /* namespace get */
