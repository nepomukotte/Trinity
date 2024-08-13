/* =====================================================================================================================
 * GetTest.cpp
 * ---------------------------------------------------------------------------------------------------------------------
 * class GetTest
 * Created on: 24 August 2010 at Irfu/Sedi/Ldef, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Frederic Druillole (frederic.druillole@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 *               Patrick Sizun (patrick.sizun@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GET (General Electronics for TPC) data acquisition framework @ CEA-Irfu.
 * Find out more on: http://www-actar-get.cea.fr
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#include "get/GetTest.h"
#include "get/utl/SlotsLoggingBackend.h"
#include "get/ChannelMask.h"
#include "get/Console.h"
#include "get/viewgraph.h"
#include "get/viewgraph.h"
#include "get/DataBloc.h"
#include "get/EccClient.h"
#include "get/CatchAllExceptions.h"
#include "get/CoBoEvent.h"
#include "get/ControllerOptions.h"
#include "get/EditTestDialog.h"
#include "get/GetControllerException.h"
#include "get/SelectTestDialog.h"
#include "get/TestController.h"
#include "get/TestSummaryWidget.h"
#include "get/WorkspaceManager.h"
#include "get/WorkspaceModel.h"
#include "get/WorkspaceWizard.h"
#include "get/EditChannelMaskDialog.h"
#include "get/daq/DataRouter.h"
#include "get/daq/DataRouterWizard.h"
#include "get/daq/DiscardStorage.h"
#include "get/daq/CustomStorage.h"
#include "get/QwtSymbolDialog.h"
#include "get/mt/AlarmLogger.h"
#include "mdaq/utl/ConsoleLoggingBackend.h"
#include "mdaq/utl/Server.h"
#include "mfm/Frame.h"
#include "CCfg/Ui/EditConfigDialog.h"
#include "utl/net/SocketAddress.h"
#include "utl/Logging.h"

using namespace get;
using namespace ::utl::net;
using namespace get::utl;
using mdaq::utl::Server;

#include <QDockWidget>
#include <QFileInfo>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QMessageBox>
#include <QSettings>
#include <QSignalMapper>
#include <QStatusBar>
#include <QStringList>
#include <QTableView>
#include <QTimer>
#include <QCoreApplication>

#include <qwt_plot.h>
#include <boost/algorithm/string/case_conv.hpp>
namespace ba = boost::algorithm;
#include <vector>
#include <stdexcept>

#undef MSG_IF
#define MSG_IF(cond)  if(!(cond)) ; else LOG_DEBUG()

/** Constructs GetController main window.
 */
GetTest::GetTest(QWidget* parentWidget) :
	QMainWindow(parentWidget),
	subWindowMenu(0),
	testController(0),
	XOscillo(10000), YOscillo(10000),
	isHeaderDisplayEnabled_(true)
{
	setWindowTitle(tr("GET Controller"));
	setWindowIcon(QPixmap(":/get/logo-highres.png"));

	QObject::connect(qApp, SIGNAL(aboutToQuit()), this, SLOT(cleanUp()));

	// Create central widget
	ui.setupUi(this);
	centralZone = new QMdiArea(this);
	centralZone->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	centralZone->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
	centralZone->setViewMode(QMdiArea::SubWindowView);
	centralZone->setMinimumSize(600, 480);
	setCentralWidget(centralZone);
	QObject::connect(centralZone, SIGNAL(subWindowActivated(QMdiSubWindow *)), this, SLOT(updateSubWindowMenu()));
	windowMapper = new QSignalMapper(this);
	QObject::connect(windowMapper, SIGNAL(mapped(QWidget *)), this, SLOT(setActiveSubWindow(QWidget *)));

	// Create command panel
	QWidget *commandPanel = new QWidget(this);
	uicW.setupUi(commandPanel);
	setControlStatus(Ctrl::Off);
	ui.actionNewTest->setEnabled(true);
	ui.actionOpenTest->setEnabled(true);
	controlDock = new QDockWidget(tr("Control"), this);
	addDockWidget(Qt::RightDockWidgetArea, controlDock);
	controlDock->setWidget(commandPanel);

	// Create test info panel
	TestSummaryWidget* summaryPanel = new TestSummaryWidget(this);
	QObject::connect(this, SIGNAL(testSelected()), summaryPanel, SLOT(updateData()));
	QObject::connect(this, SIGNAL(testModified()), summaryPanel, SLOT(updateData()));
	QObject::connect(this, SIGNAL(dataRouterSettingsChanged()), summaryPanel, SLOT(setDataRouterInfo()));
	QObject::connect(this, SIGNAL(workspaceLoaded()), summaryPanel, SLOT(clear()));
	testInfoDock = new QDockWidget(tr("Information"), this);
	testInfoDock->setWidget(summaryPanel);
	addDockWidget(Qt::LeftDockWidgetArea, testInfoDock);

	// Instantiate log console
	msgConsole = new Console(this);
	consWin = new QMdiSubWindow(this, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	consWin->setWidget(msgConsole);
	consWin->setAttribute(Qt::WA_Disabled);
	consWin->setWindowTitle(tr("Log"));
	consWin = centralZone->addSubWindow(consWin, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);

	// Instantiate error log console
	/*errorConsole = new Console(this);
	errorLog = new QMdiSubWindow(this, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	errorLog->setWidget(errorConsole);
	errorLog->setAttribute(Qt::WA_Disabled);
	errorLog->setWindowTitle(tr("Errors"));
	errorLog = centralZone->addSubWindow(errorLog, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	*/

	// Instantiate header view console
	Console *headerView = new Console(this);
	HeaderLog = new QMdiSubWindow(this, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	HeaderLog->setWidget(headerView);
	HeaderLog->setAttribute(Qt::WA_Disabled);
	HeaderLog->setWindowTitle(tr("Headers"));
	HeaderLog = centralZone->addSubWindow(HeaderLog, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	connect(this, SIGNAL(logHd(const QString &)), headerView, SLOT(displayMessage(const QString &)));

	// Initialize channels to display
	channels = ChannelMask::fromString("0-67");

	// Create and start data acquisition thread
	try
	{
		dataBloc.reset(new DataBloc());
	}
	CATCH_ALL_EXCEPTIONS()

	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
		busyDisplaying[asadIdx] = false;

	// Check preferred curve symbol
	curveSymbolStyle = QSettings().value("symbolStyle", "Cross").toString();

	// Register user types to be used as signal parameters accross threads
	qRegisterMetaType< get::CoBoEvent >("get::CoBoEvent");

	// Connect actions for menus, toolbar and widgets
	createActions();
	connectActions();
	updateSubWindowMenu();
	setTabbedViewMode(true);
	centralZone->setActiveSubWindow(consWin);

	// Use GUI as logging backend
	SlotsLoggingBackend* loggingBackend = new SlotsLoggingBackend();
	QObject::connect(loggingBackend, SIGNAL(log(const QString &)),
			msgConsole, SLOT(displayMessage(const QString &)));
	QObject::connect(loggingBackend, SIGNAL(error(const QString &)),
			msgConsole, SLOT(displayMessage(const QString &)));
	QObject::connect(loggingBackend, SIGNAL(error(const QString &)),
			this, SLOT(showErrors()));
	::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(loggingBackend));

	connect(this, SIGNAL(eventReady(const get::CoBoEvent &)), this, SLOT(displayEvent(const get::CoBoEvent &)));
	for (size_t asadIdx = 0; asadIdx < 4; ++asadIdx)
		displayTime[asadIdx].start();
}
//_________________________________________________________________________________________________
void GetTest::setEditTestEnabled(bool enabled)
{
	ui.actionNewTest->setEnabled(enabled);
	ui.actionOpenTest->setEnabled(enabled);
	ui.actionEditTest->setEnabled(enabled);
	ui.actionEditHardware->setEnabled(enabled);
	ui.actionEditConditions->setEnabled(enabled);
	emit editTestEnabled(enabled);
}
//_________________________________________________________________________________________________
void GetTest::cleanUp()
{
    // Set logging backend
	::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ::mdaq::utl::ConsoleLoggingBackend));

	if (testController.get())
	{
		stopDaq();
		testController->wait(5000); // Wait 5s for thread to terminate
	}
	destroyDataRouter();
}
//_________________________________________________________________________________________________
void GetTest::showMessages()
{
	centralZone->setActiveSubWindow(consWin);
}
//_________________________________________________________________________________________________
void GetTest::showErrors()
{
	//centralZone->setActiveSubWindow(errorLog);
	centralZone->setActiveSubWindow(consWin);
}
//_________________________________________________________________________________________________
bool GetTest::selectDataRouter()
{
	daq::DataRouterWizard dialog(this);
	if (QDialog::Accepted == dialog.exec())
	{
		return createDataRouter();
	}
	return false;
}
//__________________________________________________________________________________________________
bool GetTest::createDataRouter()
{
	try
	{
		emit dataRouterSettingsChanged();
		if (Ctrl::Off != controlStatus())
			setControlStatus(Ctrl::Idle);
		// In all cases, we first destroy embedded router
		destroyDataRouter();

		const ControllerOptions & options = ControllerOptions::instance();

		// If needed, we create an internal data router
		if (not options.isDataRouterExternal())
		{
			// Get router parameters
			SocketAddress flowEndpoint = options.dataFlowEndpoint();
			SocketAddress ctrlEndpoint = options.dataCtrlEndpoint();
			std::string senderFlowType = options.dataFlowType();

			std::string routerFlowType = "TCP";
			if (senderFlowType == "ICE")
				routerFlowType = "ICE";
			else if (senderFlowType == "FDT")
				routerFlowType = "FDT";
			else if (senderFlowType == "Debug")
				return true;

			LOG_INFO() << "Starting data router";
			LOG_INFO() << "Control endpoint: " << ctrlEndpoint;
			LOG_INFO() << "Flow endpoint: " << flowEndpoint << "; flow type: " << senderFlowType;

			// Create data processor
			std::auto_ptr < mdaq::daq::DataProcessorCore > dataProcessor;
			if (dataBloc.get() != 0)
			{
				dataProcessor.reset(new daq::CustomStorage(*dataBloc.get()));
				const daq::FramePublisher* publisher =
						reinterpret_cast< daq::CustomStorage* >(dataProcessor.get())->publisher();
				QObject::connect(publisher, SIGNAL(frameReceived()),
					this, SLOT(onFrameReceived()));
			}
			else
				dataProcessor.reset(new daq::DiscardStorage());

			// Create router
			mdaq::utl::Server & ctrlServer = mdaq::utl::Server::create(ctrlEndpoint);
			::IceUtil::Handle< daq::DataRouter > dataRouter(
					new daq::DataRouter(flowEndpoint, routerFlowType, dataProcessor));
			ctrlServer.addServant("DataRouter", dataRouter);
			ctrlServer.start();
			startRunProcessor();
		}
		return true;
	}
	CATCH_ALL_EXCEPTIONS()
	destroyDataRouter();
	return false;
}
//__________________________________________________________________________________________________
void GetTest::destroyDataRouter()
{
	LOG_DEBUG() << "Destroying data router server...";
	mdaq::utl::Server::destroy(ControllerOptions::instance().dataCtrlEndpoint());
	// Reset frame counter
	uicW.frameCountBox->setValue(0);
	daqrun = false;
}
//__________________________________________________________________________________________________
bool GetTest::createAlarmLogger(const ::utl::net::SocketAddress & address)
{
	LOG_DEBUG() << "Creating alarm logger with address " << address;
	try
	{
		try
		{
			Server::find(address).addServant("AlarmLogger", new mt::AlarmLogger);
		}
		catch (const Server::NotCreatedException & e)
		{
			Server::create(address).addServant("AlarmLogger", new mt::AlarmLogger).start();
		}
		ControllerOptions::instance().setAlarmLoggerEndpoint(address);
	}
	catch (const Ice::LocalException & e)
	{
		LOG_ERROR() << "Could not create alarm logger: " << e.ice_name();
		return false;
	}
	return true;
}
//__________________________________________________________________________________________________
mdaq::daq::RunProcessorPrx GetTest::runProcessor()
{
	std::string proxyString = ControllerOptions::instance().dataCtrlEndpoint().iceProxyString("DataRouter");
	Ice::CommunicatorPtr commmunicator = ::mdaq::utl::Server::ic();
	Ice::ObjectPrx proxy = commmunicator->stringToProxy(proxyString)->ice_timeout(5000); // Timeout in milliseconds
	return mdaq::daq::RunProcessorPrx::checkedCast(proxy);
}
//__________________________________________________________________________________________________
void GetTest::startRunProcessor()
{
	try
	{
		runProcessor()->runStart();
	}
	catch (const ::Ice::LocalException & e)
	{
		LOG_WARN() << "Error starting run processor: " << e.ice_name();
	}
}
//__________________________________________________________________________________________________
void GetTest::stopRunProcessor()
{
	LOG_INFO() << "Stopping run processor";
	try
	{
		runProcessor()->runStop();
	}
	catch (const ::Ice::LocalException & e)
	{
		LOG_WARN() << "Error stopping run processor: " << e.ice_name();
	}
}
//__________________________________________________________________________________________________
void GetTest::resetRunProcessor()
{
	startRunProcessor();
}
//__________________________________________________________________________________________________
void GetTest::frameAvailable(mfm::Frame& getFrame)
{
	const size_t asadIdx = getFrame.headerField(27u, 1u).value<uint8_t>();
	if (asadIdx >= 4) return;
	if (displayTime[asadIdx].secsTo(QTime::currentTime()) < 0) // If current time has crossed midnight
	{
		displayTime[asadIdx].start();
	}
	else if (not busyDisplaying[asadIdx] and displayTime[asadIdx].secsTo(QTime::currentTime()) > 1)
	{
		displayTime[asadIdx].start();
		displayHeader(getFrame);
		//coboEvent[asadIdx].fromFrame_hardcoded(getFrame);
		coboEvent[asadIdx].fromFrame(getFrame);
		busyDisplaying[asadIdx] = true;
		emit eventReady(coboEvent[asadIdx]);
	}
}
//__________________________________________________________________________________________________
void GetTest::createActions()
{
#if QT_VERSION >= 0x040600
	ui.actionNewTest->setIcon(QIcon::fromTheme("document-new"));
	ui.actionOpenTest->setIcon(QIcon::fromTheme("document-open"));
	ui.actionEditTest->setIcon(QIcon::fromTheme("document-properties"));
	ui.actionQuit->setIcon(QIcon::fromTheme("application-exit"));
	ui.actionAbout->setIcon(QIcon::fromTheme("help-about"));
#endif
#if QT_VERSION >= 0x040200
	ui.actionAbout->setMenuRole(QAction::AboutRole);
#endif
	ui.actionNewTest->setShortcuts(QKeySequence::New);
	ui.actionOpenTest->setShortcuts(QKeySequence::Open);
	ui.actionQuit->setShortcuts(QKeySequence::Quit);
	ui.actionAbout->setShortcuts(QKeySequence::HelpContents);
	ui.actionEditHardware->setIcon(QIcon(":/application/cconfiged/logo.svg"));
	ui.actionEditConditions->setIcon(QIcon(":/application/cconfiged/logo.svg"));
	ui.actionSaveRawFrames->setChecked(ControllerOptions::instance().isSaveRawFrameEnabled());
	ui.actionDisplayHeader->setChecked(ControllerOptions::instance().isHeaderDisplayEnabled());
}

void GetTest::connectActions()
{
	connect(ui.actionNewTest, SIGNAL(triggered()), this, SLOT(createTest()));
	connect(ui.actionOpenTest, SIGNAL(triggered()), this, SLOT(selectTest()));
	connect(ui.actionEditTest, SIGNAL(triggered()), this, SLOT(editTest()));
	connect(ui.actionEditConditions, SIGNAL(triggered()), this, SLOT(editConditionsCConfig()));
	connect(ui.actionEditHardware, SIGNAL(triggered()), this, SLOT(editHardwareCConfig()));
	connect(ui.actionSwitchWorkspace, SIGNAL(triggered()), this, SLOT(switchWorkspace()));
	connect(ui.actionQuit, SIGNAL(triggered()), qApp, SLOT(quit()));
	connect(ui.actionSetTabbedViewMode, SIGNAL(triggered(bool)), this, SLOT(setTabbedViewMode(bool)));
	connect(ui.actionTile_Windows, SIGNAL(triggered()), centralZone, SLOT(tileSubWindows()));
	connect(ui.actionCascaded_Windows, SIGNAL(triggered()), centralZone, SLOT(cascadeSubWindows()));
	connect(ui.actionSaveRawFrames, SIGNAL(triggered(bool)), this, SLOT(setSaveRawFramesEnabled(bool)));
	connect(ui.actionSetChannelDisplay, SIGNAL(triggered()), this, SLOT(selectChannels()));
	connect(ui.actionSelectSymbol, SIGNAL(triggered()), this, SLOT(selectCurveSymbol()));
	connect(ui.actionDisplayHeader, SIGNAL(triggered(bool)), this, SLOT(setHeaderDisplayEnabled(bool)));
	connect(ui.actionSelectDataRouter, SIGNAL(triggered()), this, SLOT(selectDataRouter()));
	connect(ui.actionAbout, SIGNAL(triggered()), this, SLOT(about()));

	TestSummaryWidget* summaryWidget = qobject_cast<TestSummaryWidget*>(testInfoDock->widget());
	summaryWidget->setConfigActions(ui.actionEditHardware, ui.actionEditConditions);
	summaryWidget->setWorkspaceActions(ui.actionNewTest, ui.actionOpenTest, ui.actionEditTest);

	// State machine control widget
	connect(uicW.pbLoadHw, SIGNAL(clicked()), this, SLOT(showMessages()));
	connect(uicW.pbLoadHw, SIGNAL(clicked()), this, SLOT(loadHw()));
	connect(uicW.pbConfig, SIGNAL(clicked()), this, SLOT(configure()));
	connect(uicW.pbLaunch, SIGNAL(clicked()), this, SLOT(startDaq()));
	connect(uicW.pbStop, SIGNAL(clicked()), this, SLOT(stopDaq()));
	connect(uicW.pbReset, SIGNAL(clicked()), this, SLOT(Reset()));
}
//__________________________________________________________________________________________________
void GetTest::updateSubWindowMenu()
{
	// Create sub-menu
	if (not subWindowMenu)
	{
		subWindowMenu = new QMenu(tr("Windows"), this);
		ui.actionShowWindows->setMenu(subWindowMenu);
	}

	// Fill menu
	subWindowMenu->clear();

	QList<QMdiSubWindow *> subWindows = centralZone->subWindowList();
	for (int i = 0; i < subWindows.size(); ++i) {
		QMdiSubWindow *subWindow = subWindows.at(i);
		QAction *action  = subWindowMenu->addAction(subWindow->windowTitle());
		action->setCheckable(true);
		action->setChecked(subWindow == centralZone->activeSubWindow());
		QObject::connect(action, SIGNAL(triggered()), windowMapper, SLOT(map()));
		windowMapper->setMapping(action, subWindows.at(i));
	}
	subWindowMenu->addSeparator();
	subWindowMenu->addAction(testInfoDock->toggleViewAction());
	subWindowMenu->addAction(controlDock->toggleViewAction());
}
//__________________________________________________________________________________________________
void GetTest::setActiveSubWindow(QWidget *window)
{
	if (!window)
		return;
	centralZone->setActiveSubWindow(qobject_cast<QMdiSubWindow *>(window));
}
//__________________________________________________________________________________________________
/**
 * Slot to edit compound configuration describing conditions of current test.
 */
void GetTest::editConditionsCConfig()
{
	try
	{
		std::string hardwarePath = WorkspaceManager::instance().currentTest()->paramConfigPath();
		CCfg::Ui::EditConfigDialog::editConfig(this, hardwarePath,
				tr("Edit test conditions"));
	}
	catch (const CCfg::Exception & e)
	{
		QMessageBox::warning(this, tr("Error"), QString::fromStdString(e.what()));
	}
	catch (const GetControllerException & e)
	{
		QMessageBox::warning(this, tr("Error"), QString::fromStdString(e.what()));
	}
}
//__________________________________________________________________________________________________
/**
 * Slot to edit compound configuration describing hardware of current test.
 */
void GetTest::editHardwareCConfig()
{
	try
	{
		std::string hardwarePath = WorkspaceManager::instance().currentTest()->hwConfigPath();
		CCfg::Ui::EditConfigDialog::editConfig(this, hardwarePath,
				tr("Edit hardware configuration"));
	}
	catch (const CCfg::Exception & e)
	{
		QMessageBox::warning(this, tr("Error"), QString::fromStdString(e.what()));
	}
	catch (const GetControllerException & e)
	{
		QMessageBox::warning(this, tr("Error"), QString::fromStdString(e.what()));
	}
}
//__________________________________________________________________________________________________
void GetTest::switchWorkspace()
{
	setEditTestEnabled(false);

	QString workspacePath = WorkspaceWizard::selectWorkspace(0, false);
	if (workspacePath.isEmpty())
		setEditTestEnabled(true);
	WorkspaceManager::instance().loadWorkspace(workspacePath.toStdString());
	ui.actionNewTest->setEnabled(true);
	ui.actionOpenTest->setEnabled(true);
	emit workspaceLoaded();
}
//__________________________________________________________________________________________________
/**
 * Slot to switch between tab window view and standard window view.
 */
void GetTest::setTabbedViewMode(bool enabled)
{
	centralZone->setViewMode(enabled?QMdiArea::TabbedView:QMdiArea::SubWindowView);
	ui.actionSetTabbedViewMode->setChecked(enabled);
	ui.actionTile_Windows->setDisabled(enabled);
	ui.actionCascaded_Windows->setDisabled(enabled);
}
//__________________________________________________________________________________________________
void GetTest::setSaveRawFramesEnabled(bool enabled)
{
	ControllerOptions::instance().setSaveRawFrameEnabled(enabled);
}
//__________________________________________________________________________________________________
void GetTest::setHeaderDisplayEnabled(bool enabled)
{
	ControllerOptions::instance().setHeaderDisplayEnabled(enabled);
}
//__________________________________________________________________________________________________
void GetTest::setMonitoringEnabled(bool enabled)
{
	ControllerOptions::instance().setDataDisplayEnabled(enabled);
}
//_________________________________________________________________________________________________
void GetTest::about()
{
	QMessageBox::about(this, tr("About GetController"),
		tr("<h2> GetController </h2>"
		"<br><br>"
		"&copy; 2010-2017 <b>CEA, <a href=%1>IRFU</a></b><br><a href=%2>CeCILL</a> license"
		"<br><br>"
		"For help, refer to the wiki of the GET collaboration at <a href=%3>%3</a>."
		).arg("http://irfu.cea.fr").arg("http://cecill.info/licences/Licence_CeCILL_V2-en.html").arg("https://dsm-trac.cea.fr/get"));
}
//__________________________________________________________________________________________________
/**
 * Asks user to select channels to display.
 */
void GetTest::selectChannels()
{
	bool ok;
	ChannelMask newMask = EditChannelMaskDialog::selectChannels(this, channels, &ok);
	if (ok)
	{
		channels = newMask;
	}
}
//__________________________________________________________________________________________________
/**
 * Opens a dialog to create a new test.
 */
void GetTest::createTest()
{
	EditTestDialog dialog(this);
	if (QDialog::Accepted == dialog.exec())
	{
		WorkspaceManager::instance().saveTest(dialog.getTestSummary());
		LOG_INFO() << "Created test '" << dialog.getTestSummary().name() << '\'';
		selectTest(dialog.getTestSummary().name());
	}
}
//__________________________________________________________________________________________________
/**
 * Edit current test.
 */
void GetTest::editTest()
{
	EditTestDialog dialog(this);
	try
	{
		std::auto_ptr< TestSummary > currentTest = WorkspaceManager::instance().currentTest();
		dialog.setTestSummary(*currentTest.get());
		LOG_INFO() << "Editing test '" << dialog.getTestSummary().name() << "'";
		if (QDialog::Accepted == dialog.exec())
		{
			WorkspaceManager::instance().saveTest(dialog.getTestSummary());
			emit testModified();
		}
	}
	catch (const GetControllerException & e)
	{
		LOG_ERROR() << e.what();
	}
}

void GetTest::selectTest()
{
	SelectTestDialog dialog(this);
	if (QDialog::Accepted == dialog.exec())
	{
		selectTest(dialog.selectedTest());
	}
}

void GetTest::selectTest(const std::string & name)
{
	if (name.empty()) return;
	WorkspaceManager::instance().setCurrentTest(name);
	setEditTestEnabled(true);
	std::auto_ptr< TestSummary > test = WorkspaceManager::instance().currentTest();
	QFileInfo harwareInfo(QString::fromStdString(test->hwConfigPath()));
	ui.actionEditHardware->setEnabled(harwareInfo.exists());
	ui.actionEditConditions->setEnabled(true);
	setControlStatus(Ctrl::Idle);
	emit testSelected();
	LOG_INFO() << "Selected test '" << name << '\'';
}

void GetTest::monitorInjection()
{
	LOG_INFO() << "Test monitor injection";
	if ((testController->iInject < 1000) and (daqrun == true)) QTimer::singleShot(1000, this, SLOT(monitorInjection()));
	LOG_INFO() << "Injection up to: " << testController->iInject;

}

void GetTest::beginTimer()
{
	QTimer::singleShot(1000, this, SLOT(endTimer()));
}

void GetTest::endTimer()
{
	emit endTime();
}

void GetTest::destroyTestController()
{
	if (testController.get())
	{
		testController->wait();
		testController.reset(0);
	}
}
//__________________________________________________________________________________________________
void GetTest::createTestController(const std::string & testId)
{
	// Delete previous instance
	destroyTestController();

	// Create new instance
	TestController* controller = new TestController(testId, *dataBloc.get(), this, this);
	testController.reset(controller);
	connect(controller, SIGNAL(finished()), this, SLOT(clearData()));
	connect(controller, SIGNAL(startTps()), this, SLOT(beginTimer()));
	connect(controller, SIGNAL(logHd(QString)), this, SIGNAL(logHd(QString)));
	connect(controller, SIGNAL(runFinished()), this, SLOT(onRunFinished()));
	connect(this, SIGNAL(stopAcq()), controller, SLOT(setStop()));
	connect(this, SIGNAL(endTime()), controller, SLOT(setEndTimer()));
	connect(this, SIGNAL(asicNum(int)), controller, SLOT(saveAsicNb(int)));
	verbose = testController->isVerbose();
}
//__________________________________________________________________________________________________
/**
 * Qt Slot to load hardware config file to ECC.
 */
void GetTest::loadHw()
{
	msgConsole->clear();

	setControlStatus(Ctrl::Off);

	try
	{
		// Get identifier of current test
		std::string testId = WorkspaceManager::instance().currentTestName();
		if (testId.empty())
		{
			LOG_ERROR() << "Could not load hardware configuration: no test has been selected!";
			return;
		}

		std::string hardwareConfigPath = WorkspaceManager::instance().currentTest()->hwConfigPath();
		LOG_INFO() << "Loading hardware config. from '" << hardwareConfigPath << '\'';
		QCoreApplication::processEvents();

		// Create test controller
		createTestController(testId);

		// Initialize variable for data synchronization between threads.
		daqrun = false;

		// Initialize test
		LOG_INFO() << "Initializing test '" << testId << '\'';
		testController->loadAndPrepare();

		// Connect board to data router
		MSG_IF(verbose) << "Waiting for connection to data router...";
		testController->daqConnect();

		setControlStatus(Ctrl::Loaded);
		return;
	}
	CATCH_ALL_EXCEPTIONS()

	setControlStatus(Ctrl::Idle);

}
//_________________________________________________________________________________________________
/**
 * Displays message in the status bar during 2 seconds.
 */
void GetTest::showStatusMessage(const QString & msg)
{
	statusBar()->showMessage(msg, 2000);
}
//_________________________________________________________________________________________________
/**
 * Creates widget for displaying data from AsAd board.
 * @param asadIdx Index of the AsAd board.
 */
void GetTest::createDisplay(const size_t & asadIdx)
{
	removeDisplay(asadIdx);
	LOG_DEBUG() << "Creating display widgets for AsAd " << asadIdx;
	// Create display for reconstructed channels
	displayTabWidgets[asadIdx] = new QTabWidget(this);
	QTabWidget* tabWidget = displayTabWidgets[asadIdx];
	tabWidget->setTabPosition(QTabWidget::South);
	for (size_t agetIdx=0; agetIdx < 4; agetIdx++)
	{
		displayGraphs[asadIdx*4+agetIdx] = new ViewGraph(tabWidget);
		ViewGraph* graph = displayGraphs[asadIdx*4+agetIdx];
		graph->setTitle(tr("AsAd no. %1 - AGET no. %2").arg(asadIdx).arg(agetIdx));
		graph->setAxisTitle(QwtPlot::xBottom, tr("Time Bucket Index"));
		connect(graph, SIGNAL(statusMessage(const QString &)), this, SLOT(showStatusMessage(const QString &)));
		tabWidget->addTab(graph, tr("AGET no. &%1").arg(agetIdx));
	}
	displayMdiSubWindows[asadIdx] = new QMdiSubWindow(this, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	QMdiSubWindow* subWindow = displayMdiSubWindows[asadIdx];
	subWindow->setWidget(tabWidget);
	subWindow->setAttribute(Qt::WA_Disabled);
	subWindow->setWindowTitle(tr("AsAd no. %1").arg(asadIdx));
	subWindow = centralZone->addSubWindow(subWindow, Qt::CustomizeWindowHint | Qt::WindowTitleHint | Qt::WindowMinimizeButtonHint | Qt::WindowMaximizeButtonHint);
	subWindow->setEnabled(true);
	subWindow->setVisible(true);
	centralZone->setActiveSubWindow(subWindow);
}
//_________________________________________________________________________________________________
/**
 * Deletes widget for displaying data from AsAd board.
 * @param asadIdx Index of the AsAd board.
 */
void GetTest::removeDisplay(const size_t & asadIdx)
{
	if (displayMdiSubWindows.count(asadIdx) != 0)
	{
		LOG_DEBUG() << "Deleting display for AsAd " << asadIdx;
		centralZone->setActiveSubWindow(consWin);
		QMdiSubWindow* subWindow = displayMdiSubWindows[asadIdx];
		centralZone->removeSubWindow(subWindow);
		delete subWindow;
		displayMdiSubWindows.erase(asadIdx);
		// QMdiSubWindow takes temporary ownership of widget; you do not have to delete it.
		displayTabWidgets.erase(asadIdx);
		updateSubWindowMenu();
	}
}
//_________________________________________________________________________________________________
/**
 * Creates widget for displaying data from all 4 AsAd boards.
 */
void GetTest::createDisplays()
{
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		createDisplay(asadIdx);
	}
	centralZone->setActiveSubWindow(displayMdiSubWindows[0]);
}
//_________________________________________________________________________________________________
/**
 * Creates widget for displaying data from certain AsAd board(s).
 * @param Mask of AsAd boards.
 */
void GetTest::createDisplays(const uint16_t & asadMask)
{
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		if (::utl::BitFieldHelper< uint16_t>::getBit(asadMask, asadIdx))
		{
			createDisplay(asadIdx);
			centralZone->setActiveSubWindow(displayMdiSubWindows[asadIdx]);
		}
	}
}
//_________________________________________________________________________________________________
/**
 * Deletes all data displays.
 */
void GetTest::removeDisplays()
{
	for (size_t asadIdx=0; asadIdx < 4u; ++asadIdx)
	{
		removeDisplay(asadIdx);
	}
}
//_________________________________________________________________________________________________
/**
 * Resets DAQ server and clears all data from queue.
 */
void GetTest::clearData()
{
	//resetDataRouter();
	dataBloc->clear();
}
//__________________________________________________________________________________________________
/**
 * Qt Slot to configure the GET Hardware system via the ECC.
 **/
void GetTest::configure()
{
	// Get Identifier of current test
	std::string testId = WorkspaceManager::instance().currentTestName();
	if (testId.empty())
	{
		LOG_ERROR() << "Could not load hardware configuration: no test has been selected!";
		return;
	}

	setControlStatus(Ctrl::Off);

	try
	{
		// Reset DAQ server
		clearData();
		// Remove pre-existing graphs
		removeDisplays();

		// (Re-)create test controller
		createTestController(testId);

		emit logHd(QString("Display of headers is %1 enabled.\n").arg(ControllerOptions::instance().isHeaderDisplayEnabled()?"":"not"));

		std::ostringstream oss;
		oss << "Channels to display: " << std::dec << channels.toString();
		LOG_DEBUG() << oss.str();

		// Configure board(s)
		testController->configure();

		setControlStatus(Ctrl::Configured);
		return;
	}
	CATCH_ALL_EXCEPTIONS()

	setControlStatus(Ctrl::Loaded);
}
//_________________________________________________________________________________________________
/**
 * Qt Slot to start acquisition.
 **/
void GetTest::startDaq()
{

	try
	{
		uicW.frameCountBox->setValue(0);
		setControlStatus(Ctrl::Running);
		clearData();

		if (ControllerOptions::instance().isDataDisplayEnabled() and not ControllerOptions::instance().isDataRouterExternal())
		{
			createDisplays(testController->asadMask());
		}
		// Check whether to print frame headers
		isHeaderDisplayEnabled_ = ControllerOptions::instance().isHeaderDisplayEnabled();

		if (testController->testSummary().type() == std::string("InjectTest"))
			QTimer::singleShot(1000, this, SLOT(monitorInjection()));

		LOG_INFO() << "Starting test '" << testController->testSummary().name() << '\'';
		testController->start();
		daqrun = true;
	}
	CATCH_ALL_EXCEPTIONS()
}
//_________________________________________________________________________________________________
/**
 * Qt Slot to stop acquisition
 **/
void GetTest::stopDaq()
{
	setControlStatus(Ctrl::Off);
	try
	{
		if (daqrun)
		{
			daqrun = false;
			MSG_IF(verbose) << "Stopping data acquisition. Wait...";
			emit stopAcq();
		}
		while (not testController->isRunFinished())
		{
			std::cout << "Waiting for test controller to finish... \r" << std::flush;
			QCoreApplication::processEvents();
		}
	}
	CATCH_ALL_EXCEPTIONS()
	setControlStatus(Ctrl::Configured);
}

void GetTest::onRunFinished()
{
	setControlStatus(Ctrl::Configured);
}

void GetTest::setControlStatus(Ctrl::State state)
{
	ctrlStatus_ = state;
	switch (state)
	{
		case Ctrl::Off:
			uicW.pbLoadHw->setEnabled(false);
			uicW.pbConfig->setEnabled(false);
			uicW.pbLaunch->setEnabled(false);
			uicW.pbStop->setEnabled(false);
			uicW.pbReset->setEnabled(false);
			setEditTestEnabled(false);
			break;
		case Ctrl::Loaded:
			uicW.pbLoadHw->setEnabled(false);
			uicW.pbConfig->setEnabled(true);
			uicW.pbLaunch->setEnabled(false);
			uicW.pbStop->setEnabled(false);
			uicW.pbReset->setEnabled(true);
			setEditTestEnabled(true);
			break;
		case Ctrl::Configured:
			uicW.pbLoadHw->setEnabled(false);
			uicW.pbConfig->setEnabled(false);
			uicW.pbLaunch->setEnabled(true);
			uicW.pbStop->setEnabled(false);
			uicW.pbReset->setEnabled(true);
			setEditTestEnabled(false);
			break;
		case Ctrl::Running:
			uicW.pbLoadHw->setEnabled(false);
			uicW.pbConfig->setEnabled(false);
			uicW.pbLaunch->setEnabled(false);
			uicW.pbStop->setEnabled(true);
			uicW.pbReset->setEnabled(false);
			setEditTestEnabled(false);
			break;
		case Ctrl::Idle:
		default:
			uicW.pbLoadHw->setEnabled(true);
			uicW.pbConfig->setEnabled(false);
			uicW.pbLaunch->setEnabled(false);
			uicW.pbStop->setEnabled(false);
			uicW.pbReset->setEnabled(false);
			setEditTestEnabled(true);
			break;
	}
}

/**
 * Qt Slot to reset the hardware configuration
 **/
void GetTest::Reset()
{
	//errorConsole->clear();
	uicW.frameCountBox->setValue(0);
	switch(controlStatus())
	{
	case Ctrl::Loaded:
		try
		{
			if (testController.get())
				testController->unprepare();
			resetRunProcessor();
		}
		CATCH_ALL_EXCEPTIONS()
		setControlStatus(Ctrl::Idle);
		break;
	case Ctrl::Configured:
		setControlStatus(Ctrl::Loaded);
		break;
	default:
		daqrun = false;
		try
		{
			stopDaq();
			clearData();
			setControlStatus(Ctrl::Idle);
		}
		CATCH_ALL_EXCEPTIONS()
		destroyTestController();
		break;
	}
}

/**
 * Updates GUI after the reception of a frame.
 */
void GetTest::onFrameReceived()
{
	uicW.frameCountBox->setValue(uicW.frameCountBox->value() + 1);
}
//_________________________________________________________________________________________________
/**
 * Selects symbol for curves.
 */
void GetTest::selectCurveSymbol()
{
	QString currentStyle = curveSymbolStyle;
	QString newStyle = gui::QwtSymbolDialog::getSymbolStyle(currentStyle, this);
	if (not newStyle.isEmpty() and newStyle != currentStyle)
	{
		QSettings().setValue("symbolStyle", newStyle);
		curveSymbolStyle = newStyle;
	}
}
//_________________________________________________________________________________________________
/**
 * Displays event data channel by channel in a Qwt window.
 */
void GetTest::displayEvent(const CoBoEvent & event)
{
	if (not ControllerOptions::instance().isDataDisplayEnabled()) return;
	try
	{
		const size_t asadIdx = event.asadIdx();
		if (not displayTabWidgets[asadIdx]) return;
		typedef std::set< short > ChannelSet;

		LOG_INFO() << "Displaying event no. " << event.eventIdx();
		// Set contents of Qwt graph
		for (size_t agetIdx=0; agetIdx < 4u; ++agetIdx)
		{
			ViewGraph* graph = displayGraphs[asadIdx*4+agetIdx];
			graph->clear();
			if (not event.containsAGetChip(agetIdx)) continue;

			ChannelSet::const_iterator chanIter;
			for (chanIter = channels.begin(); chanIter != channels.end(); ++chanIter)
			{
				try
				{
					short chanIdx = *chanIter;
					std::vector< double > xData(event.channel(chanIdx, agetIdx).buckIndexes().begin(),
							event.channel(chanIdx, agetIdx).buckIndexes().end());
					// Add readout offset
					if (event.readOffset() > 0)
					{
						transform(xData.begin(), xData.end(), xData.begin(), bind2nd(std::plus< double > (), event.readOffset()));
					}
					const std::vector< double > yData(event.channel(chanIdx, agetIdx).sampleValues().begin(),
						event.channel(chanIdx, agetIdx).sampleValues().end());
					graph->addCurve(QObject::tr("%1").arg(chanIdx), curveSymbolStyle, true,
							event.channel(chanIdx, agetIdx).sampleCount(), &xData.front(), &yData.front());
					graph->updateAxes();
					graph->validatePlot();
				}
				catch (const CoBoEvent::ChannelNotFound & e)
				{
					; //LOG_WARN() << e.what();
				}
			}
			graph->setTitle(tr("Event no. %1 - AsAd no. %2 - AGET no. %3").arg(event.eventIdx()).arg(asadIdx).arg(agetIdx));
		}
		busyDisplaying[asadIdx] = false;
	}
	CATCH_ALL_EXCEPTIONS()
}
//_________________________________________________________________________________________________
/**
 * Prints header.
 */
void GetTest::displayHeader(mfm::Frame & frame)
{
	// Print frame header
	if (not isHeaderDisplayEnabled_) return;
	std::ostringstream oss;
	TestController::writeHeader(frame, oss);
	emit logHd(QString::fromStdString(oss.str()));
}
//_________________________________________________________________________________________________
/**
 * Destructor.
 */
GetTest::~GetTest()
{
    // Set logging backend
	::utl::BackendLogger::setBackend(::utl::LoggingBackendPtr(new ::mdaq::utl::ConsoleLoggingBackend));
	destroyTestController();
}

#undef CATCH_ALL_EXCEPTIONS
