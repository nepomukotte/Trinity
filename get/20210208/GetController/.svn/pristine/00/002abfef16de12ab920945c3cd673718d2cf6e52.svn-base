/* =============================================================================
 * GetTest.h
 * -----------------------------------------------------------------------------
 * class get::GetTest
 * Created on: 24 August 2010 at Irfu/Sedi/Ldef, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Druillole (frederic.druillole@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the GetController project dedicated to the GET
 * data acquisition test bench. "GET" stands for General Electronics for TPC.
 * http://www-actar-get.cea.fr
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */
#ifndef get_GetTest_h_INCLUDED
#define get_GetTest_h_INCLUDED

/**
 * \file GetTest.h
 * \brief Main class to create the main GUI of GetController project
 * \author Frederic Druillole (frederic.druillole@cea.fr)
 * \version 0.1
 * \date 24 August 2010
 *
 * Class contains all graphical interface management for the GetController application
 *
 */
#include "get/GetTest.ui.h"
#include "get/ControlWidget.ui.h"
#include "get/DataBloc.h"
#include "get/plot.h"
#include "get/viewgraph.h"
#include "get/ChannelMask.h"
#include "get/CoBoEvent.h"
#include "get/TestController.h"

#include "mdaq/DefaultPortNums.h"
#include "mdaq/daq/RunProcessor.h"
#include "utl/net/SocketAddress.h"

#include <QMainWindow>
#include <QVector>
#include <QTime>
#include <IceUtil/Handle.h>
#include <iostream>
#include <memory>
#include <map>

// Forward declarations
class Plot;
class Console;
namespace get
{
	class CoBoEvent;
}
namespace mfm
{
	class Frame;
}

class QDockWidget;
class QMdiArea;
class QMdiSubWindow;
class QSignalMapper;
class QTextEdit;
class QTreeView;
class QwtSymbol;

namespace Ctrl
{
	enum State { Off=0, Idle, Loaded, Configured, Running};
}

namespace get
{
//__________________________________________________________________________________________________
/**
 * A class to create GUI and manage the GetController application.
 */
class GetTest: public QMainWindow, public IFrameProcessor
{
	Q_OBJECT /*!< Belong to Qt signal/slot interface*/
public:
	GetTest(QWidget* parentWidget=0);
	~GetTest();
	bool createDataRouter();
	void destroyDataRouter();
	bool createAlarmLogger(const ::utl::net::SocketAddress & address);
	void startRunProcessor();
	void stopRunProcessor();
	void resetRunProcessor();
	void frameAvailable(mfm::Frame& getFrame); ///< @override IFrameProcessor::frameAvailable
signals:
	void logHd(const QString &);
	void stopAcq();
	void requestFrameMon();
	void requestEventMon();
	void isAgetReady();
	void endTime();
	void startT();
	void asicNum(int);
	void testSelected();
	void testModified();
	void workspaceLoaded();
	void dataRouterSettingsChanged();
	void editTestEnabled(bool enabled);
	void eventReady(const get::CoBoEvent&);

public slots:
	/// Prompts the user for data router parameters and creates data router if necessary
	bool selectDataRouter();
private slots:
	void about();
	void onRunFinished();
	void onFrameReceived();
	void cleanUp();
	void selectCurveSymbol();
	void editConditionsCConfig();
	void editHardwareCConfig();
	void switchWorkspace();
	void setTabbedViewMode(bool enabled);
	void updateSubWindowMenu();
	void setActiveSubWindow(QWidget *window);
	void createTest();
	void editTest();
	void selectTest();
	void selectTest(const std::string & name);
	void setSaveRawFramesEnabled(bool enabled);
	void selectChannels();
	void setHeaderDisplayEnabled(bool enabled);
	void setMonitoringEnabled(bool enabled);
	void displayEvent(const get::CoBoEvent & coboEvent);
	void displayHeader(mfm::Frame & frame);
	void createDisplay(const size_t & asadIdx);
	void removeDisplay(const size_t & asadIdx);
	void createDisplays();
	void createDisplays(const uint16_t & asadMask);
	void removeDisplays();
	void loadHw();
	void clearData();
	void configure();
	void startDaq();
	void stopDaq();
	void Reset();
	void monitorInjection();
	void beginTimer();
	void endTimer();
	void showStatusMessage(const QString & msg);
	void showMessages();
	void showErrors();
	void setEditTestEnabled(bool enabled);
private:
	mdaq::daq::RunProcessorPrx runProcessor();
	// void createThreads();
	void destroyTestController();
	void createTestController(const std::string & testId);
	void createActions();
	void connectActions();
	void setControlStatus(Ctrl::State state);
	Ctrl::State controlStatus() const { return ctrlStatus_; }
private:
	::Ui::GetController ui; /*!<GUI Qtdesigner ressources mainWindow*/
	::Ui::ControlWidgetClass uicW; /*!<GUI Qtdesigner ressources Control buttons*/
	QMdiArea *centralZone; ///< Multiple Document Interface area.
	QMenu* subWindowMenu;
	QSignalMapper* windowMapper;
	QMdiSubWindow *consWin; ///< Standard output console window.
	Console *msgConsole;
	//Console *errorConsole;
	QDockWidget* controlDock;
	QDockWidget* testInfoDock;
	std::map< size_t, QMdiSubWindow* > displayMdiSubWindows; ///< MDI sub-windows for the data display of the 4 AsAd boards
	std::map< size_t, QTabWidget* > displayTabWidgets; ///< Tab widgets for the data display of the 4 AsAd boards
	std::map< size_t, ViewGraph* > displayGraphs; ///< Display graphs
	//QMdiSubWindow *errorLog;
	QMdiSubWindow *HeaderLog;
	QTextEdit *zT1;
	QTreeView *vue;
	std::auto_ptr<DataBloc> dataBloc;
	bool daqrun;
	std::auto_ptr<TestController> testController;
	double agetData[34816];
	QVector<double> XOscillo;
	QVector<double> YOscillo;
	bool verbose;
	bool isDataTrash;
	ChannelMask channels; ///< Channels to display
	/// @name Data router
	/// @{
	/// @}
	Ctrl::State ctrlStatus_;
	QTime displayTime[4];
	CoBoEvent coboEvent[4]; /// Events to display
	bool busyDisplaying[4];
	QString curveSymbolStyle;
	bool isHeaderDisplayEnabled_; ///< Whether to print frame headers
};
//__________________________________________________________________________________________________
} // namespace get

#endif // get_GetTest_h_INCLUDED
