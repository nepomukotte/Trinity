/*
 * $Id: CoBoFrameViewer.cpp 1667 2017-07-06 12:58:35Z psizun $
 * @file CoBoFrameViewer.cpp
 * @created 14 juin 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: patrick.sizun@cea.fr
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of
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

#include "get/CoBoFrameViewer.h"
#include "get/viewgraph.h"
#include "get/EditChannelMaskDialog.h"
#include "get/CoBoStandardDisplay.h"
#ifdef HAVE_QWT_PLOT_3D
#include "get/CoBoSurfDisplay.h"
#endif
#include "QwtSymbolDialog.h"

#include <mfm/Exception.h>
#include <mfm/Frame.h>
#include <mfm/Field.h>
#include <mfm/BitField.h>
#include <utl/Logging.h>
using namespace mfm;

#include <QAction>
#include <QApplication>
#include <QComboBox>
#include <QMenu>
#include <QMenuBar>
#include <QStatusBar>
#include <QCoreApplication>
#include <QTabWidget>
#include <QToolBar>
#include <QMessageBox>
#include <QFile>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <QPainter>
#include <QPixmap>
#include <QSettings>
#include <QSignalMapper>
#include <QSpinBox>
#include <QInputDialog>
#include <QTimer>
#include <QLabel>
#include <QHBoxLayout>
#include <QtGlobal>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
//#include <boost/timer.hpp>

#include <cerrno>
using std::strerror;
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

class IntSeqGenerator
{
private:
	int value;
public:
	IntSeqGenerator(int initialValue = 0) : value(initialValue)
	{
	}
	int operator()()
	{
		return value++;
	}
};

namespace get {
//_________________________________________________________________________________________________
CoBoFrameViewer::CoBoFrameViewer(QWidget* parentWidget)
	: QMainWindow(parentWidget), autoscanTimer(0), scanPeriod_ms(2000), mode3d_(false), numRawFrames_(1u)
{
	// Central widget
	graphTabWidget = new QTabWidget(this);
	graphTabWidget->setTabPosition(QTabWidget::South);
#if QT_VERSION >= 0x040200
	graphTabWidget->setObjectName("coboTabWidget");
	//graphTabWidget->setStyleSheet("background-image: url(:/get/logo.png)");
	//graphTabWidget->setStyleSheet("border-image: url(:/get/logo.png)");
#endif
	// Actions
	createActions();

	// Default channels
	IntSeqGenerator gen(0);
	std::generate_n(std::inserter(channels, channels.begin()), 68, gen);

	connect(this, SIGNAL(frameChanged()), this, SLOT(onFrameChanged()));

	setCentralWidget(graphTabWidget);
	setWindowTitle(tr("CoBo Frame Viewer"));
	setWindowIcon(QPixmap(":/get/logo-highres.png"));
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::createDisplays(bool mode3d)
{
	// Remember current tab
	int currentDisplayIndex = 0;
	if (graphTabWidget->count() > 0)
		currentDisplayIndex = graphTabWidget->currentIndex();


	// Delete previous widgets
	while (graphTabWidget->count() > 0)
	{
		QWidget *display = graphTabWidget->widget(0);
		graphTabWidget->removeTab(0);
		if (display)
			delete display;
	}
	displays.clear();

	// Create displays
	for (size_t frameId=0; frameId < numRawFrames_; ++frameId)
	{
		for (size_t agetIdx=0; agetIdx < 4; agetIdx++)
		{
			QString displayLabel = (numRawFrames_ > 1)?tr("%1").arg(frameId*4+agetIdx):tr("AGET no. %1").arg(agetIdx);
#ifdef HAVE_QWT_PLOT_3D
			if (not mode3d)
			{
#endif
				CoBoStandardDisplay* display = new CoBoStandardDisplay(this);
				connect(display, SIGNAL(statusMessage(const QString &)), this, SLOT(showStatusMessage(const QString &)));
				graphTabWidget->addTab(display, displayLabel);
				displays[frameId*4+agetIdx] = display;
#ifdef HAVE_QWT_PLOT_3D
			}
			else
			{
				CoBoSurfDisplay* display = new CoBoSurfDisplay(this);
				graphTabWidget->addTab(display, displayLabel);
				displays[frameId*4+agetIdx] = display;
			}
#endif
		}
	}
	selectCurrentDisplay(currentDisplayIndex);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::createActions()
{
	// Menu
	QMenu* menuFile_ = menuBar()->addMenu(tr("&File"));
	QMenu* menuNavigate_ = menuBar()->addMenu(tr("&Navigate"));
	QMenu* menuDisplay_ = menuBar()->addMenu(tr("&Display"));
	QMenu* menuAnalysis_ = menuBar()->addMenu(tr("&Analysis"));
	QMenu* menuHelp_ = menuBar()->addMenu(tr("&Help"));

	// Toolbar
	QToolBar* toolBar_ = addToolBar(tr("File"));

	// Actions
	QAction* actionOpen = new QAction(tr("Open"), this);
	actionOpen->setShortcuts(QKeySequence::Open);
#if QT_VERSION >= 0x040600
	actionOpen->setIcon(QIcon::fromTheme("document-open"));
#endif
	connect(actionOpen, SIGNAL(triggered()), this, SLOT(openFile()));
	menuFile_->addAction(actionOpen);
	toolBar_->addAction(actionOpen);

	QAction* actionSaveAll = new QAction(tr("Save All As"), this);
#if QT_VERSION >= 0x040600
	actionSaveAll->setIcon(QIcon::fromTheme("document-save-as", QIcon(":/icons/oxygen/16x16/actions/document-save-all.png")));
#else
	actionSaveAll->setIcon(QIcon(":/icons/oxygen/16x16/actions/document-save-all.png"));
#endif
	connect(actionSaveAll, SIGNAL(triggered()), this, SLOT(saveAll()));
	menuFile_->addAction(actionSaveAll);
	toolBar_->addAction(actionSaveAll);

	actionPrevious = new QAction(tr("Previous"), this);
	actionPrevious->setDisabled(true);
	actionPrevious->setShortcuts(QKeySequence::FindPrevious);
	#if QT_VERSION >= 0x040600
	actionPrevious->setIcon(QIcon::fromTheme("go-previous"));
	#endif
	connect(actionPrevious, SIGNAL(triggered()), this, SLOT(selectPreviousFrame()));
	menuNavigate_->addAction(actionPrevious);
	toolBar_->addAction(actionPrevious);

	actionNext = new QAction(tr("Next"), this);
	actionNext->setDisabled(true);
	actionNext->setShortcuts(QKeySequence::FindNext);
	#if QT_VERSION >= 0x040600
	actionNext->setIcon(QIcon::fromTheme("go-next"));
	#endif
	connect(actionNext, SIGNAL(triggered()), this, SLOT(selectNextFrame()));
	menuNavigate_->addAction(actionNext);
	toolBar_->addAction(actionNext);

	actionGoTo = new QAction(tr("Go To Frame..."), this);
	actionGoTo->setDisabled(true);
	actionGoTo->setShortcuts(QKeySequence::Find);
	#if QT_VERSION >= 0x040600
	actionGoTo->setIcon(QIcon::fromTheme("go-jump"));
	#endif
	connect(actionGoTo, SIGNAL(triggered()), this, SLOT(selectFrame()));
	menuNavigate_->addAction(actionGoTo);
	toolBar_->addAction(actionGoTo);

	actionPlay = new QAction(tr("Autoscan"), this);
	actionPlay->setDisabled(true);
	actionPlay->setCheckable(true);
	actionPlay->setShortcut(QKeySequence("Space"));
#if QT_VERSION >= 0x040600
	QIcon playIcon;
	playIcon.addPixmap(QIcon::fromTheme("media-playback-start").pixmap(16), QIcon::Normal, QIcon::Off);
	playIcon.addPixmap(QIcon::fromTheme("media-playback-stop").pixmap(16), QIcon::Normal, QIcon::On);
	actionPlay->setIcon(playIcon);
#endif
	autoscanTimer = new QTimer(this);
	QObject::connect(autoscanTimer, SIGNAL(timeout()), this, SLOT(updateAutoscan()));
	connect(actionPlay, SIGNAL(toggled(bool)), this, SLOT(toggleAutoscan(bool)));
	menuNavigate_->addAction(actionPlay);
	toolBar_->addAction(actionPlay);
	graphTabWidget->addAction(actionPlay);

	QSettings settings;
	scanPeriod_ms = settings.value("scanPeriod_ms", 2000).toInt();
	actionEditPeriod = new QAction(tr("Edit Period"), this);
	connect(actionEditPeriod, SIGNAL(triggered()), this, SLOT(editPeriod()));
	menuNavigate_->addAction(actionEditPeriod);

	action3D = new QAction(tr("3D"), this);
	action3D->setCheckable(true);
	action3D->setShortcut(QKeySequence("D"));
#ifdef HAVE_QWT_PLOT_3D
	connect(action3D, SIGNAL(toggled(bool)), this, SLOT(toggleDisplay3d(bool)));
#else
	action3D->setDisabled(true);
#endif
	menuDisplay_->addAction(action3D);
	toolBar_->addAction(action3D);

	actionSelectChannels = new QAction(tr("Select Channels"), this);
	actionSelectChannels->setShortcuts(QKeySequence::Preferences);
	#if QT_VERSION >= 0x040600
	actionSelectChannels->setIcon(QIcon::fromTheme("document-properties"));
	#endif
	connect(actionSelectChannels, SIGNAL(triggered()), this, SLOT(selectChannels()));
	toolBar_->addAction(actionSelectChannels);

	QAction* actionQuit = menuFile_->addAction(tr("Quit"));
	actionQuit->setShortcuts(QKeySequence::Quit);
#if QT_VERSION >= 0x040200
	actionQuit->setMenuRole(QAction::QuitRole);
#endif
#if QT_VERSION >= 0x040600
	actionQuit->setIcon(QIcon::fromTheme("application-exit"));
#endif
	connect(actionQuit, SIGNAL(triggered()), qApp, SLOT(quit()));

	// Full screen mode
	QAction* actionFullScreen = new QAction(tr("&Fullscreen"), this);
	actionFullScreen->setCheckable(true);
	actionFullScreen->setShortcut(QKeySequence(Qt::Key_F));
	connect(actionFullScreen, SIGNAL(toggled(bool)), this, SLOT(showFullscreen(bool)));
	menuDisplay_->addAction(actionFullScreen);
	graphTabWidget->addAction(actionFullScreen);

	// Symbol style
	QAction* actionSelectSymbol = new QAction(tr("&Symbol style"), this);
	connect(actionSelectSymbol, SIGNAL(triggered()), this, SLOT(selectSymbol()));
	menuDisplay_->addAction(actionSelectSymbol);

	// Color of line and symbol
	QAction* actionSelectColorMode = new QAction(tr("&Colors..."), this);
	connect(actionSelectColorMode, SIGNAL(triggered()), this, SLOT(selectRandomColorMode()));
	menuDisplay_->addAction(actionSelectColorMode);

	// FWHM
	actionCalcFwhm = new QAction(tr("&FHWM"), this);
	actionCalcFwhm->setToolTip(tr("Compute Full-Width at Half-Maximum"));
	connect(actionCalcFwhm, SIGNAL(triggered()), this, SLOT(calcFwhm()));
	menuAnalysis_->addAction(actionCalcFwhm);

	// Help
	QAction* actionAbout = menuHelp_->addAction(tr("About"));
	actionAbout->setShortcuts(QKeySequence::HelpContents);
#if QT_VERSION >= 0x040200
	actionAbout->setMenuRole(QAction::AboutRole);
#endif
#if QT_VERSION >= 0x040600
	actionAbout->setIcon(QIcon::fromTheme("help-about"));
#endif
	connect(actionAbout, SIGNAL(triggered()), this, SLOT(about()));
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::about()
{
	QMessageBox::about(this, tr("About CoBoFrameViewer"),
		tr("<h2> CoBoFrameViewer </h2>"
		"<br><br>"
		"<b>CoBoFrameViewer can be used to display the contents of files containing raw CoBo data frames in the MFM format."
		"<br><br>"
		"&copy; 2012-2014 <b>CEA, <a href=%1>IRFU</a></b><br><a href=%2>CeCILL</a> license"
		"<br><br>"
		"For further details, refer to the wiki of the GET collaboration at <a href=%3>%3</a>."
		).arg("http://irfu.cea.fr").arg("http://www.cecill.info").arg("https://dsm-trac.cea.fr/get"));
}
//_________________________________________________________________________________________________
QSize CoBoFrameViewer::sizeHint() const
{
	return QSize(1280, 800);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::onFrameChanged()
{
	actionPrevious->setEnabled(frameIndex_ > 0u);
	actionNext->setEnabled(true);
	actionGoTo->setEnabled(true);
	actionPlay->setEnabled(true);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::showStatusMessage(const QString & msg)
{
	statusBar()->showMessage(msg);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectCurrentDisplay(const size_t idx)
{
	graphTabWidget->setCurrentIndex(idx);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::setChannels(const ChannelMask & channels)
{
	this->channels = channels;
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::openFile()
{
	QString newFilename = QFileDialog::getOpenFileName(this, tr("Select CoBo binary frame file"), QString(),
			QString("%1 (*.graw *.gmer);;%2 (*.gmat);;%3 (*)")
			.arg(tr("Raw frame files")).arg(tr("Matrix files")).arg(tr("All files")));

	if (newFilename.isEmpty())
			return;

	if (not setFrame(newFilename.toStdString(), 0u))
	{
		QMessageBox::warning(this, tr("Error"), tr("Error loading data from file %1.").arg(newFilename));
	}
}
//_________________________________________________________________________________________________
/**
 * Saves all displays into a single image.
 */
void CoBoFrameViewer::saveAll()
{
	QString defaultBaseName  = QString("%1_frame%2")
			.arg(QFileInfo(QString::fromStdString(filename_)).completeBaseName())
			.arg(frameIndex_);
	QString defaultFileName = defaultBaseName + ".png";

	// Get filename
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save All Tabs As"), defaultFileName);
	if (fileName.isEmpty()) return;

	// Create individual images
	bool ok = true;
	QPixmap images[4];
	for (int agetIdx=0; agetIdx < 4 and ok; agetIdx++)
	{
		QString agetFileName = QString("%1_aget%2.%3").arg(QFileInfo(fileName).completeBaseName())
				.arg(agetIdx).arg(QFileInfo(fileName).suffix());
		ok = displays[agetIdx]->saveAs(agetFileName);
		if (ok) images[agetIdx] = QPixmap(agetFileName);
		// Remove file
		QFile::remove(agetFileName);
	}
	if (not ok) return;

	// Join images
	const int firstColWidth = qMax(images[0].width(), images[2].width());
	const int firstRowHeight = qMax(images[0].height(), images[1].height());
	const int secondColWidth = qMax(images[1].width(), images[3].width());
	const int secondRowHeight = qMax(images[2].height(), images[3].height());
	QPixmap mosaic(firstColWidth + secondColWidth, firstRowHeight + secondRowHeight);
	QPainter painter(&mosaic);
	painter.drawPixmap(QPoint(0, 0), images[0]);
	painter.drawPixmap(QPoint(firstColWidth, 0), images[1]);
	painter.drawPixmap(QPoint(0, firstRowHeight), images[2]);
	painter.drawPixmap(QPoint(firstColWidth, firstRowHeight), images[3]);

	// Save mosaic
	LOG_INFO() << "Saving mosaic as '" << fileName.toStdString() << '\'';
	mosaic.save(fileName);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectFrame()
{
	bool ok;
#if QT_VERSION >= 0x040500
	int newIndex = QInputDialog::getInt(this, tr("Go to frame ..."), tr("Frame no.:"), frameIndex_, 0, 1000000, 1, &ok);
#else
	int newIndex = QInputDialog::getInteger(this, tr("Go to frame ..."), tr("Frame no.:"), frameIndex_, 0, 1000000, 1, &ok);
#endif

	if (not ok) return;

	if (not setFrame(filename_, newIndex))
	{
		QMessageBox::warning(this, tr("Error"), tr("Error loading frame no. %1 from file %2.").arg(newIndex).arg(QString::fromStdString(filename_)));
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::toggleAutoscan(bool enable)
{
	if (enable)
	{
		LOG_DEBUG() << "Starting autoscan ...";
		actionEditPeriod->setDisabled(true);
		autoscanTimer->start(scanPeriod_ms);
	}
	else
	{
		LOG_DEBUG() << "Stopping autoscan ...";
		autoscanTimer->stop();
		actionEditPeriod->setEnabled(true);
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::startAutoscan(bool enable)
{
	actionPlay->setChecked(enable);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::editPeriod()
{
#if QT_VERSION >= 0x040500
	scanPeriod_ms = QInputDialog::getInt(this, tr("Edit scan period"),
			tr("Period [ms]:"), scanPeriod_ms, 1, 300000, 200);
#else
	scanPeriod_ms = QInputDialog::getInteger(this, tr("Edit scan period"),
			tr("Period [ms]:"), scanPeriod_ms, 1, 300000, 200);
#endif

	QSettings settings;
	settings.setValue("scanPeriod_ms", scanPeriod_ms);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::set3dDisplayMode(bool enable)
{
	action3D->setChecked(enable);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::toggleDisplay3d(bool enable)
{
	createDisplays(enable);
	reloadFrame();
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::updateAutoscan()
{
	if (not setFrame(filename_, frameIndex_ + 1))
	{
		if (not setFrame(filename_, 0u))
		{
			actionPlay->setChecked(false);
		}
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::showFullscreen(bool enable)
{
	if (enable)
	{
		graphTabWidget->setParent(0);
		graphTabWidget->setWindowFlags(graphTabWidget->windowFlags() | Qt::CustomizeWindowHint | Qt::WindowStaysOnTopHint | Qt::WindowMaximizeButtonHint | Qt::WindowCloseButtonHint);
		graphTabWidget->setWindowState(graphTabWidget->windowState() | Qt::WindowFullScreen);
		graphTabWidget->show();
	}
	else
	{
		//graphTabWidget->setWindowState(Qt::WindowNoState | Qt::WindowActive);
		graphTabWidget->showNormal();
		graphTabWidget->setParent(this);
		setCentralWidget(graphTabWidget);
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::reloadFrame()
{
	if (not filename_.empty())
	{
		setFrame(filename_, frameIndex_);
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectNextFrame()
{
	if (not setFrame(filename_, frameIndex_ + 1))
	{
		QMessageBox::warning(this, tr("Error"), tr("Error loading next frame from file %2.").arg(QString::fromStdString(filename_)));
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectSymbol()
{
	QString currentStyle = QSettings().value("symbolStyle", "NoSymbol").toString();
	QString newStyle = gui::QwtSymbolDialog::getSymbolStyle(currentStyle, this);
	if (not newStyle.isEmpty() and newStyle != currentStyle)
	{
		QSettings().setValue("symbolStyle", newStyle);
		createDisplays(mode3d_);
		reloadFrame();
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectRandomColorMode()
{
	bool currentColorMode = QSettings().value("randomColorsEnabled", false).toBool();
	bool colorMode = QMessageBox::Yes == QMessageBox::question(this, tr("Select color mode"),
			tr("Do you want to use random colors?"), QMessageBox::Yes | QMessageBox::No,
			currentColorMode?QMessageBox::Yes:QMessageBox::No);

	if (colorMode != currentColorMode)
	{
		QSettings().setValue("randomColorsEnabled", colorMode);
		createDisplays(mode3d_);
		reloadFrame();
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::calcFwhm()
{
	const get::CoBoEvent & event = coboEvent;
	//const size_t asadIdx = coboEvent.asadIdx();
	try
	{
		ChannelMap::const_iterator chanIt;
		for (chanIt = event.channels().begin(); chanIt != event.channels().end(); ++chanIt)
		{
			const uint16_t  aget = (uint16_t) chanIt->first.second;
			//const size_t fpnChannel = chanIt->second.fpnNeighbourChannel();
			const double offset = 0; //event.channel(fpnChannel, aget).meanSampleValue();
			LOG_INFO() << "aget=" << aget
					<< " chan=" << (uint16_t) chanIt->first.first
					//<< " fpnNeighbour=" << fpnChannel
					<< " FWHM=" << chanIt->second.fwhm(0.1, 0.5, offset)
					<< " peakingTime=" << chanIt->second.peakingTime(0.1, 0.05, offset);
		}
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectPreviousFrame()
{
	if (frameIndex_ <= 0u) return;
	if (not setFrame(filename_, frameIndex_ - 1))
	{
		QMessageBox::warning(this, tr("Error"), tr("Error loading next frame from file %2.").arg(QString::fromStdString(filename_)));
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::toggleAllChannels()
{
	qobject_cast<ViewGraph*>(graphTabWidget->currentWidget())->toggleAllCurves();
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::selectChannels()
{
	bool ok;
	ChannelMask newMask = EditChannelMaskDialog::selectChannels(this, channels, &ok);
	if (ok)
	{
		channels = newMask;
		if (not filename_.empty())
			setFrame(filename_, frameIndex_);
	}
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::setRawFrame(std::ifstream & file, const size_t frameIndex)
{
	LOG_INFO() << "Reading frame no. '" << frameIndex << "' from file...";
	std::auto_ptr<Frame> frame;
	try
	{
		Frame::seekFrame(file, frameIndex);
		frame = Frame::read(file);
	}
	catch (const std::ifstream::failure & e)
	{
		if (file.rdstate() & std::ifstream::eofbit)
		{
			LOG_ERROR() << "Error reading frame: EOF reached.";
		}
		else
		{
			LOG_ERROR() << "Error reading frame: " << e.what();
		}
		return false;
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Error reading frame: " << e.what();
		return false;
	}
	file.close();

	frame->header().print();
	std::cout << std::endl;

	if (frame->header().isBlobFrame())
	{
		numRawFrames_ = 0;
		if (0x7 == frame->header().frameType())
			decodeCoBoTopologyFrame(*frame.get());
		else
			decodeMuTanTFrame(*frame.get());
		return true;
	}
	else if (frame->header().isLayeredFrame())
	{
		LOG_INFO() << "Decoding layered frame...";
		numRawFrames_ = frame->itemCount();
		for (size_t subFrameId=0; subFrameId < numRawFrames_; subFrameId++)
		{
			std::auto_ptr< Frame > subFrame = frame->frameAt(subFrameId);
			if (decodeRawFrame(*subFrame.get()))
			{
				setView(subFrameId);
			}
		}
		return true;
	}
	else
	{
		LOG_INFO() << "Decoding raw CoBo frame...";
		numRawFrames_ = 1;
		if (decodeRawFrame(*frame.get()))
		{
			setView();
			return true;
		}
	}
	return false;
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::setMatrix(std::ifstream & file, const size_t frameIndex)
{
	LOG_INFO() << "Reading matrix no. '" << frameIndex << "' from file...";
	try
	{
		size_t const matrixSize = 4*68*512*sizeof(double);
		file.seekg(matrixSize*frameIndex);
	}
	catch (const std::ifstream::failure & e)
	{
		if (file.rdstate() & std::ifstream::eofbit)
		{
			LOG_ERROR() << "Error reading matrix: EOF reached.";
		}
		else
		{
			LOG_ERROR() << "Error reading matrix: " << e.what();
		}
		return false;
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << "Error reading matrix: " << e.what();
		return false;
	}

	if (decodeMatrix(file))
	{
		numRawFrames_ = 1;
		setView();
		return true;
	}
	return false;
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::setFrame(const std::string & filename, const size_t frameIndex, const size_t & byteOffset)
{
	coboEvent.clear();
	LOG_INFO() << "Opening file '" << filename << "'...";
	std::ifstream file;
	file.exceptions(std::ifstream::failbit | std::ifstream::badbit );
	try
	{
		file.open(filename.c_str(), std::ios::in | std::ios::binary);
	}
	catch (const std::ifstream::failure & e)
	{
		LOG_ERROR() << "Could not open file '" << filename << "': " << strerror(errno);
		return false;
	}

	// Skip bytes
	file.seekg(byteOffset, std::ios_base::cur);

	bool success = false;
	try
	{
#if BOOST_FILESYSTEM_VERSION >= 3
		if (fs::path(filename).extension().string() == ".gmat")
#else
		if (fs::path(filename).extension() == ".gmat")
#endif
			success = setMatrix(file, frameIndex);
		else
			success = setRawFrame(file, frameIndex);

		if (file.is_open())
			file.close();
	}
	catch (std::exception & e)
	{
		LOG_ERROR() << e.what();
		success = false;
	}


	if (success)
	{
		filename_ = filename;
		frameIndex_ = frameIndex;
		setLabels();
		emit frameChanged();
	}
	return success;
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::decodeRawFrame(mfm::Frame & frame)
{
	int status = EXIT_SUCCESS;

	try
	{
		//boost::timer timer;
		if (not	coboEvent.fromFrame(frame))
			return false;
		//coboEvent.fromFrame_hardcoded(frame);
		//LOG_INFO() << "Frame decoded in " << timer.elapsed() << " s";

		// Print meta-data
		LOG_INFO() << "eventTime=" << coboEvent.eventTime();
		LOG_INFO() << "eventIdx=" << coboEvent.eventIdx();
		LOG_INFO() << "coboIdx=" << (short) coboEvent.coboIdx();
		LOG_INFO() << "asadIdx=" << (short) coboEvent.asadIdx();
		LOG_INFO() << "readOffset=" << coboEvent.readOffset();
		LOG_INFO() << "status=" << (short) coboEvent.status();
		// Print hit patterns
		for (short agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			std::ostringstream patternName;
			patternName << "hitPat_" << agetIdx;
			LOG_INFO() << patternName.str() << ": " << coboEvent.hitPattern(agetIdx);
		}
		// Print multiplicity
		for (short agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			LOG_INFO() << "multip_" << agetIdx << '=' << coboEvent.multiplicity(agetIdx);
		}
		// Print sliding window
		LOG_INFO() << "windowOut=" << coboEvent.windowOut();
		// Print last cell index
		for (short agetIdx=0; agetIdx < 4; ++agetIdx)
		{
			LOG_INFO() << "lastCell_" << agetIdx << '=' << coboEvent.lastCell(agetIdx);
		}
	}
	catch (const mfm::Exception &e)
	{
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
	}
	return (EXIT_SUCCESS == status);
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::decodeCoBoTopologyFrame(mfm::Frame & frame)
{
	if (0x7 == frame.header().frameType())
	{
		LOG_INFO() << "Decoding CoBo topology frame...";
		try
		{
			// Print meta-data
			std::ostringstream oss;
			oss << "Topology:";
			// Skip 8 Byte header
			std::istream & stream = frame.serializer().inputStream(8u);
			uint8_t value;
			stream.read(reinterpret_cast<char *>(&value), 1lu);
			oss << " coboIdx=" << (uint16_t) value;
			stream.read(reinterpret_cast<char *>(&value), 1lu);
			oss << " asadMask=" << std::hex << std::showbase << (uint16_t) value << std::dec;
			stream.read(reinterpret_cast<char *>(&value), 1lu);
			oss << " 2pMode=" << std::boolalpha << (bool) value;
			LOG_INFO() << oss.str();
		}
		catch (const mfm::Exception &)
		{}
		return true;
	}
	return false;
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::decodeMuTanTFrame(mfm::Frame & frame)
{
	if (0x8 == frame.header().frameType())
	{
		LOG_INFO() << "Decoding MuTanT data frame...";
		try
		{
			// Print meta-data
			std::ostringstream oss;
			oss << "Event:";
			oss << " TIMESTAMP=" << frame.headerField(8u, 6u).value< uint32_t >();
			oss << " EVENT_NUMBER=" << frame.headerField(14u, 4u).value< uint32_t >();
			oss << " TRIGGER_INFO=" << std::hex << std::showbase << frame.headerField(18u, 2u).value< uint16_t >() << std::dec;
			oss << " MULTIPLICITY_A_MEM="  << frame.headerField(20u, 2u).value< uint16_t >();
			oss << " MULTIPLICITY_B_MEM="  << frame.headerField(22u, 2u).value< uint16_t >();
			oss << " L0_EVT_COUNTER="  << frame.headerField(24u, 4u).value< uint32_t >();
			oss << " L1A_EVT_COUNTER="  << frame.headerField(28u, 4u).value< uint32_t >();
			oss << " L1B_EVT_COUNTER="  << frame.headerField(32u, 4u).value< uint32_t >();
			oss << " L2_EVT_COUNTER="  << frame.headerField(36u, 4u).value< uint32_t >();
			oss << " SCALER1_REG="  << frame.headerField(40u, 4u).value< uint32_t >();
			oss << " SCALER2_REG="  << frame.headerField(44u, 4u).value< uint32_t >();
			oss << " SCALER3_REG="  << frame.headerField(48u, 4u).value< uint32_t >();
			oss << " SCALER4_REG="  << frame.headerField(52u, 4u).value< uint32_t >();
			oss << " SCALER5_REG="  << frame.headerField(56u, 4u).value< uint32_t >();
			oss << " D2P_TIME="  << frame.headerField(60u, 4u).value< uint32_t >();
			LOG_INFO() << oss.str();
		}
		catch (const mfm::Exception &)
		{}
		return true;
	}
	return false;
}
//_________________________________________________________________________________________________
bool CoBoFrameViewer::decodeMatrix(std::ifstream & file)
{
	int status = EXIT_SUCCESS;
	LOG_INFO() << "Decoding matrix...";
	try
	{
		coboEvent.fromMatrix(file);
	}
	catch (const std::exception & e)
	{
		LOG_FATAL() << e.what();
		status = EXIT_FAILURE;
	}
	return (EXIT_SUCCESS == status);
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::setView(const size_t & subFrameId)
{
	if (numRawFrames_*4 != displays.size()) createDisplays(mode3d_);

	const get::CoBoEvent & event = coboEvent;
	const size_t coboIdx = coboEvent.coboIdx();
	const size_t asadIdx = coboEvent.asadIdx();
	try
	{
		for (size_t agetIdx=0; agetIdx < 4; agetIdx++)
		{
			displays[subFrameId*4+agetIdx]->setData(event, agetIdx, channels);
			QString label = (numRawFrames_> 1)?tr("Sub-frame no. %1 - ").arg(subFrameId):QString();
			label += tr("CoBo no. %3 - AsAd no. %1 - AGET no. %2").arg(asadIdx).arg(agetIdx).arg(coboIdx);
			displays[subFrameId*4+agetIdx]->setLabel(label);
		}
	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}
}
//_________________________________________________________________________________________________
void CoBoFrameViewer::setLabels()
{
	QString title = QFileInfo(QString::fromStdString(filename_)).completeBaseName();
	QString eventIdx = QString("Event no. %1").arg(coboEvent.eventIdx());
	for (size_t displayId=0; displayId < 4*numRawFrames_; displayId++)
	{
		QString label = displays[displayId]->getLabel();
		label.prepend(tr("%1 - Frame no. %2 - %3 \n").arg(title).arg(frameIndex_).arg(eventIdx));
		displays[displayId]->setLabel(label);
	}
}
//_________________________________________________________________________________________________
} // namespace get
