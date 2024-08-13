/*
 * $Id: CoBoFrameViewer.h 1667 2017-07-06 12:58:35Z psizun $
 * @file CoBoFrameViewer.h
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

#ifndef get_CoBoFrameViewer_h_INCLUDED
#define get_CoBoFrameViewer_h_INCLUDED

#include "get/CoBoEvent.h"
#include "get/ChannelMask.h"
#include <QMainWindow>
#include <QSize>
#include <QString>
#include <QWidgetAction>
#include <cstddef>

// Forward declarations
namespace mfm {
	class Frame;
}
class QTimer;
class QAction;
class QTabWidget;
class QToolBar;
class QMenu;
class ViewGraph;

namespace get {
//_________________________________________________________________________________________________
class CoBoDisplay;
//_________________________________________________________________________________________________
class CoBoFrameViewer: public QMainWindow
{
	Q_OBJECT
public:
	CoBoFrameViewer(QWidget* parentWidget=0);
	bool setFrame(const std::string & filename, const size_t frameIndex,  const size_t & byteOffset=0);
	void selectCurrentDisplay(const size_t idx);
	void setChannels(const ChannelMask & channels);
	void setView(const size_t & subFrameId=0);
	void setLabels();
	void startAutoscan(bool enable = true);
	void set3dDisplayMode(bool enable = true);
	QSize sizeHint() const;
signals:
	void frameChanged();
private slots:
	void about();
	void showStatusMessage(const QString &);
	void openFile();
	void selectFrame();
	void toggleAutoscan(bool);
	void toggleDisplay3d(bool);
	void showFullscreen(bool);
	void reloadFrame();
	void selectNextFrame();
	void selectSymbol();
	void selectRandomColorMode();
	void calcFwhm();
	void updateAutoscan();
	void selectPreviousFrame();
	void toggleAllChannels();
	void selectChannels();
	void onFrameChanged();
	void editPeriod();
	void saveAll();
private:
	void createDisplays(bool mode3d=false);
	void createActions();
	bool setRawFrame(std::ifstream & file, const size_t frameIndex);
	bool setMatrix(std::ifstream & file, const size_t frameIndex);
	bool decodeRawFrame(mfm::Frame & frame);
	bool decodeMatrix(std::ifstream & file);
	bool decodeCoBoTopologyFrame(mfm::Frame & frame);
	bool decodeMuTanTFrame(mfm::Frame & frame);
	QTabWidget* graphTabWidget;
	QAction* actionPrevious;
	QAction* actionNext;
	QAction* actionGoTo;
	QAction* actionPlay;
	QAction* actionEditPeriod;
	QAction* action3D;
	QAction* actionSelectChannels;
	QAction* actionCalcFwhm;
	std::map< size_t, CoBoDisplay* > displays; ///< Displays indexed by subFrameId*4+agetIdx
	ChannelMask channels;
	get::CoBoEvent coboEvent;
	std::string filename_;
	size_t frameIndex_;
	QTimer* autoscanTimer;
	int scanPeriod_ms;
	bool mode3d_;
	size_t numRawFrames_;
};
//_________________________________________________________________________________________________
} // namespace get

#endif /* get_CoBoFrameViewer_h_INCLUDED */
