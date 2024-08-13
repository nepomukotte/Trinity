/**
 * @file CoBoFrameBrowser.h
 * @date Apr 9, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: CoBoFrameBrowser.h 1435 2015-04-16 08:42:44Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the CoBoFrameViewer software project.
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

#ifndef CoBoFrameBrowser_h_INCLUDED
#define CoBoFrameBrowser_h_INCLUDED

#include "GDataFrame.h"

#include <TGClient.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TGNumberEntry.h>
#include <TH2F.h>
#include <TMultiGraph.h>
#include <TRootEmbeddedCanvas.h>
#include <TTree.h>
#include <RQ_OBJECT.h>
#include <RVersion.h> // ROOT_VERSION_CODE

#include <string>

#if ROOT_VERSION_CODE > 393216
namespace GET {
#else
namespace get {
#endif
//__________________________________________________________________________________________________
class CoBoFrameBrowser
{
	RQ_OBJECT("CoBoFrameBrowser")
public:
	CoBoFrameBrowser(const TGWindow *p, UInt_t w=1280, UInt_t h=1024);
	virtual ~CoBoFrameBrowser();
	virtual void CloseWindow();
	void LoadFile(const std::string & filename);
	// slots
	void DrawPreviousFrame();
	void DrawCurrentFrame();
	void DrawNextFrame();
	void DrawFrame(Long_t frameId);
	void DrawFrame2D(Long_t frameId);
	void DrawFrame3D(Long_t frameId);
	void GoToAction(Long_t);
private:
	void CreateCanvas();
	void DefineButtons();
private:
	TGMainFrame *fMainFrame;
	TCanvas *fCanvas;
    TRootEmbeddedCanvas *fEcanvas;  ///< Browser canvas.
    TGTextButton        *fButPrev;  ///< Previous event button.
    TGTextButton        *fButNext;  ///< Next event button.
    TGNumberEntry       *fButGoTo;  ///< GoTo event button.
    TGCheckButton       *fBut3D;
    TGTextButton        *fButExit;  ///< Exit button.
    TFile* fFile;
    TTree* fTree;
    Long64_t fEntries;
    Long64_t fCurrentId;
    TMultiGraph* multiGraph[4];
    TH2F* histo[4];
	GDataFrame* fFrame;
	bool isUpToDate;
	TString drawOption;
	ClassDef(CoBoFrameBrowser, 1);
};
//__________________________________________________________________________________________________
} // namespace get
#endif /* CoBoFrameBrowser_h_INCLUDED */
