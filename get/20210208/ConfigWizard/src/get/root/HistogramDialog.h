/**
 * @file HistogramDialog.h
 * @date Jan 27, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: HistogramDialog.h 1105 2014-01-29 08:47:35Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the ConfigWizard software project.
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

#ifndef get_root_HistogramDialog_h_INCLUDED
#define get_root_HistogramDialog_h_INCLUDED

#include "get/root/HistogramWidget.h"

#include <QObject>
#ifdef HAVE_QTROOT
#include <QDialog>
class QAbstractButton;
class QCheckBox;
class QDialogButtonBox;
class QMdiArea;
class QPushButton;
#include <TH2.h>
#endif // HAVE_QTROOT

namespace get
{
namespace root
{
//_________________________________________________________________________________________________
#ifndef HAVE_QTROOT
class HistogramDialog : public QObject
{
    Q_OBJECT
};
//_________________________________________________________________________________________________
#else // HAVE_QTROOT
class HistogramDialog : public QDialog
{
    Q_OBJECT
public:
    HistogramDialog(QWidget* parentWidget=0);
	QSize sizeHint () const { return QSize(860, 620); }
    void addHistogram(TH1*, const char* title="", const char* canvasName="canvas");
protected:
    void closeEvent(QCloseEvent *event);
public slots:
	void setMdiModeEnabled(bool mdiModeEnabled);
private slots:
    void on_buttonBox_clicked(QAbstractButton* button);
    void onHistogramClosed();
private:
    void addHistogram(HistogramWidget* mapWidget);
private:
    QMdiArea* mdiArea;
    QDialogButtonBox* buttonBox;
    QPushButton* mdiModeButton;
    size_t histogramCount;
};
#endif // HAVE_QTROOT
//_________________________________________________________________________________________________
} /* namespace root */
} /* namespace get */
#endif /* get_root_HistogramDialog_h_INCLUDED */
