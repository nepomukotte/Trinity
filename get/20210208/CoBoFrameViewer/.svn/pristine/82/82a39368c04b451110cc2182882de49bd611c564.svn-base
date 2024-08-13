/**
 * @file EditChannelMaskDialog.h
 * @date 26 oct. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the CoBoFrameViewer software project.
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

#ifndef get_EditChannelMaskDialog_h_INCLUDED
#define get_EditChannelMaskDialog_h_INCLUDED

#include "get/ChannelMask.h"
#include <QDialog>
class QLineEdit;
class QDialogButtonBox;

namespace get
{
//__________________________________________________________________________________________________
class EditChannelMaskDialog :  public QDialog
{
	Q_OBJECT
public:
	EditChannelMaskDialog(QWidget* parentWidget=0);
	void setChannels(const ChannelMask & mask);
	ChannelMask channels() const;
	static ChannelMask selectChannels(QWidget* widget, const ChannelMask & value = ChannelMask(), bool* ok=0);
	private slots:
	void onMaskChanged();
private:
	QLineEdit* maskEdit;
	QDialogButtonBox* buttonBox;
};
//__________________________________________________________________________________________________
} /* namespace get */
#endif /* get_EditChannelMaskDialog_h_INCLUDED */
