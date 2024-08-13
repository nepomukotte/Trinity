/**
 * @file EditChannelMaskDialog.cpp
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

#include "get/EditChannelMaskDialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QRegExpValidator>
#include <QString>
#include <QPushButton>

namespace get
{
//__________________________________________________________________________________________________
EditChannelMaskDialog::EditChannelMaskDialog(QWidget* parentWidget)
	: QDialog(parentWidget)
{
	QLabel* label = new QLabel(tr("Selected channels:"), this);
	maskEdit = new QLineEdit(this);

	QLayout* hLayout = new QHBoxLayout;
	hLayout->addWidget(label);
	hLayout->addWidget(maskEdit);

	QRegExp rx("^\\d{1,2}(-\\d{1,2})?(,\\d{1,2}(-\\d{1,2})?)*$");
	QValidator *validator = new QRegExpValidator(rx, this);
	maskEdit->setValidator(validator);
	connect(maskEdit, SIGNAL(textChanged(const QString &)), this, SLOT(onMaskChanged()));

	buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok|QDialogButtonBox::Cancel,
			Qt::Horizontal, this);
	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

	QBoxLayout* vLayout = new QVBoxLayout;
	vLayout->addLayout(hLayout);
	vLayout->addWidget(buttonBox);

	setLayout(vLayout);
	setWindowTitle(tr("Select Channels..."));
}
//__________________________________________________________________________________________________
void EditChannelMaskDialog::setChannels(const ChannelMask & mask)
{
	maskEdit->setText(QString::fromStdString(mask.toString()));
}
//__________________________________________________________________________________________________
ChannelMask EditChannelMaskDialog::channels() const
{
	return ChannelMask::fromString(maskEdit->text().toStdString());
}
//__________________________________________________________________________________________________
void EditChannelMaskDialog::onMaskChanged()
{
	buttonBox->button(QDialogButtonBox::Ok)->setEnabled(maskEdit->hasAcceptableInput());
}
//__________________________________________________________________________________________________
ChannelMask EditChannelMaskDialog::selectChannels(QWidget* widget, const ChannelMask & value, bool* ok)
{
	EditChannelMaskDialog dialog(widget);
	dialog.setChannels(value);

	if (dialog.exec())
	{
		if (ok) *ok=true;
		return dialog.channels();
	}
	else
	{
		if (ok) *ok=false;
	}
	return ChannelMask();
}
//__________________________________________________________________________________________________
} /* namespace get */
