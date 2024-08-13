/**
 * @file AttributeFinderDialog.cpp
 * @date 30 nov. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the ConfigWidgets software project.
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

#include "CCfg/Ui/AttributeFinderDialog.h"
#include "CCfg/Ui/AttributeFinder.ui.h"
#include <QLineEdit>
#include <QStringList>
#include <QPushButton>
#include <QCheckBox>

namespace CCfg
{
namespace Ui
{
//_________________________________________________________________________________________________
AttributeFinderDialog::AttributeFinderDialog(QWidget* parentWidget)
	: Utils::Ui::DesignedWidget<Ui_AttributeFinder, QDialog>(parentWidget)
{
	setupUi();
}
//_________________________________________________________________________________________________
AttributeFinderDialog::~AttributeFinderDialog()
{
}
//_________________________________________________________________________________________________
void AttributeFinderDialog::setupUi()
{
	QPushButton* findButton = ui_->buttonBox->addButton(tr("Find"), QDialogButtonBox::AcceptRole);
	findButton->setDefault(true);

	ui_->columnBox->addItems(QStringList() << tr("Name") << tr("Value") << tr("Unit") << tr("Range"));
}
//_________________________________________________________________________________________________
void AttributeFinderDialog::accept()
{
	findAttribute();
	if (not ui_->keepOpenBox->isChecked())
		hide();
}
//_________________________________________________________________________________________________
void AttributeFinderDialog::reject()
{
	hide();
}
//_________________________________________________________________________________________________
void AttributeFinderDialog::findAttribute()
{
	emit searchAttribute(ui_->lineEdit->text(),
			ui_->caseSensitiveCheckBox->isChecked(),
			ui_->exactMatchCheckBox->isChecked(),
			ui_->columnBox->currentIndex());
}
//_________________________________________________________________________________________________
} /* namespace Ui */
} /* namespace CCfg */
