/**
 * @file QwtSymbolDialog.cpp
 * @date Oct 15, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: QwtSymbolDialog.cpp 989 2013-10-30 12:11:11Z psizun $
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

#include "QwtSymbolDialog.h"

#include <QComboBox>
#include <QDialogButtonBox>
#include <QBoxLayout>
#include <QLabel>

namespace get
{
namespace gui
{
//_________________________________________________________________________________________________
QStringList QwtSymbolDialog::symbolStyleNames = QStringList()
	        << "NoSymbol"	// No Style. The symbol cannot be drawn.
	        << "Ellipse"	// Ellipse or circle
	        << "Rect"		// Rectangle
	        << "Diamond"	// Diamond
	        << "Triangle"	// Triangle pointing upwards
	        << "DTriangle"	// Triangle pointing downwards
	        << "UTriangle"	// Triangle pointing upwards
	        << "LTriangle"	// Triangle pointing left
	        << "RTriangle"	// Triangle pointing right
	        << "Cross"		// Cross (+)
	        << "XCross"		// Diagonal cross (X)
	        << "HLine"		// Horizontal line
	        << "VLine"		// Vertical line
	        << "Star1"		// X combined with +
	        << "Star2"		// Six-pointed star
	        << "Hexagon";	// Hexagon
//_________________________________________________________________________________________________
QwtSymbolDialog::QwtSymbolDialog(QWidget* parentWidget)
	: QDialog(parentWidget)
{
	QLabel* label = new QLabel("Style of symbol:", this);
	styleComboBox = new QComboBox(this);
	styleComboBox->addItems(symbolStyleNames);

	QBoxLayout* hlayout = new QHBoxLayout;
	hlayout->addWidget(label);
	hlayout->addWidget(styleComboBox);

	QDialogButtonBox* buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok|QDialogButtonBox::Cancel,
			Qt::Horizontal, this);

	QBoxLayout* vlayout = new QVBoxLayout(this);
	vlayout->addLayout(hlayout);
	vlayout->addWidget(buttonBox);

	setLayout(vlayout);
	setWindowTitle(tr("Select symbol"));

	connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
	connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));
}
//_________________________________________________________________________________________________
void QwtSymbolDialog::setSymbolStyle(const QString & style)
{
	styleComboBox->setCurrentIndex(styleComboBox->findText(style));
}
//_________________________________________________________________________________________________
QString QwtSymbolDialog::symbolStyle() const
{
	return styleComboBox->currentText();
}
//_________________________________________________________________________________________________
QString QwtSymbolDialog::getSymbolStyle(const QString & defaultStyle, QWidget* parentWidget)
{
	QwtSymbolDialog dialog(parentWidget);
	dialog.setSymbolStyle(defaultStyle);
	if (QDialog::Accepted == dialog.exec())
	{
		return dialog.symbolStyle();
	}
	return QString();
}
//_________________________________________________________________________________________________
} /* namespace gui */
} /* namespace get */
