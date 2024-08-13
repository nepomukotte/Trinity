/**
 * @file DataRouterWizard.cpp
 * @date 21 nov. 2012
 * @author sizun
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Patrick Sizun
 *
 * This file is part of the GetController software project.
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

#include "get/daq/DataRouterWizard.h"
#include "get/ControllerOptions.h"
#include "mdaq/DefaultPortNums.h"
#include "utl/Logging.h"

#include <QApplication>
#include <QIcon>
#include <QGroupBox>
#include <QRadioButton>
#include <QLabel>
#include <QLineEdit>
#include <QSpinBox>
#include <QComboBox>
#include <QMessageBox>
#include <QVBoxLayout>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QSpacerItem>
#include <QSettings>
#include <QRegExp>
#include <QRegExpValidator>

namespace get
{
namespace daq
{
//_________________________________________________________________________________________________
DataRouterWizard::DataRouterWizard(QWidget* parentWidget) : QWizard(parentWidget)
{
	configPage = new ConfigPage(this);
	addPage(configPage);
	setPixmap(QWizard::WatermarkPixmap, QApplication::windowIcon().pixmap(QSize(50, 50)));
	setWindowTitle(tr("Data router configuration"));
	setOptions(QWizard::NoBackButtonOnStartPage);

	connect(this, SIGNAL(accepted()), this, SLOT(saveSettings()));
}
//_________________________________________________________________________________________________
void DataRouterWizard::saveSettings()
{
	ControllerOptions & options = ControllerOptions::instance();
	options.setDataRouterExternal(field("external").toBool());
	options.setDataCtrlEndpoint(::utl::net::SocketAddress(field("ctrlAddress").toString().toStdString(), field("ctrlPort").toInt()));
	options.setDataFlowAddress(::utl::net::IpAddress(field("flowAddress").toString().toStdString()));
	options.setDataFlowPort(field("flowPort").toUInt());
	options.setDataFlowType(ConfigPage::flowTypes.at(field("flowType").toInt()).toStdString());
}
//_________________________________________________________________________________________________
const QStringList DataRouterWizard::ConfigPage::flowTypes = QStringList() << "ZBUF" << "TCP" << "ICE" << "FDT" << "Debug";
//_________________________________________________________________________________________________
DataRouterWizard::ConfigPage::ConfigPage(QWidget* parentWidget) : QWizardPage(parentWidget)
{
	setTitle(tr("Setup data router"));
	setSubTitle(tr("Specify the protocol, I.P. address and port of the data router, and whether it is managed by %1 or externally.")
		.arg(QCoreApplication::applicationName()));

	QLabel* label=0;

	QGroupBox* flowGroup = new QGroupBox(tr("Data flow"), this);
	QGridLayout* flowLayout = new QGridLayout(this);
	label = new QLabel(tr("Type:"), this);
	flowLayout->addWidget(label, 0, 0);
	QComboBox* flowTypeBox = new QComboBox(this);
	flowTypeBox->addItems(flowTypes);
	flowLayout->addWidget(flowTypeBox, 0, 2);
	label = new QLabel(tr("IP Address:"), this);
	flowLayout->addWidget(label, 1, 0);
	flowAddressEdit = new QLineEdit(this);
	flowLayout->addWidget(flowAddressEdit, 1, 2);
	label = new QLabel(tr("Port:"), this);
	flowLayout->addWidget(label, 2, 0);
	flowPortBox = new QSpinBox(this);
	flowLayout->addWidget(flowPortBox, 2, 2);
	flowLayout->setColumnStretch(1, 1);
	flowGroup->setLayout(flowLayout);

	ctrlGroup = new QGroupBox(tr("Control"), this);
	QGridLayout* ctrlLayout = new QGridLayout(this);
	label = new QLabel(tr("IP Address:"), this);
	ctrlLayout->addWidget(label, 0, 0);
	QLineEdit* ctrlAddressEdit = new QLineEdit(this);
	ctrlLayout->addWidget(ctrlAddressEdit, 0, 2);
	label = new QLabel(tr("Port:"), this);
	ctrlLayout->addWidget(label, 1, 0);
	QSpinBox* ctrlPortBox = new QSpinBox(this);
	ctrlLayout->addWidget(ctrlPortBox, 1, 2);
	ctrlLayout->setColumnStretch(1, 1);
	ctrlGroup->setLayout(ctrlLayout);

	modeGroup = new QGroupBox(tr("Mode"), this);
	QHBoxLayout* modeLayout = new QHBoxLayout(this);
	label = new QLabel(tr("Router:"), this);
	modeLayout->addWidget(label);
	modeLayout->addStretch(1);
	QRadioButton* internalButton = new QRadioButton(tr("Internal"), this);
	internalButton->setChecked(true);
	QRadioButton* externalButton = new QRadioButton(tr("External"), this);
	modeLayout->addWidget(internalButton);
	modeLayout->addWidget(externalButton);
	modeGroup->setLayout(modeLayout);

	QBoxLayout* layout = new QVBoxLayout(this);
	layout->addWidget(flowGroup);
	layout->addWidget(modeGroup);
	layout->addWidget(ctrlGroup);
	setLayout(layout);

	// Constrain inputs
#if QT_VERSION >= 0x040700
	ctrlAddressEdit->setPlaceholderText("255.255.255.255");
	flowAddressEdit->setPlaceholderText("255.255.255.255");
#endif
	QRegExp ipRegExp("(\\d\\d*)\\.(\\d\\d*)\\.(\\d\\d*)\\.(\\d\\d*)");
	QRegExpValidator* ipValidator = new QRegExpValidator(ipRegExp, this);
	ctrlAddressEdit->setValidator(ipValidator);
	flowAddressEdit->setValidator(ipValidator);
	ctrlPortBox->setMaximum(65535);
	flowPortBox->setMaximum(65535);

	connect(flowTypeBox, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(onFlowTypeChanged(const QString &)));
	connect(flowTypeBox, SIGNAL(currentIndexChanged(const QString &)), this, SIGNAL(completeChanged()));

	// Alignment
	ctrlAddressEdit->setAlignment(Qt::AlignRight);
	flowAddressEdit->setAlignment(Qt::AlignRight);
	ctrlPortBox->setAlignment(Qt::AlignRight);
	flowPortBox->setAlignment(Qt::AlignRight);

	// Register fields
	registerField("external", externalButton);
	registerField("ctrlAddress*", ctrlAddressEdit);
	registerField("ctrlPort*", ctrlPortBox);
	registerField("flowAddress*", flowAddressEdit);
	registerField("flowPort*", flowPortBox);
	registerField("flowType*", flowTypeBox);
}
//_________________________________________________________________________________________________
void DataRouterWizard::ConfigPage::initializePage()
{
	ControllerOptions & options = ControllerOptions::instance();
	setField("external", options.isDataRouterExternal());
	setField("ctrlAddress", QString::fromStdString(options.dataCtrlEndpoint().ipAddress().toString()));
	setField("ctrlPort", options.dataCtrlEndpoint().port());
	setField("flowAddress", QString::fromStdString(options.dataFlowAddress().toString()));
	setField("flowPort", options.dataFlowPort());
	setField("flowType", flowTypes.indexOf(QString::fromStdString(options.dataFlowType())));
}
//_________________________________________________________________________________________________
bool DataRouterWizard::ConfigPage::isComplete () const
{
	return (field("flowType").toInt() >= 0);
}
//_________________________________________________________________________________________________
bool DataRouterWizard::ConfigPage::validatePage()
{
	if (field("flowAddress").toString() == "127.0.0.1")
	{
	        QMessageBox::warning(this, tr("Invalid IP address"), tr("127.0.0.1 is not a valid IP address for the data flow.\nThe data router cannot receive data on the loopback interface."));
	        return false;
	}
	return true;
}
//_________________________________________________________________________________________________
void DataRouterWizard::ConfigPage::onFlowTypeChanged(const QString & flowType)
{
	bool isDebug = (flowType == "Debug");
	flowAddressEdit->setDisabled(isDebug);
	flowPortBox->setDisabled(isDebug);
	ctrlGroup->setDisabled(isDebug);
	modeGroup->setDisabled(isDebug);
}
//_________________________________________________________________________________________________
} /* namespace daq */
} /* namespace get */
