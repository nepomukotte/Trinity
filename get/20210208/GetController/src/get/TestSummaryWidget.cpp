/**
 * $Id: TestSummaryWidget.cpp 1399 2015-01-30 15:29:52Z psizun $
 * @file TestSummaryWidget.cpp
 * @date 22 mai 2012
 * @author sizun
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors:
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

#include "get/TestSummaryWidget.h"
#include "get/TestSummary.h"
#include "get/WorkspaceManager.h"
#include "get/GetControllerException.h"
#include "get/ControllerOptions.h"
#include <QPushButton>
#include <QFileInfo>
#include <memory>

namespace get {
//_________________________________________________________________________________________________
TestSummaryWidget::TestSummaryWidget(QWidget* parentWidget)
	: QWidget(parentWidget)
{
	setupUi(this);
	newTestButton->setIcon(QIcon::fromTheme("document-new"));
	selectTestButton->setIcon(QIcon::fromTheme("document-open"));
	editTestButton->setIcon(QIcon::fromTheme("document-properties"));
	editHardwareCConfigButton->setIcon(QIcon(":/application/cconfiged/logo.svg"));
	editConditionsCConfigButton->setIcon(QIcon(":/application/cconfiged/logo.svg"));

	clear();
}
//_________________________________________________________________________________________________
/**
 * Sets actions triggered by configuration file buttons.
 */
void TestSummaryWidget::setConfigActions(QAction* actionEditHardware, QAction* actionEditConditions)
{
	editHardwareCConfigButton->setDefaultAction(actionEditHardware);
	editConditionsCConfigButton->setDefaultAction(actionEditConditions);
}
//_________________________________________________________________________________________________
/**
 * Sets actions triggered by workspace buttons.
 */
void TestSummaryWidget::setWorkspaceActions(QAction* actionNewTest, QAction* actionSelectTest, QAction* actionEditTest)
{
	newTestButton->setDefaultAction(actionNewTest);
	selectTestButton->setDefaultAction(actionSelectTest);
	editTestButton->setDefaultAction(actionEditTest);
}
//_________________________________________________________________________________________________
/**
 * Resets contents.
 */
void TestSummaryWidget::clear()
{
	testNameLabel->setText("-");
	testTypeLabel->setText("-");
	eccEndpointLabel->setText("-");
	targetEndpointLabel->setText("-");
	hardwareConfigPathLabel->setText("");
	editHardwareCConfigButton->setDisabled(true);
	conditionsConfigPathLabel->setText("");
	editConditionsCConfigButton->setDisabled(true);
	dataPathLabel->setText("");
	editTestButton->setDisabled(true);
	setDataRouterInfo();
}
//_________________________________________________________________________________________________
/**
 * Sets test information.
 */
void TestSummaryWidget::setData(::get::TestSummary const & test)
{
	testNameLabel->setText(QString::fromStdString(test.name()));
	testTypeLabel->setText(QString::fromStdString(test.type()));
	eccEndpointLabel->setText(QString::fromStdString(test.eccEndpoint().toString()));
	targetEndpointLabel->setText(QString::fromStdString(test.targetEndpoint().toString()));
	QFileInfo hardwareInfo(QString::fromStdString(test.hwConfigPath()));
	hardwareConfigPathLabel->setText(hardwareInfo.fileName());
	hardwareConfigPathLabel->setToolTip(hardwareInfo.filePath());
	editHardwareCConfigButton->setEnabled(hardwareInfo.exists());
	QFileInfo conditionsInfo(QString::fromStdString(test.paramConfigPath()));
	conditionsConfigPathLabel->setText(conditionsInfo.fileName());
	conditionsConfigPathLabel->setToolTip(conditionsInfo.filePath());
	editConditionsCConfigButton->setEnabled(true);
	dataPathLabel->setText(QString::fromStdString(test.dataPath()));
}
//_________________________________________________________________________________________________
void TestSummaryWidget::setDataRouterInfo()
{
	daqEndpointLabel->setText(QString::fromStdString(ControllerOptions::instance().dataFlowEndpoint().toString()));
}
//_________________________________________________________________________________________________
/**
 * Updates test information.
 */
void TestSummaryWidget::updateData()
{
	try
	{
		std::auto_ptr< ::get::TestSummary > test = WorkspaceManager::instance().currentTest();
		setData(*test.get());
	}
	catch (const GetControllerException & e)
	{
		clear();
	}
}
//_________________________________________________________________________________________________
} /* namespace get */
