/**
 * @file get-config-wizard.h
 * @date Nov 19, 2013
 * @author sizun
 *
 * @note SVN tag: $Id: get-config-wizard.cpp 1103 2014-01-27 16:29:15Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note
 * @note This file is part of the GetConfigWizard software project.
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
#include "get/ConfigWizardMainWindow.h"
#include "utl/Logging.h"

#include <QApplication>
#include <cstdlib>

int main(int argc, char* argv[])
{
	int status = EXIT_SUCCESS;

	try
	{
		QApplication app(argc, argv);
		//Q_INIT_RESOURCE(CoBoFrameViewer);
		//Q_INIT_RESOURCE(GetConfigWizard);
		//QApplication::setWindowIcon(QPixmap(":/get/logo.png"));
		QApplication::setWindowIcon(QPixmap(":/application/cconfiged/logo.svg"));
		QApplication::setOrganizationName("actar-get");
		QApplication::setOrganizationDomain("www-actar-get.cea.fr");
		QApplication::setApplicationName(QObject::tr("get-config-wizard"));

		get::ConfigWizardMainWindow wizard;
		wizard.show();

		status = app.exec();
	}
	catch (const std::exception & e)
	{
		LOG_FATAL() << e.what();
		return EXIT_FAILURE;
	}
	return status;
}
