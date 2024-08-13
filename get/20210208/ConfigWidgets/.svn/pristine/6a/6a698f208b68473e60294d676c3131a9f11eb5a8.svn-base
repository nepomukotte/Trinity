#include "config.h"
#include "CCfg/Ui/EditorMainWindow.h"
#include "CCfg/Exception.h"
#include "Utils/Ui/fromArgs.h"
#include "Utils/Error/ExcepHandling.h"
#include <QApplication>
#include <QFileInfo>
#include <QDir>
#include <iostream>

using namespace CCfg::Ui;
using namespace Utils::Ui;

int main(int argc, char *argv[])
{
	try {
		// Looking for args
		QStringList cfgFileNames;
		fromArgs(argc, argv, cfgFileNames);

		QApplication app(argc, argv);
		Q_INIT_RESOURCE(ConfigWidgets);
	    QApplication::setOrganizationName("cea");
	    QApplication::setOrganizationDomain("irfu.cea.fr");
	    QApplication::setApplicationName("cconfiged");
	    QApplication::setApplicationVersion(VERSION);
		QApplication::setWindowIcon(QPixmap(":/application/cconfiged/logo.svg"));

		CCfg::Ui::EditorMainWindow editor;
		if (not cfgFileNames.empty())
		{
			editor.openFiles(cfgFileNames);
		}

		editor.show();
		return app.exec();
	}
	catch (const char* e)
	{
		std::cerr << "ERROR: " << e << std::endl;
	}
	catch (const std::string& e)
	{
		std::cerr << "ERROR: " << e << std::endl;
	}
	catch (const CCfg::Exception& e )
	{
		std::cerr << "ERROR: " << e.what() << std::endl;
	}
	catch(...)
	{
		std::cerr << "ERROR: UNKNOWN" << std::endl;
	}

	return 1;
}
