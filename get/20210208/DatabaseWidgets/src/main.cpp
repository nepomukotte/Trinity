#include "config.h"
#include "CCfg/Ui/EditorMainWindow2.h"
#include "Utils/Ui/fromArgs.h"
#include "Utils/Error/ExcepHandling.h"
#include <QApplication>
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
		Q_INIT_RESOURCE(DatabaseWidgets);
	    QApplication::setOrganizationName("cea");
	    QApplication::setOrganizationDomain("irfu.cea.fr");
	    QApplication::setApplicationName("dbconfiged");
	    QApplication::setApplicationVersion(VERSION);
		QApplication::setWindowIcon(QPixmap(":/application/dbconfiged/logo.svg"));

		CCfg::Ui::EditorMainWindow2 editor;
		if(not cfgFileNames.empty())
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
