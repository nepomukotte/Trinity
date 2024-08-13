/**
 * @file fromArgs.cpp
 * @created Aug 22, 2013
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This file is part of the Mordicus software project.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free
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

#include "fromArgs.h"
#include <QString>
#include <QFileInfo>
#include <QDir>

namespace Utils {
namespace Ui {

const char* fromArgs(int argc, char* argv[], QStringList & cfgFileNames)
{
	QString filePath(argv[0]);
	QFileInfo fileInfo(filePath);
	std::string usage =
		"Usage: " +
		fileInfo.fileName().toStdString() +
		" [--help] [--path <pathName>] [<fileName> ...]";

	// Looking for args
	char* cfgFileName = 0;
	cfgFileNames.clear();

	for (int i=1 ; i < argc ; i++)
	{
		if ((QString(argv[i]) == "--file") and (i+1 < argc))
		{
			cfgFileName = argv[++i];
			cfgFileNames.push_back(argv[i]);
		}
		else if ((QString(argv[i]) == "--path") and (i+1 < argc))
		{
			if (not QDir::setCurrent(QString(argv[++i])))
				throw "Path " + std::string(argv[i]) + " does not exists!";
		}
		else if (QString(argv[i]) == "--help")
		{
			throw usage;
		}
		else if (argv[i][0] != '\0' and argv[i][0] != '-')
		{
			cfgFileName = argv[i];
			cfgFileNames.push_back(argv[i]);
		}
		else
		{
			throw usage;
		}
	}
	if (not cfgFileNames.empty())
	{
		for (int i=0; i < cfgFileNames.size(); ++i)
		{
			filePath = cfgFileNames[i];
			fileInfo.setFile(filePath);
			if (not fileInfo.exists())
				throw "Cannot find file " + fileInfo.fileName().toStdString();
		}
	}
	return cfgFileName;
}

const char* fromArgs(int argc, char* argv[])
{
	QStringList cfgFileNames;
	return fromArgs(argc, argv, cfgFileNames);
}

} // namespace Ui
} // namespace Utils
