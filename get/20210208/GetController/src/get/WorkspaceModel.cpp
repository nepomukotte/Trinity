/**
 * $Id: WorkspaceModel.cpp 403 2012-10-05 14:47:10Z psizun $
 * @file WorkspaceModel.cpp
 * @date 30 avr. 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * -----------------------------------------------------------------------------
 */

#include "get/WorkspaceModel.h"
#include "get/TestSummary.h"
#include "get/WorkspaceManager.h"
#include <QFileInfo>

namespace get {
//_____________________________________________________________________________
WorkspaceModel::WorkspaceModel(QObject * parent)
	: QAbstractTableModel(parent)
{
	std::set< std::string > names = WorkspaceManager::instance().testNames();
	std::set< std::string >::const_iterator it = names.begin();
	for (; it != names.end(); ++it)
		testNames_ << QString::fromStdString(*it);
}
//_____________________________________________________________________________
QVariant WorkspaceModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role != Qt::DisplayRole or orientation != Qt::Horizontal)
		return QVariant();
    switch (section)
    {
    	case Name:
    		return tr("Name");
    		break;
    	case Type:
    		return tr("Type");
    		break;
    	case Hardware:
    		return tr("Hardware");
    		break;
    	case Conditions:
    		return tr("Conditions");
    		break;
    	case DataPath:
    		return tr("Data");
    		break;
    	case Ecc:
    		return tr("ECC");
    		break;
    	case Daq:
    		return tr("DAQ");
    		break;
    	case Target:
    		return tr("TARGET");
    		break;
    	default:
    		break;
    }
    return QVariant();
}
//_____________________________________________________________________________
int WorkspaceModel::rowCount(const QModelIndex& /* parent */) const
{
	return testNames_.size();
}
//_____________________________________________________________________________
int WorkspaceModel::columnCount(const QModelIndex& /* parent */) const
{
	return ParamCount;
}
//_____________________________________________________________________________
QVariant WorkspaceModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    std::auto_ptr< get::TestSummary > test = WorkspaceManager::instance().readTest(testNames_.at(index.row()).toStdString());
    switch (index.column())
    {
    	case Name:
    		if (role == Qt::DisplayRole)
    			return QString::fromStdString(test->name());
    		break;
    	case Type:
    		if (role == Qt::DisplayRole)
    			return QString::fromStdString(test->type());
    		break;
    	case Hardware:
    		if (role == Qt::DisplayRole)
    			return QFileInfo(QString::fromStdString(test->hwConfigPath())).fileName();
    		else if (role == Qt::ToolTipRole)
    			return QString::fromStdString(test->hwConfigPath());
    		break;
    	case Conditions:
    		if (role == Qt::DisplayRole)
    			return QFileInfo(QString::fromStdString(test->paramConfigPath())).fileName();
    		else if (role == Qt::ToolTipRole)
    			return QString::fromStdString(test->paramConfigPath());
    		break;
    	case DataPath:
    		if (role == Qt::DisplayRole)
    			return QString::fromStdString(test->dataPath());
    		break;
    	case Ecc:
    		if (role == Qt::DisplayRole)
    			return QString::fromStdString(test->eccEndpoint().toString());
    		break;
    	case Target:
    		if (role == Qt::DisplayRole)
    			return QString::fromStdString(test->targetEndpoint().toString());
    		break;
    	default:
    		break;
    }
    return QVariant();
}
//_____________________________________________________________________________
 /* namespace get */
}
