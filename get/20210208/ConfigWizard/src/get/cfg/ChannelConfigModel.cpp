/**
 * @file ChannelConfigModel.cpp
 * @date Jan 22, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: ChannelConfigModel.cpp 1105 2014-01-29 08:47:35Z psizun $
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

#include "ChannelConfigModel.h"
#include "CCfg/View/ClassIterator.h"
#include "utl/Logging.h"
#include <QStringList>

using ::get::utl::ChannelUID;
using CCfg::Model::ConfigNode;
using namespace CCfg;

namespace get
{
namespace cfg
{
const std::string ChannelConfigModel::columnNames[COLUMN_COUNT] =
{
		"CoBo", "AsAd", "AGET", "Channel", "MSB", "LSB", "Threshold", "Gain", "Selected",
		"Inhib.", "Read", "Zero"
};
const std::string ChannelConfigModel::columnTips[COLUMN_COUNT] =
{
		"", "", "", "",
		"3 most significant bits of DAC threshold", "4 least significant bits of DAC threshold", "Total channel DAC threshold (7 bits)",
		"Gain", "Whether the channel is selected for the test mode.",
		"Whether to inhibit channel contribution to hit register or multiplicity", "When the channel is to be read", "Zero suppression threshold"
};
//_________________________________________________________________________________________________
ChannelConfigModel::ChannelConfigModel(QObject* parentObject) : QAbstractProxyModel(parentObject)
{
}
//_________________________________________________________________________________________________
CCfg::Model::ConfigModel* ChannelConfigModel::configModel() const
{
	return dynamic_cast<CCfg::Model::ConfigModel*>(QAbstractProxyModel::sourceModel());
}
//_________________________________________________________________________________________________
void ChannelConfigModel::setSourceModel(QAbstractItemModel* model)
{
    if (sourceModel())
    {
        //disconnect(sourceModel(), SIGNAL(rowsAboutToBeInserted(QModelIndex, int, int)), this, SLOT(source_modelReset()));
        disconnect(sourceModel(), SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(source_modelReset()));

        //disconnect(sourceModel(), SIGNAL(rowsAboutToBeRemoved(QModelIndex, int, int)), this, SLOT(source_modelReset()));
        disconnect(sourceModel(), SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(source_modelReset()));

        disconnect(sourceModel(), SIGNAL(dataChanged(QModelIndex,QModelIndex)), this, SLOT(source_modelReset()));
        disconnect(sourceModel(), SIGNAL(modelReset()), this, SLOT(source_modelReset()));

        disconnect(sourceModel(), SIGNAL(layoutAboutToBeChanged()), this, SIGNAL(layoutAboutToBeChanged()));
        disconnect(sourceModel(), SIGNAL(layoutChanged()), this, SIGNAL(layoutChanged()));
    }

	QAbstractProxyModel::setSourceModel(model);
	buildCache();

    if (sourceModel())
    {
        //connect(sourceModel(), SIGNAL(rowsAboutToBeInserted(QModelIndex, int, int)), this, SLOT(source_modelReset()));
        connect(sourceModel(), SIGNAL(rowsInserted(QModelIndex, int, int)), this, SLOT(source_modelReset()));

        //connect(sourceModel(), SIGNAL(rowsAboutToBeRemoved(QModelIndex, int, int)), this, SLOT(source_modelReset()));
        connect(sourceModel(), SIGNAL(rowsRemoved(QModelIndex, int, int)), this, SLOT(source_modelReset()));

        connect(sourceModel(), SIGNAL(dataChanged(QModelIndex,QModelIndex)), this, SLOT(source_modelReset()));
        connect(sourceModel(), SIGNAL(modelReset()), this, SLOT(source_modelReset()));

        connect(sourceModel(), SIGNAL(layoutAboutToBeChanged()), this, SIGNAL(layoutAboutToBeChanged()));
        connect(sourceModel(), SIGNAL(layoutChanged()), this, SIGNAL(layoutChanged()));
    }

	reset();
}
//_________________________________________________________________________________________________
void ChannelConfigModel::source_modelReset()
{
    buildCache();
    reset();
}
//_________________________________________________________________________________________________
QVariant ChannelConfigModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	QVariant data;
	if(orientation == Qt::Horizontal)
	{
		switch (role)
		{
		case Qt::DisplayRole:
			data = QString::fromStdString(columnNames[section]);
			break;
		case Qt::ToolTipRole:
			data = QString::fromStdString(columnTips[section]);
			break;
		default:
			data = QAbstractItemModel::headerData(section, orientation, role);
			break;
		}
	}
	else
	{
		switch (role)
		{
		case Qt::DisplayRole:
			data = QVariant(section+1);
			break;
		default:
			data = QAbstractItemModel::headerData(section, orientation, role);
			break;
		}
	}
	return data;
}
//_________________________________________________________________________________________________
Qt::ItemFlags ChannelConfigModel::flags(const QModelIndex & index) const
{
	Qt::ItemFlags flags = Qt::NoItemFlags;
	if (index.isValid())
	{
		flags = Qt::ItemIsEnabled | Qt::ItemIsSelectable;
		if (LsbThresholdColumn == index.column())
			flags |= Qt::ItemIsEditable;
	}
	return flags;
}
//_________________________________________________________________________________________________
QVariant ChannelConfigModel::data(const QModelIndex & index, int role) const
{
	if (not index.isValid())
		return QVariant();
	if (role == Qt::DisplayRole)
	{
		switch(index.column())
		{
		case ChannelColumn:
			return QString::fromStdString(channelNode(index)->getIndex());
			break;
		case AGETColumn:
			return QString::fromStdString(agetNode(index)->getIndex());
			break;
		case AsAdColumn:
			return QString::fromStdString(asadNode(index)->getIndex());
			break;
		case CoBoColumn:
			return QString::fromStdString(coboNode(index)->getIndex());
			break;
		case MsbThresholdColumn:
			if (not  isFPN(index))
				return msbThreshold(index);
			break;
		case LsbThresholdColumn:
			if (not  isFPN(index))
				return lsbThreshold(index);
			break;
		case ThresholdColumn:
			if (not  isFPN(index))
				return threshold(index);
			break;
		case GainColumn:
			if (not  isFPN(index))
				return gain(index);
			break;
		case SelectedForTestColumn:
			return isSelected(index);
			break;
		case TriggerInhibition:
			return triggerInhibition(index);
			break;
		case ReadIfColumn:
			return isReadIf(index);
			break;
		case ZeroSupThresholdColumn:
			return zeroSuppressionThreshold(index);
			break;
		default:
			return QVariant();
			break;
		}
	}
	else if (Qt::TextAlignmentRole == role)
	{
		return QVariant(Qt::AlignRight | Qt::AlignVCenter);
	}
	else if (Qt::ToolTipRole == role)
	{
		return QString::fromStdString(columnTips[index.column()]);
	}
	else if (Qt::ForegroundRole == role)
	{
		if ((isFPN(index) and ChannelColumn == index.column()))
			return Qt::red;
	}
	return QVariant();
}
//_________________________________________________________________________________________________
QModelIndex ChannelConfigModel::findChannel(const utl::ChannelUID & uid) const
{
	QModelIndex idx;
	int row = channelUIDs_.indexOf(uid);
	if (row >= 0)
		idx = index(row);
	return idx;
}
//_________________________________________________________________________________________________
ConfigNode* ChannelConfigModel::channelNode(const QModelIndex & index) const
{
	if (not index.isValid())
		return 0;
	return configModel()->indexToNode(channelIndexes_.at(index.row()));
}
//_________________________________________________________________________________________________
ConfigNode* ChannelConfigModel::agetNode(const QModelIndex & index) const
{
	if (not index.isValid())
		return 0;
	return channelNode(index)->getParent();
}
//_________________________________________________________________________________________________
ConfigNode* ChannelConfigModel::asadNode(const QModelIndex & index) const
{
	if (not index.isValid())
		return 0;
	return 	agetNode(index)->getParent();
}
//_________________________________________________________________________________________________
ConfigNode* ChannelConfigModel::coboNode(const QModelIndex & index) const
{
	if (not index.isValid())
		return 0;
	return 	asadNode(index)->getParent();
}
//_________________________________________________________________________________________________
/**
 * Returns value of the 4 least significant bits of the channel threshold.
 */
uint ChannelConfigModel::lsbThreshold(const QModelIndex & index) const
{
	ConfigNode* channel = channelNode(index);
	return 	(channel?channel->search("LSBThresholdValue", "")->getValue().getInt():0);
}
//_________________________________________________________________________________________________
/**
 * Returns value of the 3 most significant bits of the channel threshold (common to all the AGET chip).
 */
uint ChannelConfigModel::msbThreshold(const QModelIndex & index) const
{
	ConfigNode* aget = agetNode(index);
	return 	(aget?aget->search("Global", "")->search("Reg1", "")->search("GlobalThresholdValue", "")->getValue().getInt():0);
}
//_________________________________________________________________________________________________
/**
 * Returns value of the total 7-bit channel threshold.
 */
uint ChannelConfigModel::threshold(const QModelIndex & index) const
{
	return 	((msbThreshold(index)) << 4 | lsbThreshold(index));
}
//_________________________________________________________________________________________________
/**
 * Returns value of the channel gain.
 */
QString ChannelConfigModel::gain(const QModelIndex & index) const
{
	ConfigNode* channel = channelNode(index);
	return 	QString::fromStdString((channel?channel->search("Gain", "")->getValue().getString():""));
}
//_________________________________________________________________________________________________
/**
 * Returns true if channel is a FPN channel.
 */
bool ChannelConfigModel::isFPN(const QModelIndex & index) const
{
	static const QStringList fpnChannels = QStringList() << "11" << "22" << "45" << "56";
	ConfigNode* channel = channelNode(index);
	return 	fpnChannels.contains(QString::fromStdString(channel->getIndex()));
}
//_________________________________________________________________________________________________
/**
 * Returns true if channel is selected for the test mode (AGET registers 3 and 4).
 */
bool ChannelConfigModel::isSelected(const QModelIndex & index) const
{
	ConfigNode* channel = channelNode(index);
	return 	(channel?channel->search("isSelectedforTestMode", "")->getValue().getBool():false);
}
//_________________________________________________________________________________________________
/**
 * Returns trigger inhibition mode for channel.
 */
QString ChannelConfigModel::triggerInhibition(const QModelIndex & index) const
{
	ConfigNode* channel = channelNode(index);
	return 	QString::fromStdString(channel?channel->search("TriggerInhibition", "")->getValue().getString():"none");
}
//_________________________________________________________________________________________________
/**
 * Returns condition for the channel to be read.
 */
QString ChannelConfigModel::isReadIf(const QModelIndex & index) const
{
	ConfigNode* channel = channelNode(index);
	return 	QString::fromStdString(channel?channel->search("Reading", "")->getValue().getString():"always");
}
//_________________________________________________________________________________________________
/**
 * Returns zero suppression threshold.
 */
uint ChannelConfigModel::zeroSuppressionThreshold(const QModelIndex & index) const
{
	ConfigNode* channel = channelNode(index);
	return 	(channel?channel->search("zeroSuppressionThreshold", "")->getValue().getInt():0);
}
//_________________________________________________________________________________________________
/**
 * Returns number of AGET channels.
 */
int ChannelConfigModel::rowCount(const QModelIndex & parent) const
{
    if (parent.isValid()) return 0;
    if (!sourceModel()) return 0;
	return channelIndexes_.size();
}
//_________________________________________________________________________________________________
int ChannelConfigModel::columnCount(const QModelIndex & parent) const
{
    if (parent.isValid()) return 0;

	return COLUMN_COUNT;
}
//_________________________________________________________________________________________________
QModelIndex ChannelConfigModel::index(int row, int column, const QModelIndex &) const
{
	return createIndex(row, column);
}
//_________________________________________________________________________________________________
QModelIndex ChannelConfigModel::parent(const QModelIndex &) const
{
	return QModelIndex();
}
//_________________________________________________________________________________________________
QModelIndex ChannelConfigModel::mapToSource(const QModelIndex & proxy) const
{
    if (!sourceModel())
        return QModelIndex();

    if (not proxy.isValid())
    	return QModelIndex();

    return configModel()->nodeToIndex(channelNode(proxy));
}
//_________________________________________________________________________________________________
QModelIndex ChannelConfigModel::mapFromSource(const QModelIndex & source) const
{
    if (!sourceModel())
        return QModelIndex();

    if (!source.parent().isValid())
        return QModelIndex();

    if (not channelIndexes_.contains(source))
    	return QModelIndex();

    return QModelIndex();
    //return index(source.row(), source.column());
}
//_________________________________________________________________________________________________
void ChannelConfigModel::buildCache()
{
	if (not configModel()) return;

	channelIndexes_.clear();
	channelUIDs_.clear();

	// Distinguish between GetController configurations and ECC configurations.
	ConfigNode* rootNode = configModel()->indexToNode(configModel()->index(0, 0, QModelIndex()));
	ConfigNode* topNode = rootNode;
	instanceName_ = "CoBo";
	if (rootNode->search("Node", "CoBo") != rootNode->end())
	{
		instanceName_ = "Instance";
		topNode = &(*rootNode->search("Node", "CoBo"));
	}

	// Loop over CoBo boards, AsAd boards, AGET chips, AGET channels
	for (ConfigNode::const_iterator cobo = topNode->begin(); cobo != topNode->end(); ++cobo)
	{
		if (cobo->getName() != instanceName_ or cobo->getIndex() == DEFAULT_OBJ_INDEX) continue;
		for (ConfigNode::const_iterator asad = cobo->begin(); asad != cobo->end(); ++asad)
		{
			if (asad->getName() != "AsAd" or asad->getIndex() == DEFAULT_OBJ_INDEX) continue;
			for (ConfigNode::const_iterator aget = asad->begin(); aget != asad->end(); ++aget)
			{
				if (aget->getName() != "Aget" or aget->getIndex() == DEFAULT_OBJ_INDEX) continue;
				for (ConfigNode::const_iterator channel = aget->begin(); channel != aget->end(); ++channel)
				{
					if (channel->getName() != "channel" or channel->getIndex() == DEFAULT_OBJ_INDEX) continue;
					//LOG_DEBUG() << channel->getPath();
					channelIndexes_ << configModel()->nodeToIndex(const_cast<ConfigNode*>(&(*channel)), 0);
					channelUIDs_ << ChannelUID(cobo->getIndex(), asad->getIndex(), aget->getIndex(), channel->getIndex());
				}
			}
		}
	}
	//LOG_DEBUG() << "Found " << channelIndexes_.size() << " channels(s)";
}
//_________________________________________________________________________________________________
/*virtual bool ChannelConfigModel::setData(const QModelIndex & index, const QVariant & value, int role)
{
	if (!mapToSource(i).model()) return false;
	QAbstractItemModel* m = const_cast<QAbstractItemModel *>(mapToSource(i).model());
	return m->setData(mapToSource(i),v,r);
}*/
//_________________________________________________________________________________________________
} /* namespace cfg */
} /* namespace get */
