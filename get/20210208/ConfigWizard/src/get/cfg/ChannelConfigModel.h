/**
 * @file ChannelConfigModel.h
 * @date Jan 22, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: ChannelConfigModel.h 1105 2014-01-29 08:47:35Z psizun $
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

#ifndef get_cfg_ChannelConfigModel_h_INCLUDED
#define get_cfg_ChannelConfigModel_h_INCLUDED

#include "get/utl/uid.h"
#include "CCfg/Model/ConfigModel.h"
#include <QAbstractProxyModel>
#include <QString>

namespace CCfg
{
namespace Model
{
class ConfigNode;
}
}

namespace get
{
namespace cfg
{
//_________________________________________________________________________________________________
class ChannelConfigModel : public QAbstractProxyModel
{
	Q_OBJECT
public:
	ChannelConfigModel(QObject* parentObject=0);
	virtual ~ChannelConfigModel() {}
	virtual void setSourceModel(QAbstractItemModel *);
	CCfg::Model::ConfigModel* configModel() const;
	virtual QVariant headerData(int section, Qt::Orientation orientation, int role) const;
	virtual QVariant data(const QModelIndex & index, int role = Qt::DisplayRole) const;
	virtual Qt::ItemFlags flags(const QModelIndex & index = QModelIndex()) const;
	virtual int	rowCount(const QModelIndex & parent = QModelIndex()) const;
	virtual int columnCount(const QModelIndex & parent) const;
	virtual QModelIndex index(int, int=0, const QModelIndex & parent = QModelIndex()) const;
	virtual QModelIndex parent(const QModelIndex&) const;
	virtual QModelIndex mapToSource(const QModelIndex&) const;
	virtual QModelIndex mapFromSource(const QModelIndex&) const;
	uint lsbThreshold(const QModelIndex & index) const;
	uint msbThreshold(const QModelIndex & index) const;
	uint threshold(const QModelIndex & index) const;
	QString gain(const QModelIndex & index) const;
	bool isFPN(const QModelIndex & index) const;
	bool isSelected(const QModelIndex & index) const;
	QString triggerInhibition(const QModelIndex & index) const;
	QString isReadIf(const QModelIndex & index) const;
	uint zeroSuppressionThreshold(const QModelIndex & index) const;
	QModelIndex findChannel(const utl::ChannelUID & uid) const;
	CCfg::Model::ConfigNode* channelNode(const QModelIndex & index) const;
	CCfg::Model::ConfigNode* agetNode(const QModelIndex & index) const;
	CCfg::Model::ConfigNode* asadNode(const QModelIndex & index) const;
	CCfg::Model::ConfigNode* coboNode(const QModelIndex & index) const;
private slots:
	void source_modelReset();
private:
	enum Columns
	{
		CoBoColumn, ///< the index of the column displaying the CoBo identifier
		AsAdColumn,	///< the index of the column displaying the AsAd identifier
		AGETColumn, ///< the index of the column displaying the AGET chip identifier
		ChannelColumn,	///< the index of the column displaying the AGET channel number
		MsbThresholdColumn, ///< the index of the column displaying the MSB threshold.
		LsbThresholdColumn, ///< the index of the column displaying the LSB threshold.
		ThresholdColumn, ///< the index of the column displaying the total threshold.
		GainColumn, ///< the index of the column displaying the gain.
		SelectedForTestColumn, ///< the index of the column displaying whether the channel is selected for the test mode.
		TriggerInhibition, ///< the index of the column displaying whether the channel is inhibited.
		ReadIfColumn, ///< the index of the column displaying when the channel is to be read.
		ZeroSupThresholdColumn, ///< the index of the column displaying the zero suppression threshold.
		COLUMN_COUNT	///< the number of column in the model.
	};
	static const std::string columnNames[COLUMN_COUNT];
	static const std::string columnTips[COLUMN_COUNT];
	void buildCache();
	std::string instanceName_;
	QModelIndexList channelIndexes_; ///< Indexes of channels in source model.
	QList< utl::ChannelUID > channelUIDs_;
};
//_________________________________________________________________________________________________
} /* namespace cfg */
} /* namespace get */
#endif /* get_cfg_ChannelConfigModel_h_INCLUDED */
