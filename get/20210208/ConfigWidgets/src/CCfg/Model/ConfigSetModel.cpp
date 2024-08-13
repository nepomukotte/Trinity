/**
 * @file ConfigSetModel.cpp
 * @date 29 août 2013
 * @author fchateau
 * @note
 * SVN tag: $Id$
 *
 * Contributors: Frederic Chateau
 *
 * This file is part of the ConfigWidgets software project.
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

#include "ConfigSetModel.h"
#include <utl/Logging.h>
#include <QSettings>
#include <QMap>
#include <QVariant>
#include <QFlags>

namespace CCfg
{
namespace Model
{

// =============================================================================
// CCfg::Model::ConfigSetModel class methods body
// -----------------------------------------------------------------------------

ConfigSetModel::ConfigSetModel(const Params& params)
	: params_(params)
{
	const Params::SlotParams& slotParams = params.configSlots();
	configSlots_.reserve(slotParams.size());
	for(size_t i  = 0; i < slotParams.size(); ++i)
	{
		configSlots_.push_back(new ConfigSlot(slotParams[i], i));
	}
}

ConfigSetModel::~ConfigSetModel()
{
}

void ConfigSetModel::setRepository(const QString& repo)
{
	if(dirRepository_.cd(repo))
	{
		LOG_INFO() << "Repository set to: " << qPrintable(repo);
		repository_ = repo;
		clearConfigItems();
		emit repositoryChanged();
	}
}

bool ConfigSetModel::isDatabaseRepository() const
{
	return false;
}

void ConfigSetModel::setUrl(size_t slotIndex, const QString& url)
{
	assert(slotIndex < configSlots_.size());
	setData(index(0, configSlots_[slotIndex].config()), url, Qt::EditRole);
}

bool ConfigSetModel::setUrl(const QString& slotId, const QString& url)
{
	ConfigSlot* slt = searchById(slotId);
	if(slt != 0)
	{
		setData(index(0, slt->config()), url, Qt::EditRole);
		return true;
	}
	else
	{
		return false;
	}
}

QString ConfigSetModel::relativeUrl(const ConfigItem& item) const
{
	return dirRepository_.relativeFilePath(item.url());
}

ConfigSetModel::ConfigSlot* ConfigSetModel::searchByUrl(const QString& url)
{
	for(size_t i  = 0; i < configSlots_.size(); ++i)
	{
		if(configSlots_[i].config().url() == url)
		{
			return &configSlots_[i];
		}
	}
	return 0;
}

ConfigSetModel::ConfigSlot* ConfigSetModel::searchById(const QString& id)
{
	for(size_t i  = 0; i < configSlots_.size(); ++i)
	{
		if(id == configSlots_[i].id().c_str())
		{
			return &configSlots_[i];
		}
	}
	return 0;
}

template<>
ConfigSetModel::AbstractItem& ConfigSetModel::item(const QModelIndex& index)
{
	return *reinterpret_cast<AbstractItem*>(index.internalPointer());
}

template<class Item>
Item& ConfigSetModel::item(const QModelIndex& index)
{
	return dynamic_cast<Item&>(item<AbstractItem>(index));
}

template<class Item>
Item* ConfigSetModel::indexToItem(const QModelIndex& index)
{
	return dynamic_cast<Item*>(&item<AbstractItem>(index));
}

QModelIndex ConfigSetModel::index(int col, const ConfigItem& item)
{
	return createIndex(0, col, const_cast<ConfigItem*>(&item));
}

QModelIndex ConfigSetModel::index(int col, const ConfigSlot& slot)
{
	return createIndex(slot.index(), col, const_cast<ConfigSlot*>(&slot));
}

QModelIndex ConfigSetModel::index(int row, int column, const QModelIndex& parent) const
{
	if(not parent.isValid())
	{
		if(row >= 0 and static_cast<size_t>(row) < configSlots_.size())
		{
			return createIndex(row, column, const_cast<ConfigSlot*>(&configSlots_[row]));
		}
		else
		{
			return QModelIndex();
		}
	}
	else
	{
		if(row == 0)
		{
			return createIndex(row, column, &item<ConfigSlot>(parent).config());
		}
		else
		{
			return QModelIndex();
		}
	}
}

QModelIndex ConfigSetModel::parent(const QModelIndex& index) const
{
	if(index.isValid())
	{
		ConfigItem* config = indexToItem<ConfigItem>(index);
		if(config == 0)
		{
			return QModelIndex();
		}
		else
		{
			return createIndex(config->slot().index(), 0, &config->slot());
		}
	}
	else
	{
		return QModelIndex();
	}
}

bool ConfigSetModel::hasChildren(const QModelIndex& parent) const
{
	return !parent.isValid() or indexToItem<ConfigSlot>(parent) != 0;
}

int ConfigSetModel::rowCount(const QModelIndex& parent) const
{
	if(parent.isValid())
	{
		return 1;
	}
	else
	{
		return configSlots_.size();
	}
}

int ConfigSetModel::columnCount(const QModelIndex& /*parent*/) const
{
	return 2;
}

QVariant ConfigSetModel::data(const QModelIndex& index, int role) const
{
	if(index.isValid() and index.column() == 0)
	{
		ConfigItem* item = indexToItem<ConfigItem>(index);
		ConfigSlot* slot = (item == 0 ? &ConfigSetModel::item<ConfigSlot>(index) : 0);
		switch(role)
		{
		case Qt::DisplayRole:
		case Qt::EditRole:
		case Qt::ToolTipRole:
			{
				QString name = (item == 0 ? QString::fromStdString(slot->name()) : relativeUrl(*item));
				if(name.isEmpty())
				{
					name = "missing !";
				}
				return QVariant(name);
			}
		case Qt::TextAlignmentRole:
			return QVariant(QFlags<Qt::AlignmentFlag>(Qt::AlignLeft) |= Qt::AlignTop);
		}
	}
	return QVariant();
}

bool ConfigSetModel::setData(const QModelIndex& index, const QVariant& value, int role)
{
	ConfigItem* item = indexToItem<ConfigItem>(index);
	if(role == Qt::EditRole and item != 0)
	{
		const QString& oldUrl = item->url();
		QString url = value.toString();
		item->setUrl(url);
		emit configBindingChanged(QString::fromStdString(item->slot().id()), oldUrl, url);
		return true;
	}
	else
	{
		return false;
	}
}

void ConfigSetModel::onConfigModified(const QString& url)
{
	ConfigSlot* configSlot = searchByUrl(url);
	if(configSlot != 0)
	{
		emit configModified(QString::fromStdString(configSlot->id()));
	}
}

void ConfigSetModel::clearConfigItems()
{
	for(size_t i  = 0; i < configSlots_.size(); ++i)
	{
		configSlots_[i].config().setUrl(QString());
	}
}

template ConfigSetModel::ConfigItem* ConfigSetModel::indexToItem(const QModelIndex& index);
template ConfigSetModel::ConfigSlot* ConfigSetModel::indexToItem(const QModelIndex& index);
template ConfigSetModel::ConfigItem& ConfigSetModel::item(const QModelIndex& index);
template ConfigSetModel::ConfigSlot& ConfigSetModel::item(const QModelIndex& index);

// =============================================================================
// CCfg::Model::ConfigSetModel::Persistency class methods body
// -----------------------------------------------------------------------------

ConfigSetModel::Persistency::Persistency(ConfigSetModel& model)
	: model(model)
{
	LOG_DEBUG() << "Restoring ConfigSetModel...";
	QSettings settings;
	settings.beginGroup("ConfigSetModel");
	model.setRepository(settings.value("repository").toString());
	QMap<QString, QVariant> bindings = settings.value("configBindings").toMap();
    for(QMap<QString, QVariant>::const_iterator binding = bindings.begin(),
    	endBinding = bindings.end(); binding != endBinding; ++binding)
    {
    	model.setUrl(binding.key(), binding.value().toString());
    }
	settings.endGroup();
}

ConfigSetModel::Persistency::~Persistency()
{
	LOG_DEBUG() << "Saving ConfigSetModel...";
	QSettings settings;
	settings.beginGroup("ConfigSetModel");
    settings.setValue("repository", model.repository());
    const ConfigSetModel::ConfigSlotList& cfgSlots = model.configSlots();
    QMap<QString, QVariant> bindings;
    for(size_t i = 0; i < cfgSlots.size(); ++i)
    {
    	bindings.insert(
    		QString::fromStdString(cfgSlots[i].id()),
    		QVariant(cfgSlots[i].config().url()));
    }
    settings.setValue("configBindings", bindings);
	settings.endGroup();
}

// =============================================================================
// CCfg::Model::ConfigSetModel::AbstractItem class methods body
// -----------------------------------------------------------------------------

ConfigSetModel::AbstractItem::~AbstractItem()
{
}

// =============================================================================
// CCfg::Model::ConfigSetModel::Params class methods body
// -----------------------------------------------------------------------------

void ConfigSetModel::Params::appendSlot(const std::string& id, const std::string& name, const std::string& format)
{
	configSlots_.push_back(ConfigSlot::Params(id, name, format));
}

} /* namespace Model */
} /* namespace CCfg */
