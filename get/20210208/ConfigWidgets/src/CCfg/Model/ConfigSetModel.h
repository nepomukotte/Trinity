/**
 * @file ConfigSetModel.h
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

#ifndef CCfg_Model_ConfigSetModel_h_INCLUDED
#define CCfg_Model_ConfigSetModel_h_INCLUDED

#include <QAbstractItemModel>
#include <QDir>
#ifndef Q_MOC_RUN   // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/ptr_container/ptr_vector.hpp>
#endif
#include <vector>

namespace CCfg
{
namespace Model
{

class ConfigSetModel : public QAbstractItemModel
{
	Q_OBJECT
public:
	class Params;
	class Persistency;
	class AbstractItem;
	class ConfigSlot;
	class ConfigItem;
	typedef boost::ptr_vector<ConfigSlot> ConfigSlotList;

	ConfigSetModel(const Params& params);
	virtual ~ConfigSetModel();

	//@{ accessors
	const ConfigSlotList& configSlots() const { return configSlots_; }
	const QString& repository() const { return repository_; }
	void setRepository(const QString& repo);
	bool isDatabaseRepository() const;
	void setUrl(size_t slotIndex, const QString& url);
	bool setUrl(const QString& slotId, const QString& url);
	QString relativeUrl(const ConfigItem& item) const;
	ConfigSlot* searchByUrl(const QString& url);
	ConfigSlot* searchById(const QString& id);
	//@}

	//@{ Index <-> object conversion
	template<class Item> static Item& item(const QModelIndex& index);
	template<class Item> static Item* indexToItem(const QModelIndex& index);
	QModelIndex index(int col, const ConfigItem& item);
	QModelIndex index(int col, const ConfigSlot& slot);
	//@}

	//@{ QAbstractItemModel overriden methods
	QModelIndex index(int row, int column, const QModelIndex& parent) const;
	QModelIndex parent(const QModelIndex& index) const;

	bool hasChildren(const QModelIndex& parent) const;
	int rowCount(const QModelIndex& parent) const;
	int columnCount(const QModelIndex& parent) const;

//	Qt::ItemFlags flags(const QModelIndex& index) const;
//	QVariant headerData(int section, Qt::Orientation orientation, int role) const;
	QVariant data(const QModelIndex& index, int role) const;

	bool setData(const QModelIndex& index, const QVariant& value, int role);
	//@}

signals:
	void repositoryChanged();
	void configBindingChanged(const QString& configSlotName, const QString& oldUrl, const QString& newUrl);
	void configModified(const QString& configSlotName);

public slots:
	void onConfigModified(const QString& configSlotName);

protected:
	void clearConfigItems();

private:
	const Params& params_;
	ConfigSlotList configSlots_;
	QString repository_;
	QDir dirRepository_;
};

class ConfigSetModel::Persistency
{
public:
	explicit Persistency(ConfigSetModel& model);
	~Persistency();

private:
	ConfigSetModel& model;
};

class ConfigSetModel::AbstractItem
{
public:
	virtual ~AbstractItem() = 0;
};

class ConfigSetModel::ConfigItem : public AbstractItem
{
public:
	ConfigItem(ConfigSlot& slot) : slot_(slot) {}

	ConfigSlot& slot() { return slot_; }
	const ConfigSlot& slot() const { return slot_; }

	const QString& url() const { return url_; }

private:
	friend class ConfigSetModel;
	void setUrl(const QString& url) { url_ = url; } // use ConfigSetModel::setUrl, otherwise no signal emitted

	ConfigSlot& slot_;
	QString url_;
};

class ConfigSetModel::ConfigSlot : public AbstractItem
{
public:
	struct Params
	{
		Params(const std::string& i = std::string(),
		       const std::string& n = std::string(),
		       const std::string& f = std::string())
			: id(i), name(n), format(f)
		{
		}
		std::string id;
		std::string name;
		std::string format;
	};

	ConfigSlot(const Params& params, int index)	: params_(params), config_(*this), index_(index) {}

	const ConfigItem& config() const { return config_; }
	ConfigItem& config() { return config_; }

	const std::string& id()  const { return params_.id; }
	const std::string& name() const { return params_.name; }
	const std::string& format()  const { return params_.format; }
	int index() const { return index_; }

private:
	const Params& params_;
	ConfigItem config_;
	int index_;
};

class ConfigSetModel::Params
{
public:
	class ConfigEntry;

	typedef std::vector<ConfigSlot::Params> SlotParams;

	void appendSlot(const std::string& id, const std::string& name, const std::string& format);

	const SlotParams& configSlots() const { return configSlots_; }
	SlotParams& configs() { return configSlots_; }

private:
	SlotParams configSlots_;
};

} /* namespace Model */
} /* namespace CCfg */

#endif /* CCfg_Model_ConfigSetModel_h_INCLUDED */
