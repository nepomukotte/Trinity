/*
 * @file DatabaseFilterProxyModel.h
 * $Id$
 */
#ifndef CCFG_MODEL_DBFILTERPROXYMODEL_H_INCLUDED
#define CCFG_MODEL_DBFILTERPROXYMODEL_H_INCLUDED

#include <QSortFilterProxyModel>
#include <QModelIndex>

namespace CCfg
{
namespace Model
{
/**
 * This is a proxy model to filter database items,
 * keeping either only runs and their conditions or only configurations and their versions.
 *
 * @author sizun
 */
class DatabaseFilterProxyModel : public QSortFilterProxyModel
{
public:
	enum FilterMode { NoFilter=0x3, RunOnly=0x01, ConfigOnly=0x02};
	DatabaseFilterProxyModel(QObject* parent=0, FilterMode mode=NoFilter);
	void setFilterMode(FilterMode newMode);
protected:
	bool filterAcceptsRow(int source_row, const QModelIndex & source_parent) const;
private:
	FilterMode filterMode;
};

}
}
#endif /* CCFG_MODEL_DBFILTERPROXYMODEL_H_INCLUDED */
