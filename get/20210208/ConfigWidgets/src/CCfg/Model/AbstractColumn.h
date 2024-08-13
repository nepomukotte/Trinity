#ifndef CCFG_MODEL_ABSTRACTCOLUMN_H_INCLUDED
#define CCFG_MODEL_ABSTRACTCOLUMN_H_INCLUDED

#include <QFont>
#include <QColor>
#include <QString>

namespace CCfg
{
namespace Model
{
// forward declarations
struct AbstractRow;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * An interface that each column of a tree view should implement.
 * It enables the use of a polymorphic approach of column rendering instead
 * of the nested switch-case approach in Qt.
 * The interface takes each method of QAbstractItemModel requesting tree data
 * and makes a distinct method for each relevant role.
 *
 * @author fchateau
 */
struct AbstractColumn
{
	virtual ~AbstractColumn();

	virtual QString getHeaderDisplayName() const = 0;
	virtual QString getHeaderTooltip() const;

	virtual QString getRowDisplayData(const AbstractRow& row) const = 0;
	virtual QString getRowEditionData(const AbstractRow& row) const;
	virtual QString getRowTooltip(const AbstractRow& row) const;
	virtual QFont getRowFont(const AbstractRow& row) const;
	virtual Qt::Alignment getRowAlignment(const AbstractRow& row) const;
	virtual QColor getRowForegroundTextColor(const AbstractRow& row) const;
	virtual QColor getRowBackgroundColor(const AbstractRow& row) const;
	virtual Qt::ItemFlags getRowFlags(const AbstractRow& row) const;

	virtual bool setRowEditionData(const QString& newData, const AbstractRow& row) const;
};

}
}

#endif // CCFG_MODEL_ABSTRACTCOLUMN_H_INCLUDED
