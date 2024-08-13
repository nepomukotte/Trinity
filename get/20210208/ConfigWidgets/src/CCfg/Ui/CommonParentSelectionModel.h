#ifndef CCFG_UI_COMMONPARENTSELECTIONMODEL_H_INCLUDED
#define CCFG_UI_COMMONPARENTSELECTIONMODEL_H_INCLUDED

#include <QItemSelectionModel>

namespace CCfg
{
namespace Ui
{

////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a selection model that limits multiselection to items
 * belonging to the same parent.
 *
 * @author fchateau
 */
class CommonParentSelectionModel : public QItemSelectionModel
{
public:
	CommonParentSelectionModel(QAbstractItemModel* model, QObject* parent = 0);
	~CommonParentSelectionModel();

public slots:
	void select(const QModelIndex& index, QItemSelectionModel::SelectionFlags command);
	void select(const QItemSelection& selection, QItemSelectionModel::SelectionFlags command);

private:
	CommonParentSelectionModel();
	CommonParentSelectionModel(const CommonParentSelectionModel& right);
	CommonParentSelectionModel& operator=(const CommonParentSelectionModel& right);
};

}
}

#endif // CCFG_UI_COMMONPARENTSELECTIONMODEL_H_INCLUDED
