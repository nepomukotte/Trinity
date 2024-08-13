#ifndef CCFG_CONTROL_WORKSPACE2_H_INCLUDED
#define CCFG_CONTROL_WORKSPACE2_H_INCLUDED

#include "CCfg/Io/DatabaseStorage.h"
#include "CCfg/Control/Workspace.h"
#include "../Model/ConfigModel2.h"
#include "../Model/DatabaseFilterProxyModel.h"

namespace CCfg
{
namespace Control
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * The workspace controller integrating database features.
 *
 * @author fchateau
 */
class Workspace2 : public Workspace
{
	Q_OBJECT
public:
	// @{ construction
	Workspace2(Ui::ClosableTabWidget* workspaceUi, QWidget* parent = 0);
	virtual ~Workspace2();
	// @}
	// @{ User-level Config Management
	void checkoutConfig();
	void checkout(Model::DatabaseFilterProxyModel::FilterMode mode);
	bool commit();
	bool commitInto();
	// @}
	// @{ User-level Run Management
	void checkoutRun();
	void editRun();
	// @}
	// @{ General Config Management
	bool commit(Model::ConfigModel2* model);
	bool commitInto(Model::ConfigModel2* model);
	// @}
	// @{ Workspace Data Accessors
	Model::ConfigModel2* getCurrentModel();
	// @}
protected:
	virtual Model::ConfigModel* createModel(QObject* parent);
	void notifyCommitSuccess();

	Io::LoginCallbackPtr loginCallback;
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Model::Workspace2 inline methods body
////////////////////////////////////////////////////////////////////////////////

inline Model::ConfigModel2* Workspace2::getCurrentModel()
{
	return static_cast<Model::ConfigModel2*>(Workspace::getCurrentModel());
}

}
}

#endif // CCFG_CONTROL_WORKSPACE2_H_INCLUDED
