#ifndef CCFG_CONTROL_BASECONFIGEDITORCONTROLLER_H_INCLUDED
#define CCFG_CONTROL_BASECONFIGEDITORCONTROLLER_H_INCLUDED

#include <QObject>

// forward declarations
class QString;


namespace CCfg
{
// forward declarations
namespace Model {
	class ConfigModel;
	class AttributeBuilderModel;
}
namespace Ui { class ConfigEditorTree; }

namespace Control
{
// forward declarations
class AbstractDelegate;

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This is the base class for configuration editor controllers.
 * It provides shared attributes and methods that would be needed to create
 * other controllers.
 *
 * @author fchateau
 */
class BaseConfigEditorController : public QObject
{
	Q_OBJECT
protected:
	Model::ConfigModel* model;	///< the model
	Ui::ConfigEditorTree* view;	///< the view

	AbstractDelegate* delegate;///< the edition delegate.
public:
	BaseConfigEditorController(QObject* parent = 0);
	virtual ~BaseConfigEditorController();

	virtual void setModel(Model::ConfigModel* newModel);
	virtual void setView(Ui::ConfigEditorTree* newView);
	Model::ConfigModel* getModel();
	Ui::ConfigEditorTree* getView();

	void setDelegate(AbstractDelegate* delegate);

public slots:
	virtual void onInsertNew(Model::AttributeBuilderModel* builder, const QString& path, bool insertAsChild);

private:
	BaseConfigEditorController(const BaseConfigEditorController&);
	BaseConfigEditorController& operator=(const BaseConfigEditorController&);
};

}
}

#endif // CCFG_CONTROL_BASECONFIGEDITORCONTROLLER_H_INCLUDED
