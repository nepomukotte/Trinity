#include "BaseConfigEditorController.h"
#include "Delegate/ConfigEditorDelegate.h"
#include "../Ui/ConfigEditorTree.h"
#include "../Model/ConfigModel.h"

namespace CCfg
{
namespace Control
{

///////////////////////////////////////////////////////////////////////////////
// CCfg::Control::BaseConfigEditorController class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
BaseConfigEditorController::BaseConfigEditorController(QObject* parent)
	: QObject(parent), model(0), view(0), delegate(0)
{
}

/**
 * Destructor.
 */
BaseConfigEditorController::~BaseConfigEditorController()
{
}

/**
 * Changes the data model of the configuration editor.
 * @param newModel the new model for this controller.
 */
void BaseConfigEditorController::setModel(Model::ConfigModel* newModel)
{
	if(newModel != model)
	{
		model = newModel;
		if(view != 0)
		{
			view->setModel(model);
		}
	}
}

/**
 * Setter of the "view" property.
 * Changes the view for which user input is processed by this controller.
 * @param newView the new view for this controller.
 */
void BaseConfigEditorController::setView(Ui::ConfigEditorTree* newView)
{
	if(newView != view)
	{
		view = newView;
		if(view != 0)
		{
			view->setModel(model);
			setDelegate(delegate);
		}
	}
}

/**
 * Returns the data model of the configuration editor.
 */
Model::ConfigModel* BaseConfigEditorController::getModel()
{
	return model;
}

/**
 * Returns the view (widget) of the configuration editor.
 */
Ui::ConfigEditorTree* BaseConfigEditorController::getView()
{
	return view;
}

/**
 * Setter of the "delegate" property.
 * Changes the Qt delegate: the object responsible for direct data edition by user.
 */
void BaseConfigEditorController::setDelegate(AbstractDelegate* newDelegate)
{
	delegate = newDelegate;

	if(view != 0) // must update the view delegate
	{
		if(delegate != 0)
		{
			QAbstractItemDelegate* oldDelegate = view->itemDelegate();
			if(delegate != oldDelegate)
			{
				view->setItemDelegate(delegate);
				view->setEditTriggers(QAbstractItemView::DoubleClicked | QAbstractItemView::EditKeyPressed);

				delete oldDelegate;
			}
		}
		else
		{
			view->setEditTriggers(QAbstractItemView::NoEditTriggers);
		}
	}
}

/**
 * Reacts to an insertion order.
 * @param builder the model which contains all informations about the attribute
 * to create.
 * @param path the path where the new attribute should be created.
 * @param insertAsChild if true, the attribute should be inserted as the first
 * child of the node whose path is specified, otherwise, it should be inserted
 * as a sibling just after it.
 */
void BaseConfigEditorController::onInsertNew(Model::AttributeBuilderModel* /*builder*/, const QString& /*path*/, bool /*insertAsChild*/)
{
}

}
}
