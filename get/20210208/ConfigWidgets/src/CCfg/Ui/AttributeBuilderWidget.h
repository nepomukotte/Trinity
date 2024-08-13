#ifndef CCFG_UI_ATTRIBUTEBUILDERWIDGET_H_INCLUDED
#define CCFG_UI_ATTRIBUTEBUILDERWIDGET_H_INCLUDED

#include "CCfg/TypeInfo.h"
#include <Utils/Ui/DesignedWidget.hpp>
#include <QWidget>

// forward declarations
class Ui_AttributeBuilder;

namespace CCfg
{
// forward declarations
namespace Model { class AttributeBuilderModel; }

namespace Ui
{
///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is a widget where the user can define the properties of new
 * objects and scalars.
 *
 * @author fchateau
 */
class AttributeBuilderWidget : public Utils::Ui::DesignedWidget<Ui_AttributeBuilder, QWidget>
{
	typedef Utils::Ui::DesignedWidget<Ui_AttributeBuilder, QWidget> Base;

	Q_OBJECT

	Model::AttributeBuilderModel* model;	///< data about the new attribute.
	bool validTarget;	///< indicates whether the target node is valid for insertion.

public:
	AttributeBuilderWidget(QWidget* parent = 0);
	virtual ~AttributeBuilderWidget();

	void setModel(Model::AttributeBuilderModel* newModel);
	Model::AttributeBuilderModel* getModel();

signals:
	void insertClicked(Model::AttributeBuilderModel* model, const QString& path, bool insertAsChild);

public slots:
	void onTargetChanged(const QString& path, bool isObject, bool isRoot, bool isInherited);

protected:
	void loadModel();
	void refreshModel();

	void enableDefault(bool isEnabled);
	void enableIndex(bool isEnabled);
	void enableValue(bool isEnabled);

	void commitName();
	void commitIndex();
	void commitType(AttributeTypeId newType);
	void commitValue();

	bool isInsertionReady();
	bool checkData();

	void incrementNumericIndex();

protected slots:
	void onNameChanged(const QString& newName);
	void onIndexChanged(const QString& newIndex);
	void onDefaultToggled(bool checked);
	void onTypeChanged(int newType);
	void onValueChanged(const QString& newValue);
	void onInsertClicked();
};

}
}

#endif // CCFG_UI_ATTRIBUTEBUILDERWIDGET_H_INCLUDED
