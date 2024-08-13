#include "AttributeBuilderModel.h"
#include "CCfg/Scalar.h"
#include "CCfg/Object.h"
#include "CCfg/Attribute.h"

namespace CCfg
{
namespace Model
{
///////////////////////////////////////////////////////////////////////////////
// CCfg::Model::AttributeBuilderModel class methods body
///////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent the parent object for garbage collection.
 */
AttributeBuilderModel::AttributeBuilderModel(QObject* parent)
	: QObject(parent), name(), index(), type(OBJECT), value()
{
}

/**
 * Destructor.
 */
AttributeBuilderModel::~AttributeBuilderModel()
{
}

/**
 * Makes a copy (clone) of this model.
 * @param cloneParent the parent object of the clone.
 * @return Returns the clone.
 */
AttributeBuilderModel* AttributeBuilderModel::clone(QObject* cloneParent) const
{
	AttributeBuilderModel* model = new AttributeBuilderModel(cloneParent);
	model->setName(name);
	model->setIndex(index);
	model->setType(type);
	model->setValue(value);
	return model;
}

/**
 * Returns the name of the new configuration attribute.
 */
const QString& AttributeBuilderModel::getName() const
{
	return name;
}

/**
 * Returns the index of the new configuration attribute.
 */
const QString& AttributeBuilderModel::getIndex() const
{
	return index;
}

/**
 * Returns the type of new attribute.
 */
AttributeTypeId AttributeBuilderModel::getType() const
{
	return type;
}

/**
 * Returns the value of the new attribute.
 */
const QString& AttributeBuilderModel::getValue() const
{
	return value;
}

/**
 * Changes the name of the new attribute.
 */
void AttributeBuilderModel::setName(const QString& newName)
{
	name = newName;
}

/**
 * Changes the index of the new attribute.
 */
void AttributeBuilderModel::setIndex(const QString& newIndex)
{
	index = newIndex;
}

/**
 * Changes the type of the new attribute.
 */
void AttributeBuilderModel::setType(AttributeTypeId newType)
{
	type = newType;
}

/**
 * Changes the value of the new attribute.
 */
void AttributeBuilderModel::setValue(const QString& newValue)
{
	value = newValue;
}

/**
 * Instanciates and returns the new attribute.
 */
Attribute* AttributeBuilderModel::createAttribute() const
{
	Attribute* attribute;
	if(type == OBJECT)
	{
		attribute = new Object(name.toStdString(), index.toStdString());
	}
	else
	{
		attribute = new Scalar(name.toStdString(), index.toStdString());
		switch(type)
		{
		case INT_SCALAR:attribute->setIntValue(value.toStdString());
			break;
		case HEX_SCALAR:attribute->setHexValue(value.toStdString());
			break;
		case REAL_SCALAR:attribute->setRealValue(value.toStdString());
			break;
		case BOOL_SCALAR:attribute->setBoolValue(value.toStdString());
			break;
		case STRING_SCALAR:attribute->setStringValue(value.toStdString());
			break;
		case OBJECT:// impossible, but it prevents the compiler from complaining.
			break;
		}
	}
	return attribute;
}

}

}
