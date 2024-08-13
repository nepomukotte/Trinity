#ifndef CCFG_MODEL_ATTRIBUTEBUILDERMODEL_H_INCLUDED
#define CCFG_MODEL_ATTRIBUTEBUILDERMODEL_H_INCLUDED

#include "CCfg/TypeInfo.h"
#include <QString>
#include <QObject>

namespace CCfg
{
// forward declarations
class Attribute;

namespace Model
{

///////////////////////////////////////////////////////////////////////////////
// classes declarations
///////////////////////////////////////////////////////////////////////////////

/**
 * This class is the data model of the attribute builder.
 * It contains all the necessary data to instanciate a new configuration
 * attribute.
 *
 * @author fchateau
 */
class AttributeBuilderModel : public QObject
{
	Q_OBJECT

	QString name;	///< the name of the new configuration attribute.
	QString index;	///< the index of the new configuration attribute.
	AttributeTypeId type;	///< the type of the new attribute (object, integer, string, etc.).
	QString value;	///< the value of the new attribute (if it is a scalar).

public:
	AttributeBuilderModel(QObject* parent = 0);
	~AttributeBuilderModel();

	AttributeBuilderModel* clone(QObject* cloneParent = 0) const;

	const QString& getName() const;
	const QString& getIndex() const;
	AttributeTypeId getType() const;
	const QString& getValue() const;

	void setName(const QString& newName);
	void setIndex(const QString& newIndex);
	void setType(AttributeTypeId newType);
	void setValue(const QString& newValue);

	Attribute* createAttribute() const;

private:
	AttributeBuilderModel(const AttributeBuilderModel& model);
	AttributeBuilderModel& operator=(const AttributeBuilderModel& model);
};

}
}

#endif // CCFG_MODEL_ATTRIBUTEBUILDERMODEL_H_INCLUDED
