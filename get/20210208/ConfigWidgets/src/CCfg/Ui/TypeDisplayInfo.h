#ifndef CCFG_UI_TYPEDISPLAYINFO_H_INCLUDED
#define CCFG_UI_TYPEDISPLAYINFO_H_INCLUDED

#include "CCfg/TypeInfo.h"
#include <Qt>
#include <memory>

// forward declarations
class QWidget;
class QLineEdit;

namespace CCfg
{
// forward declarations
class Range;
class Interval;

namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class provides informations about how to display and edit configuration
 * values according to their type.
 *
 * @author fchateau
 */
class TypeDisplayInfo
{
	typedef QWidget* (*TypeEditionWidgetCreator)(QWidget*, const Interval*);
	static TypeEditionWidgetCreator creators[VALUETYPE_COUNT];

	TypeDisplayInfo();
public:
	~TypeDisplayInfo();

	static const char* getValidationRegexp(ValueTypeId id);
	static const char* getInputMask(ValueTypeId id);
	static const char* getEditorPropertyName(ValueTypeId id);
	static QWidget* createEditor(ValueTypeId id, const Range* range, QWidget* parent);
	static Qt::AlignmentFlag getDisplayHorizAlignment(ValueTypeId id);

private:
	template<typename T>
	static QWidget* createValueEditor(QWidget* parent, const Interval* interval);

	template<typename T>
	static void setRestrictions(QLineEdit* editor);
};

}
}

#endif // CCFG_UI_TYPEDISPLAYINFO_H_INCLUDED
