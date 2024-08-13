#include "TypeDisplayInfo.h"
#include "CCfg/Ui/MyComboBox.h"
#include "CCfg/Interval.h"
#include "CCfg/Enum.h"
#include <QSpinBox>
#include <QLineEdit>
#include <limits>

namespace CCfg
{
using Ui::MyComboBox;


namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::TypeDisplayInfo class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * This private structure declares informations describing a type display informations.
 */
struct TypeDisplayData
{
	const char* const regexp;
	const char* const mask;
	const char* const editorProperty;
	Qt::AlignmentFlag alignment;
};

const TypeDisplayData types[VALUETYPE_COUNT] =
{
	{"", "", "text", Qt::AlignRight},
	{"0x[0-9A-F]{1,8}\\s*", ">\\0\\xHhhhhhhh", "text", Qt::AlignRight},
	{"[-+]?([0-9]*\\.[0-9]+|[0-9]+)([eE][-+]?[0-9]+)?", "", "text", Qt::AlignRight},
	{"", "", "text", Qt::AlignRight},
	{"", "", "text", Qt::AlignLeft}
};

////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::TypeDisplayInfo class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Returns the regular expression that validates an edited configuration value.
 * It is meant to be used by a QRegexpValidator for example.
 * @param id the id of the value type.
 * @return Returns the regular expression.
 */
const char* TypeDisplayInfo::getValidationRegexp(ValueTypeId id)
{
	return types[id].regexp;
}

/**
 * Returns the input mask used to enter a value.
 * It complies to the format defined by QLineEdit::setInputMask.
 * @param id the id of the value type.
 * @return Returns the mask.
 */
const char* TypeDisplayInfo::getInputMask(ValueTypeId id)
{
	return types[id].mask;
}

/**
 * Returns the name of the property used to get and set a value of a specific
 * type in its editor.
 * @param id the id of the value type.
 * @return Returns the name of the property.
 */
const char* TypeDisplayInfo::getEditorPropertyName(ValueTypeId id)
{
	return types[id].editorProperty;
}

/**
 * Creates and returns a suitable widget for editing values of a specific type.
 * @param id the id of the value type.
 * @param range the range  of authorized values for the editor.
 * @param parent the parent of the created widget.
 */
QWidget* TypeDisplayInfo::createEditor(ValueTypeId id, const Range* range, QWidget* parent)
{
	const Enum* const enumeration = dynamic_cast<const Enum*>(range);
	if(enumeration != 0)
	{
		// uses a combo box, whatever the type
		MyComboBox* comboBox = new MyComboBox(parent);
		const Enum::Set& values = enumeration->getValues();
		for(Enum::Set::const_iterator iVal = values.begin(); iVal != values.end(); ++iVal)
		{
			comboBox->addItem(iVal->getAsString().c_str());
		}
		return comboBox;
	}
	else
	{
		return (*creators[id])(parent, dynamic_cast<const Interval*>(range));
	}
}

/**
 * Returns the horizontal alignment that should be used to display a specific
 * type in a spreadsheet-like table.
 * @param id the id of the value type.
 * @return Returns the alignment.
 */
Qt::AlignmentFlag TypeDisplayInfo::getDisplayHorizAlignment(ValueTypeId id)
{
	return types[id].alignment;
}

template<>
QWidget* TypeDisplayInfo::createValueEditor<Int>(QWidget* parent, const Interval* interval)
{
	QSpinBox* editor = new QSpinBox(parent);
	Int min, max;
	if(interval != 0)
	{
		if(interval->getLowerBoundMode() == Interval::Unbounded)
		{
			min = std::numeric_limits<Int>::min();
		}
		else
		{
			min = interval->getLowerBound();
			if(interval->getLowerBoundMode() == Interval::ExcludeBound)		++min;
		}

		if(interval->getUpperBoundMode() == Interval::Unbounded)
		{
			max = std::numeric_limits<Int>::max();
		}
		else
		{
			max = interval->getUpperBound();
			if(interval->getUpperBoundMode() == Interval::ExcludeBound)		--max;
		}
	}
	else
	{
		min = std::numeric_limits<Int>::min();
		max = std::numeric_limits<Int>::max();
	}
	editor->setRange(min, max);
	editor->setValue(atoi(TypeInfo::getDefaultValue(INTEGER)));
	return editor;
}

template<>
QWidget* TypeDisplayInfo::createValueEditor<Hex>(QWidget* parent, const Interval* /*interval*/)
{
	QLineEdit* editor = new QLineEdit(parent);
	setRestrictions<Hex>(editor);
	editor->setText(TypeInfo::getDefaultValue(HEXADECIMAL));
	return editor;
}

template<>
QWidget* TypeDisplayInfo::createValueEditor<Real>(QWidget* parent, const Interval* /*interval*/)
{
	QLineEdit* editor = new QLineEdit(parent);
	setRestrictions<Real>(editor);
	editor->setText(TypeInfo::getDefaultValue(REAL));
	return editor;
}

template<>
QWidget* TypeDisplayInfo::createValueEditor<Bool>(QWidget* parent, const Interval* interval)
{
	MyComboBox* editor = new MyComboBox(parent);
	if(interval == 0 or interval->contains(Value(false)))
	{
		editor->addItem("false", false);
	}
	if(interval == 0 or interval->contains(Value(true)))
	{
		editor->addItem("true", true);
	}
	return editor;
}

template<>
QWidget* TypeDisplayInfo::createValueEditor<String>(QWidget* parent, const Interval* /*interval*/)
{
	QLineEdit* editor = new QLineEdit(parent);
	editor->setText(TypeInfo::getDefaultValue(STRING));
	return editor;
}

/**
 * Defines validators, input masks, etc. on types that defines such restrictions.
 */
template<typename T>
void TypeDisplayInfo::setRestrictions(QLineEdit* editor)
{
	const char* const regexp = TypeDisplayInfo::getValidationRegexp(Id<T>::val);
	const char* const mask = TypeDisplayInfo::getInputMask(Id<T>::val);

	if(regexp != 0 and regexp[0] != '\0')
	{
		editor->setValidator(new QRegExpValidator(QRegExp(regexp), editor));
	}
	if(mask != 0 and mask[0] != '\0')
	{
		editor->setInputMask(mask);
	}
}

TypeDisplayInfo::TypeEditionWidgetCreator
TypeDisplayInfo::creators[VALUETYPE_COUNT] =
{
		&TypeDisplayInfo::createValueEditor<Int>,
		&TypeDisplayInfo::createValueEditor<Hex>,
		&TypeDisplayInfo::createValueEditor<Real>,
		&TypeDisplayInfo::createValueEditor<Bool>,
		&TypeDisplayInfo::createValueEditor<String>
};


}
}
