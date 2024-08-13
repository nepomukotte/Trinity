#include "MyComboBox.h"

namespace CCfg
{
namespace Ui
{

////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::MyComboBox class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 * @param parent the parent widget.
 */
MyComboBox::MyComboBox(QWidget* parent)
 : QComboBox(parent)
{
}

/**
 * Destructor.
 */
MyComboBox::~MyComboBox()
{
}

/**
 * Returns the text of the current item.
 */
QString MyComboBox::getText() const
{
	return itemText(currentIndex());
}

/**
 * Select the item whose text is the same as the specified string.
 */
void MyComboBox::setText(const QString& text)
{
	setCurrentIndex(findText(text));
}

}
}
