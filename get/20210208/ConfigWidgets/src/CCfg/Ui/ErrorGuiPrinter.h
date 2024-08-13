#ifndef CCFG_UI_ERRORGUIPRINTER_H_INCLUDED
#define CCFG_UI_ERRORGUIPRINTER_H_INCLUDED

#include "Utils/Error/AbstractPrinter.h"
#include <QWidget>

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// classes declarations
////////////////////////////////////////////////////////////////////////////////

/**
 * This class displays error messages in a dialog box.
 *
 * @author fchateau
 */
class ErrorGuiPrinter : public Utils::Error::AbstractPrinter
{
	QWidget* parentWidget;	///< the parent widget for message boxes.
public:
	ErrorGuiPrinter(QWidget* parent = 0);

	void print(std::exception& e, Utils::Error::AbstractFormatter& formatter);

private:
	ErrorGuiPrinter(const ErrorGuiPrinter& rvalue);
	ErrorGuiPrinter& operator=(const ErrorGuiPrinter& rvalue);
};

}
}

#endif // CCFG_UI_ERRORGUIPRINTER_H_INCLUDED
