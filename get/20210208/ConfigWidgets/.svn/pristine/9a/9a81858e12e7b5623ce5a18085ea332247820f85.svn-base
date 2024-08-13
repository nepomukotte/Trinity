#include "ErrorGuiPrinter.h"
#include <QMessageBox>
#include <sstream>

namespace CCfg
{
namespace Ui
{
////////////////////////////////////////////////////////////////////////////////
// CCfg::Ui::ErrorGuiPrinter class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Default constructor.
 */
ErrorGuiPrinter::ErrorGuiPrinter(QWidget* parent) : parentWidget(parent)
{
}

/*
 * (see inherited doc)
 */
void ErrorGuiPrinter::print(std::exception& e, Utils::Error::AbstractFormatter& formatter)
{
	std::string title = formatter.formatTitle(e);
	std::string brief = formatter.formatBrief(e);
	std::string details = formatter.formatDetails(e);
	std::string suggestion = formatter.formatSuggestion(e);

	std::ostringstream stream;
	if(not brief.empty())
	{
		stream << "<b>" << brief << "</b>";
	}
	if(not details.empty())
	{
		if(not brief.empty())
		{
			stream << "<br>";
		}
		stream << details;
	}
	if(not suggestion.empty())
	{
		stream << "<br><br><i>" << suggestion << "</i>";
	}

	QMessageBox::warning(parentWidget,
		QString::fromStdString(title),
		QString::fromStdString(stream.str()));
}

}
}
