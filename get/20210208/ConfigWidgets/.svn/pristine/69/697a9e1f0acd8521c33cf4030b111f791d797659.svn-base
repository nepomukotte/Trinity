#include "AbstractFormatter.h"
#include "Utils/Error/Exception.h"

namespace Utils
{
namespace Error
{
////////////////////////////////////////////////////////////////////////////////
// Utils::Error::AbstractFormatter class methods body
////////////////////////////////////////////////////////////////////////////////

/**
 * Analyses the exception to create a title for the error.
 * If can return an empty string if the available informations are unsufficient
 * to provide something meaningful.
 */
std::string AbstractFormatter::formatTitle(std::exception& e)
{
	Utils::Error::Exception* e2 = dynamic_cast<Utils::Error::Exception*>(&e);
	if(e2 != 0)
	{
		return e2->getTitle();
	}
	else
	{
		return std::string();
	}
}

/**
 * Analyses the exception to create a brief (one line) description of the error.
 * If can return an empty string if the available informations are unsufficient
 * to provide something meaningful.
 */
std::string AbstractFormatter::formatBrief(std::exception& e)
{
	Utils::Error::Exception* e2 = dynamic_cast<Utils::Error::Exception*>(&e);
	if(e2 != 0)
	{
		return e2->getBrief();
	}
	else
	{
		return std::string();
	}
}

/**
 * Analyses the exception to create a detailed description of the error.
 * It can return an empty string if details are not available, or if the brief
 * description is already telling everything there is to say about the error.
 */
std::string AbstractFormatter::formatDetails(std::exception& e)
{
	Utils::Error::Exception* e2 = dynamic_cast<Utils::Error::Exception*>(&e);
	if(e2 != 0)
	{
		return e2->getDetails();
	}
	else
	{
		return e.what();
	}
}

/**
 * Analyses the exception to create a suggestion for the user to resolve the
 * error.
 * It can return an empty string if we don't know how to resolve the error.
 */
std::string AbstractFormatter::formatSuggestion(std::exception& e)
{
	Utils::Error::Exception* e2 = dynamic_cast<Utils::Error::Exception*>(&e);
	if(e2 != 0)
	{
		return e2->getSuggestion();
	}
	else
	{
		return std::string();
	}
}

}
}
