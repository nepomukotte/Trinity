#ifndef UTILS_ERROR_EXCEPHANDLING_H_INCLUDED
#define UTILS_ERROR_EXCEPHANDLING_H_INCLUDED

#include "Reporter.h"

#define CATCH_AND_PRINT	\
	catch(std::exception& e)\
	{\
		std::cerr << e.what() << std::endl;\
	}\
	catch(...)\
	{\
		std::cerr << "Unknown exception." << std::endl;\
	}

#define CATCH_AND_REPORT \
	catch(std::exception& e)\
	{\
		Utils::Error::Reporter::report(e);\
	}

#endif // UTILS_ERROR_EXCEPHANDLING_H_INCLUDED
