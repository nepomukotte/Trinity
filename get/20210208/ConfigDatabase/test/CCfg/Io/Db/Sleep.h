#ifndef UTILS_DATE_SLEEP_HPP_INCLUDED
#define UTILS_DATE_SLEEP_HPP_INCLUDED

#ifdef _WIN32
	#include <windows.h>

	inline unsigned int sleep(unsigned int sec)
	{
		Sleep(sec * 1000);
		return 0;
	}
#else
	#include <unistd.h>
#endif



#endif // UTILS_DATE_SLEEP_HPP_INCLUDED
