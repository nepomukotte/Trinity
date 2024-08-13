/*
 * demoMonitor.cpp
 *
 *  Created on: Jun 21, 2010
 *      Author: jchavas
 */

#ifdef VXWORKS

#include "mdaq/mt/PollingMonitorAdapter.hpp"
#include "mdaq/mt/InterruptMonitorAdapter.hpp"
#include <iostream>
#include <cstdio>
#include <stdexcept>
#include <logLib.h>
#include <taskLib.h>
#include <unistd.h>

class Algo
{
public:
	Algo(int i) : intAlgo(i) {};
	bool isItTrue() { logMsg("Algo: isItTrue launched\n",0,0,0,0,0,0); return true;};
	void handler() {logMsg("Algo: Handler\n",0,0,0,0,0,0);}
	void actionBefore() {std::cout << "Algo: actionBefore: intAlgo = " << intAlgo << std::endl;}
	void actionAfter() {std::cout << "Algo: action: intAlgo = " << intAlgo << std::endl;}
private:
	int intAlgo;
};

Algo a(1);
/*
mdaq::mt::PollingMonitorAdapter<> daqLoop(2000000);
mdaq::mt::PollingMonitorAdapter<Algo> daqLoop(2000000, &a, &Algo::actionAfter);
mdaq::mt::PollingMonitorAdapter<Algo> daqLoop(2000000, &a, &Algo::actionBefore, &Algo::actionAfter);
mdaq::mt::PollingMonitorAdapter<Algo> daqLoop(2000000, &a, 0, &Algo::actionAfter);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(2000000, &a, &Algo::actionAfter,
		&a, &Algo::isItTrue);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(2000000, &a, &Algo::actionBefore, &Algo::actionAfter,
		&a, &Algo::isItTrue);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(2000000, &a, 0, &Algo::actionAfter,
		&a, &Algo::isItTrue);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(2000000, &a, &Algo::actionAfter,
		&a, &Algo::isItTrue, &Algo::handler);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(2000000, &a, &Algo::actionAfter,
		&a, &Algo::isItTrue, 0);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(2000000, &a, &Algo::actionBefore, &Algo::actionAfter,
		&a, &Algo::isItTrue);
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(1000000, &a, 0, &Algo::actionAfter,
		&a, &Algo::isItTrue, O);
		*/
mdaq::mt::PollingMonitorAdapter<Algo, Algo> daqLoop(1000000, &a, &Algo::actionBefore, &Algo::actionAfter,
		&a, &Algo::isItTrue, &Algo::handler);

Algo a2(2);
/*
mdaq::mt::InterruptMonitorAdapter<> daqLoop2(0x60040000, 0);
mdaq::mt::InterruptMonitorAdapter<Algo> daqLoop2(0x60040000, 0, &a2, &Algo::actionAfter);
mdaq::mt::InterruptMonitorAdapter<Algo> daqLoop2(0x60040000, 0, &a2, &Algo::actionBefore, &Algo::actionAfter);
mdaq::mt::InterruptMonitorAdapter<Algo> daqLoop2(0x60040000, 0, &a2, 0, &Algo::actionAfter);
mdaq::mt::InterruptMonitorAdapter<Algo,Algo> daqLoop2(0x60040000, 0, &a2, &Algo::actionAfter, &a2, &Algo::handler);
*/
mdaq::mt::InterruptMonitorAdapter<Algo,Algo> daqLoop2(0x60040000, 0, &a2, &Algo::actionBefore, &Algo::actionAfter,
		&a2, &Algo::handler);

int pollingMonStart()
{
	try
	{
		daqLoop.start(3000);
	}
	catch (std::runtime_error& e)
	{
		std::cout << "*** " << e.what() << std::endl;
		daqLoop.stop();
	}
	catch (...)
	{
		std::cout << "error" << std::endl;
		daqLoop.stop();
	}
	return 0;
}

int pollingMonStop()
{
	try
	{
		daqLoop.stop();
	}
	catch (std::runtime_error& e)
	{
		std::cout << "*** " << e.what() << std::endl;
	}
	catch (...)
	{
		std::cout << "error" << std::endl;
	}
	return 0;
}


int interruptMonStart()
{
	try
	{
		daqLoop2.start(20000);
	}
	catch (std::runtime_error& e)
	{
		std::cout << "*** " << e.what() << std::endl;
		daqLoop2.stop();
	}
	catch (...)
	{
		std::cout << "*** error" << std::endl;
		daqLoop2.stop();
	}
	return 0;
}

int interruptMonStop()
{
	try
	{
		daqLoop2.stop();
	}
	catch (std::runtime_error& e)
	{
		std::cout << "*** " << e.what() << std::endl;
	}
	catch (...)
	{
		std::cout << "error" << std::endl;
	}
	return 0;
}


#endif

