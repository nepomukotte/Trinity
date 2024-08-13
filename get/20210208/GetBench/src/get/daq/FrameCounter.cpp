/*
 * $Id$
 * @file FrameCounter.cpp
 * Created on: 29 oct. 2012
 */

#include "get/daq/FrameCounter.h"
#include "mfm/Serializer.h"
#include "mfm/Frame.h"
#include "mfm/FrameDictionary.h"
#include "utl/Logging.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <climits>

using std::cout;
using std::endl;
using std::fixed;
using std::string;

namespace get {
namespace daq {

const char* FrameCounter::units[] = {" b/s", " Kb/s", " Mb/s", " Gb/s"};

void PerfFrameRecvTimer::handler()
{
	netPerf.resetData();
}

FrameCounter::FrameCounter(int printPeriod_s)
	: DataProcessorCore(false),
	  printPeriod_ms(1000*printPeriod_s),
	  perfRcveTimer(*this)
{
	LOG_DEBUG() << "Creating FrameCounter " << std::hex << this << std::dec;
	resetData();
}

FrameCounter::~FrameCounter()
{
	LOG_DEBUG() << "Destroying FrameCounter";
}

void FrameCounter::startTimer()
{
	perfRcveTimer.schedule(printPeriod_ms, printPeriod_ms);
	resetData();
}

void FrameCounter::stopTimer()
{
	perfRcveTimer.cancel();
	resetData();
}

void FrameCounter::resetData()
{
	if (isActive())
	{
		double perf = bytesReceived / (IceUtil::Time::now() - startTime).toSecondsDouble();
		perf *= CHAR_BIT; // nb of bits per second
		int iunit = 3;
		if      (perf >= 1e15) { perf *= 1e-9;          }
		else if (perf >= 1e+6) { perf *= 1e-6; --iunit; }
		else if (perf >= 1e+3) { perf *= 1e-3; --iunit; }
		else                   {               --iunit; }
		if (perf > 0.0)
		{
			cout << ">>> "
				<< std::fixed << std::setw(9) << std::setprecision(1) << perf << units[iunit] << ' '
				<< std::fixed << std::setw(9) << std::setprecision(1) << framesReceived/(printPeriod_ms*0.001) << " fps" << " <<< " << endl;
		}
	}

	bytesReceived = 0.;
	framesReceived = 0;
	startTime = IceUtil::Time::now();
}

void FrameCounter::setActive(bool _isActive)
{
	if (_isActive)
	{
		startTimer();
		DataProcessorCore::setActive(true);
	}
	else
	{
		DataProcessorCore::setActive(false);
		stopTimer();
	}
}

void FrameCounter::processData(const uint8_t* buffer, uint32_t nBytes)
{
	const mfm::Byte* chunkBeg =  reinterpret_cast<const mfm::Byte*>(buffer);
	const mfm::Byte* chunkEnd =  reinterpret_cast<const mfm::Byte*>(buffer + nBytes);

	try
	{
		FrameBuilder::addDataChunk(chunkBeg, chunkEnd);
		return;
	}
	catch (const mfm::Exception& e)
	{
		// This is probably due to a data format error: dump some of it for debugging:
		LOG_ERROR() << e.what();
		LOG_ERROR() << "Dumping data bytes in hexadecimal (but no more than 8 lines, i.e. 256 Bytes): ";
		std::ios::fmtflags flags(std::cerr.flags()); // saving format flags
		std::cerr << std::hex << std::right << std::setfill('0');
		for (size_t i = 0; (i < nBytes) and (i < 8*32); ++i)
		{
			if (i % 32 == 0) std::cerr << '\n' << std::setw(3) << i << ": ";
			if (i % 8 == 0) std::cerr << ' ';
			std::cerr << std::setw(2) << (int) buffer[i] << ' ' ;
		}
		std::cerr << std::endl;
		std::cerr.flags(flags); // restoring format flags

	}
	catch (const std::exception & e)
	{
		LOG_ERROR() << e.what();
	}
	// Reset frame builder
	LOG_WARN() << "Resetting frame builder...";
	FrameBuilder::reset();
}

void FrameCounter::processFrame(mfm::Frame& frame)
{
	bytesReceived += frame.header().frameSize_B();
	framesReceived++;
	//LOG_DEBUG() << "Received " << bytesReceived/1024. << " kiB of frame(s)";
}

} /* namespace daq */
} /* namespace get */
