/*
 * $Id$
 * @file FrameCounter.h
 *
 * Created on: 29 oct. 2012
 */

#ifndef get_daq_FrameCounter_h_INCLUDED
#define get_daq_FrameCounter_h_INCLUDED

#include "mdaq/daq/DataProcessorCore.h"
#include "mdaq/mt/Timer.h"
#include "mfm/FrameBuilder.h"

#ifdef VXWORKS
#include <IceE/IceE.h>
#include <IceE/Time.h>
#else
#include <Ice/Ice.h>
#include <IceUtil/Time.h>
#endif

namespace get {
namespace daq {

class FrameCounter;

class PerfFrameRecvTimer : public mdaq::mt::Timer
{
public:
	PerfFrameRecvTimer(FrameCounter& netPerf) : netPerf(netPerf) {}
	virtual void handler();

private:
	FrameCounter& netPerf;
};
/**
 * A data processor that reconstructs MFM frames and simply counts them.
 */
class FrameCounter : public virtual ::mdaq::daq::DataProcessorCore, public virtual ::mfm::FrameBuilder
{
public:
	FrameCounter(int printPeriod_s = 5);
	virtual ~FrameCounter();

	/// @name Methods overloading mdaq::daq::DataProcessorCore methods
	///@{
	virtual void resetData();
	virtual void setActive(bool _isActive = true);
	virtual void processData(const uint8_t* buffer, uint32_t nBytes); ///< Overrides method from mdaq::daq::DataProcessorCore
	///@}

	/// @name Methods overloading mfm::FrameBuilder methods
	/// @{
	virtual void processFrame(mfm::Frame& frame); ///< Overrides method inherited from mdaq::daq::DataProcessorCore
	/// @}

	void startTimer();
	void stopTimer();

protected:
	int printPeriod_ms;
	PerfFrameRecvTimer perfRcveTimer;
	double bytesReceived;
	size_t framesReceived;
	IceUtil::Time startTime;
	static const char* units[];
};

} /* namespace daq */
} /* namespace get */
#endif /* get_daq_FrameCounter_h_INCLUDED */
