/**
 * @file PeriodicTask.cpp
 * @date Jan 9, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: PeriodicTask.cpp 1140 2014-02-13 14:59:46Z psizun $
 * @note Contributor(s): Patrick Sizun
 * @note 
 * @note This file is part of the GetBench software project.
 *
 * @copyright © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 *
 * @par FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding  * by the rules of distribution of free
 * software. You can use, modify and/or redistribute the software under the terms of the CeCILL license as circulated by
 * CEA, CNRS and INRIA at the following URL: "http://www.cecill.info". As a counterpart to the access to the source code
 * and rights to copy, modify and redistribute granted by the license, users are provided only with a limited warranty
 * and the software's author, the holder of the economic rights, and the successive licensors have only limited
 * liability. In this respect, the user's attention is drawn to the risks associated with loading, using, modifying
 * and/or developing or reproducing the software by the user in light of its specific status of free software, that may
 * mean that it is complicated to manipulate, and that also therefore means that it is reserved for developers and
 * experienced professionals having in-depth computer knowledge. Users are therefore encouraged to load and test the
 * software's suitability as regards their requirements in conditions enabling the security of their systems and/or data
 * to be ensured and, more generally, to use and operate it in the same conditions as regards security. The fact that
 * you are presently reading this means that you have had knowledge of the CeCILL license and that you accept its terms.
 *
 * @par COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 *
 */

#include "get/utl/PeriodicTask.h"
#include "utl/Logging.h"

#include <boost/bind.hpp>
#include <boost/cstdint.hpp>
#include <boost/asio/deadline_timer.hpp>
#include <boost/asio/placeholders.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/exponential_distribution.hpp>
using boost::posix_time::time_duration;
using boost::posix_time::milliseconds;
using boost::posix_time::microseconds;
#ifdef VXWORKS
#include <taskLib.h>
#include <vxWorks.h>
#endif

#include <cmath>
#include <cstdlib>

namespace get
{
namespace utl
{
//__________________________________________________________________________________________________
class PeriodicTask::TaskThread : public IceUtil::Thread
{
public:
	TaskThread(PeriodicTask & task, time_duration period = milliseconds(1000), bool poisson=false)
		: period_(period), poissonMode_(poisson), task_(task),
		  dist_(1.0/period.total_microseconds()), die_(gen_, dist_),
		  timer_(io_), stop_(false)
	{
		 timer_.expires_from_now(nextDelay());
	}
	void stop()
	{
		stop_ = true;
		timer_.cancel();
	}
	time_duration nextDelay() const
	{
		if (poissonMode_)
			return nextPoissonDelay();
		return period_;
	}
	time_duration nextPoissonDelay() const
	{
		return microseconds(static_cast< boost::int64_t >(std::floor(0.5 + die_())));
	}
protected:
	void run()
	{
#ifdef VXWORKS
		taskPrioritySet(taskIdSelf(), 105); // Ensures task has priority greater than interrupt monitor but lower than slow control
#endif
		timer_.async_wait(boost::bind(&TaskThread::onTimeout, this, boost::asio::placeholders::error));
		io_.run();
	}
	void onTimeout(const boost::system::error_code & e)
	{
		if (boost::asio::error::operation_aborted == e)
		{
			stop_ = true;
		}
		if (not stop_)
		{
			timer_.expires_at(timer_.expires_at() + nextDelay());
			timer_.async_wait(boost::bind(&TaskThread::onTimeout, this, boost::asio::placeholders::error));
			try
			{
				task_.doTask();
			}
			catch (const std::exception & e)
			{
				LOG_ERROR() << e.what();
			}
		}
	}
private:
	time_duration period_; ///< Period
	bool poissonMode_; ///< Flag indicating a Poisson process instead of a strictly periodic process
	PeriodicTask & task_;
	boost::mt19937 gen_; ///< Random number generator
	boost::exponential_distribution< > dist_; ///< Intervals between occurrences of a Poisson process follow an exponential distribution
	mutable boost::variate_generator< boost::mt19937 &, boost::exponential_distribution< > > die_;
	boost::asio::io_service io_;
	boost::asio::deadline_timer timer_;
	bool stop_;
};
//__________________________________________________________________________________________________
PeriodicTask::PeriodicTask(size_t period_ms, bool poisson)
	: taskThread_(new TaskThread(*this, milliseconds(period_ms), poisson)),
	  taskThreadCtrl_(taskThread_->start())
{
}
//__________________________________________________________________________________________________
PeriodicTask::~PeriodicTask()
{
	if (taskThread_)
	{
		taskThread_->stop();
		taskThreadCtrl_.join();
	}
}
//__________________________________________________________________________________________________
} /* namespace utl */
} /* namespace get */
