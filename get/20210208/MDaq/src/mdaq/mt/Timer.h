
/* =====================================================================================================================
 * Timer.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::mt::Timer
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the Multitasking sub-package of the Mordicus
 * Real-Time Software Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
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
 * ---------------------------------------------------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mdaq_mt_Timer_h_INCLUDED
#define mdaq_mt_Timer_h_INCLUDED

#include "utl/Exception.h"
#include "SigAction.h"
#include "SigEvent.hpp"

namespace mdaq
{
namespace mt
{
/**
 * This class implements a timer that calls a handler method each time it elapses. It can be periodic or one-shot
 * depending on configuration. Any actual timer must derive from this and implement the abstract handler() method
 * inherited from SigAction.
 */
class Timer : public SigAction
{
public:
	Timer(int sigNo_ = SIGALRM);
	virtual ~Timer();

	void schedule(long _duration_ms, long _period_ms);
	void schedule(long _duration_ms, int  _period_ms) { schedule((long) _duration_ms, (long) _period_ms); }
	void schedule(int  _duration_ms, long _period_ms) { schedule((long) _duration_ms, (long) _period_ms); }
	void schedule(int  _duration_ms, int  _period_ms) { schedule((long) _duration_ms, (long) _period_ms); }
	void schedule(long _duration_ms, bool _isPeriodic) { schedule((long) _duration_ms, _isPeriodic ? (long) _duration_ms : 0L); }
	void schedule(int  _duration_ms, bool _isPeriodic) { schedule((long) _duration_ms, _isPeriodic ? (long) _duration_ms : 0L); }
	void schedule(long _duration_ms) { schedule((long) _duration_ms, false); }
	void schedule(int  _duration_ms) { schedule((long) _duration_ms, false); }
	void cancel();

	long period_ms();
	long timeLeft_ms();
	bool isPeriodic() { return period_ms() != 0; }
	bool isScheduled() { return timeLeft_ms() != 0; }

private:
	SigEvent<Timer> sigEvent;
	timer_t posixTimer;
};

} // namespace mt
} // namespace mdaq

#endif // mdaq_mt_Timer_h_INCLUDED
