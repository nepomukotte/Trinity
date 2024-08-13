/**
 * @file PeriodicTask.h
 * @date Jan 9, 2014
 * @author sizun
 * 
 * @note SVN tag: $Id: PeriodicTask.h 1118 2014-02-03 09:11:50Z psizun $
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

#ifndef get_utl_PeriodicTask_h_INCLUDED
#define get_utl_PeriodicTask_h_INCLUDED

#include ICEUTIL_HEADER(Thread.h)

namespace get
{
namespace utl
{
//__________________________________________________________________________________________________
/**
 * An abstract class for periodic tasks.
 * Actual tasks need to derive PeriodicTask and implement the doTask() method.
 */
class PeriodicTask
{
public:
	PeriodicTask(size_t period_ms = 10000u, bool poisson=false);
	virtual ~PeriodicTask();
	virtual void doTask() = 0;
protected:
	class TaskThread;
private:
	IceUtil::Handle<TaskThread> taskThread_;
	IceUtil::ThreadControl  taskThreadCtrl_;
};
//__________________________________________________________________________________________________
} /* namespace utl */
} /* namespace get */
#endif /* get_utl_PeriodicTask_h_INCLUDED */
