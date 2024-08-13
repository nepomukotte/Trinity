/* =============================================================================
 * InterruptProducer.cpp
 * -----------------------------------------------------------------------------
 * class mdaq::mt::InterruptProducer
 * Created on: 28 July 2009 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Chateau (frederic.chateau@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 *               Shebli Anvar (shebli.anvar@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the Multitasking sub-package of the Mordicus
 * Real-Time Software Development Framework.
 * -----------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding
 * by the rules of distribution of free software. You can use, modify and/or
 * redistribute the software under the terms of the CeCILL license as circulated
 * by CEA, CNRS and INRIA at the following URL: "http://www.cecill.info".
 * As a counterpart to the access to the source code and rights to copy, modify
 * and redistribute granted by the license, users are provided only with a
 * limited warranty and the software's author, the holder of the economic
 * rights, and the successive licensors have only limited liability. In this
 * respect, the user's attention is drawn to the risks associated with loading,
 * using, modifying and/or developing or reproducing the software by the user in
 * light of its specific status of free software, that may mean that it is
 * complicated to manipulate, and that also therefore means that it is reserved
 * for developers and experienced professionals having in-depth computer
 * knowledge. Users are therefore encouraged to load and test the software's
 * suitability as regards their requirements in conditions enabling the security
 * of their systems and/or data to be ensured and, more generally, to use and
 * operate it in the same conditions as regards security.
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 * -----------------------------------------------------------------------------
 * COMMERCIAL SOFTWARE LICENCING
 * You can obtain this software from CEA under other licencing terms for
 * commercial purposes. For this you will need to negotiate a specific contract
 * with a legal representative of CEA.
 * =============================================================================
 */

#ifdef VXWORKS

#include "InterruptProducer.h"
#include <stdexcept>

using std::runtime_error;

namespace mdaq
{
namespace mt
{
// =============================================================================
// mdaq::mt::InterruptProducer class methods body
// -----------------------------------------------------------------------------

/**
 * Default constructor.
 */
InterruptProducer::InterruptProducer(Semaphore& notifySema, uint32_t baseAddress, int code)
	: EventProducer(notifySema), baseAddress_(baseAddress), code_(code)
{
}

/**
 * Destructor.
 * It is maybe not a good idea to put a throw in a destructor. It may produce some memory leaks
 */
InterruptProducer::~InterruptProducer()
{
	if (code_ != codeUndefined_)
	{
		intDisable(code_);
		if ( intDisconnect(INUM_TO_IVEC(code_), reinterpret_cast<VOIDFUNCPTR>(&dispatchHandler), reinterpret_cast<int>(this) ) == ERROR )
			;
	}
}

void InterruptProducer::start()
{
	if ( (baseAddress_ == baseAddressUndefined_) || (code_ == codeUndefined_) )
		throw runtime_error("either baseAddress_ or code_ is NOT specified");

	if ( intConnect(
			INUM_TO_IVEC(code_),
			reinterpret_cast<VOIDFUNCPTR>(&dispatchHandler), reinterpret_cast<int>(this)
		) == ERROR )
		throw runtime_error("The interrupt has NOT been connected properly");
	interruptEnable();
	intEnable(code_);
}

void InterruptProducer::stop()
{
	if (code_ != codeUndefined_)
	{
		intDisable(code_);
		if ( intDisconnect(INUM_TO_IVEC(code_), reinterpret_cast<VOIDFUNCPTR>(&dispatchHandler), reinterpret_cast<int>(this) ) == ERROR )
			;
	}
}

/**
 * Dispatch the handler.
 * There is a question about where to put handler() with respect to releaseInterruptSource()
 * I have chosen for now to put it after, so that, inside handler, one gives the interruptEnable
 * to the user, after the release of the interrupt source.
 * If, on the contrary, you don't want to handle yourself the interrupt enable stuff inside your logic,
 * it may be better to have it before.
 * The only way around this dilemma would be to put the handler() before releaseInterruptSource()
 * and to add a virtual function intEnable() after releaseInterruptSource(),
 * or to render releaseInterruptSource() virtual
 */
void InterruptProducer::dispatchHandler(InterruptProducer* intProd)
{
	if (intProd->verifyInterruptSource())
	{
		intProd->releaseInterruptSource();
		intProd->handler();
		intProd->notify();
	}
}

/**
 * This is to show an example of how to toggle write back ISR to clear interrupts.
 */
void InterruptProducer::releaseInterruptSource()
{
	uint32_t IpStatus;

	/*
	* It clears the interrupt source
	* by toggle write back to the IP ISR register.
	*/
	IpStatus = INTR_mReadReg(baseAddress_, INTR_IPISR_OFFSET);
	INTR_mWriteReg(baseAddress_, INTR_IPISR_OFFSET, IpStatus);
}


/**
 *
 * It checks that the interrupt comes from the user logic
 *
 */
bool InterruptProducer::verifyInterruptSource()
{
	uint32_t InterruptStatus;

	/*
	* Get status from Device Interrupt Status Register.
	*/
	InterruptStatus = INTR_mReadReg(baseAddress_, INTR_DISR_OFFSET);

	/*
	* Verify the source of the interrupt is the user logic
	*/
	if ( (InterruptStatus & INTR_IPIR_MASK) == INTR_IPIR_MASK )
		return true;
	else
		return false;
}

/**
 *
 * Enable all possible interrupts from INTR_EXAMPLE device.
 *
 * @param   baseaddr_p is the base address of the INTR_EXAMPLE device.
 *
 * @return  None.
 *
 * @note    None.
 *
 */
void InterruptProducer::interruptEnable()
{
  /*
   * Enable all interrupt source from user logic.
   */
  INTR_mWriteReg(baseAddress_, INTR_IPIER_OFFSET, 0x00000001);

  /*
   * Enable all possible interrupt sources from device.
   */
  INTR_mWriteReg(baseAddress_, INTR_DIER_OFFSET,
    INTR_TERR_MASK
    | INTR_DPTO_MASK
    | INTR_IPIR_MASK
    );

  /*
   * Set global interrupt enable.
   */
  INTR_mWriteReg(baseAddress_, INTR_DGIER_OFFSET, INTR_GIE_MASK);
}

} //mt
} //mdaq

#endif
