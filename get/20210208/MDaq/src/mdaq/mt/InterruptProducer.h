/* =============================================================================
 * InterruptProducer.h
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

#ifndef mdaq_mt_InterruptProducer_INCLUDED
#define mdaq_mt_InterruptProducer_INCLUDED

#include "Semaphore.h"
#include "SysClock.h"
#include "EventProducer.h"
#include <intLib.h>
#include <iv.h>
#include <stdexcept>

namespace mdaq
{
namespace mt
{
/**
 * This class provides the basic services for notifying a polling timer
 */
class InterruptProducer : public EventProducer
{
public:
	InterruptProducer(Semaphore& notifySema, uint32_t baseAddress, int code);
	virtual ~InterruptProducer();
	virtual void start();
	virtual void stop();
	bool verifyInterruptSource();
	void releaseInterruptSource();
private:
	static void dispatchHandler(InterruptProducer* intProd);
	void     interruptEnable();
	uint32_t baseAddress_;
	int      code_;
	static const uint32_t baseAddressUndefined_ = (uint32_t)-1;
	static const int codeUndefined_ = (uint32_t)-1;
};

}
}

/***************************************************************************/
/***************************** Xilinx **************************************/
/***************************************************************************/

/***************************** Include Files *******************************/

//#include "xbasic_types.h"
//#include "xstatus.h"
//#include "xio.h"

/*****************************************************************************/
/***************** Macros (Inline Functions) Definitions *********************/
/*****************************************************************************/

/* The following macro is specific to the GNU compiler and PowerPC family. It
 * performs an EIEIO instruction such that I/O operations are synced correctly.
 * This macro is not necessarily portable across compilers since it uses
 * inline assembly.
 */
#if defined __GNUC__
#  define SYNCHRONIZE_IO __asm__ volatile ("eieio") /* should be 'mbar' ultimately */
#elif defined __DCC__
#  define SYNCHRONIZE_IO __asm volatile(" eieio")   /* should be 'mbar' ultimately */
#else
#  define SYNCHRONIZE_IO
#endif

#define XIo_In32(InputPtr) (*(volatile unsigned long *)(InputPtr)); SYNCHRONIZE_IO;
#define XIo_Out32(OutputPtr, Value) \
    { (*(volatile unsigned long *)(OutputPtr) = Value); SYNCHRONIZE_IO; }


/**
 *
 * Write a value to a INTR register. A 32 bit write is performed.
 * If the component is implemented in a smaller width, only the least
 * significant data is written.
 *
 * @param   BaseAddress is the base address of the INTR device.
 * @param   RegOffset is the register offset from the base to write to.
 * @param   Data is the data written to the register.
 *
 * @return  None.
 *
 * @note
 * C-style signature:
 * 	void INTR_mWriteReg(Xuint32 BaseAddress, unsigned RegOffset, Xuint32 Data)
 *
 */
#define INTR_mWriteReg(BaseAddress, RegOffset, Data) \
 	XIo_Out32((BaseAddress) + (RegOffset), (unsigned long)(Data))

/**
 *
 * Read a value from a INTR register. A 32 bit read is performed.
 * If the component is implemented in a smaller width, only the least
 * significant data is read from the register. The most significant data
 * will be read as 0.
 *
 * @param   BaseAddress is the base address of the INTR device.
 * @param   RegOffset is the register offset from the base to write to.
 *
 * @return  Data is the data from the register.
 *
 * @note
 * C-style signature:
 * 	Xuint32 INTR_mReadReg(Xuint32 BaseAddress, unsigned RegOffset)
 *
 */
#define INTR_mReadReg(BaseAddress, RegOffset) \
 	XIo_In32((BaseAddress) + (RegOffset))


/*****************************************************************************/
/***************** INTC SPECIFIC *********************************************/
/*****************************************************************************/

/**
 * User Logic Slave Space Offsets
 * -- SLV_REG0 : user logic slave module register 0
 */
#define INTR_USER_SLV_SPACE_OFFSET (0x00000000)
#define INTR_SLV_REG0_OFFSET (INTR_USER_SLV_SPACE_OFFSET + 0x00000000)

/**
 * Interrupt Controller Space Offsets
 * -- INTR_DISR  : device (peripheral) interrupt status register
 * -- INTR_DIPR  : device (peripheral) interrupt pending register
 * -- INTR_DIER  : device (peripheral) interrupt enable register
 * -- INTR_DIIR  : device (peripheral) interrupt id (priority encoder) register
 * -- INTR_DGIER : device (peripheral) global interrupt enable register
 * -- INTR_ISR   : ip (user logic) interrupt status register
 * -- INTR_IER   : ip (user logic) interrupt enable register
 */
#define INTR_CNTRL_SPACE_OFFSET (0x00000100)
#define INTR_DISR_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x00000000)
#define INTR_DIPR_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x00000004)
#define INTR_DIER_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x00000008)
#define INTR_DIIR_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x00000018)
#define INTR_DGIER_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x0000001C)
#define INTR_IPISR_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x00000020)
#define INTR_IPIER_OFFSET (INTR_CNTRL_SPACE_OFFSET + 0x00000028)

/**
 * Interrupt Controller Masks
 * -- INTR_TERR_MASK : transaction error
 * -- INTR_DPTO_MASK : data phase time-out
 * -- INTR_IPIR_MASK : ip interrupt requeset
 * -- INTR_RFDL_MASK : read packet fifo deadlock interrupt request
 * -- INTR_WFDL_MASK : write packet fifo deadlock interrupt request
 * -- INTR_IID_MASK  : interrupt id
 * -- INTR_GIE_MASK  : global interrupt enable
 * -- INTR_NOPEND    : the DIPR has no pending interrupts
 */
#define INTR_TERR_MASK (0x00000001UL)
#define INTR_DPTO_MASK (0x00000002UL)
#define INTR_IPIR_MASK (0x00000004UL)
#define INTR_RFDL_MASK (0x00000020UL)
#define INTR_WFDL_MASK (0x00000040UL)
#define INTR_IID_MASK (0x000000FFUL)
#define INTR_GIE_MASK (0x80000000UL)
#define INTR_NOPEND (0x80)


#endif // mdaq_mt_InterruptProducer_INCLUDED
