/**
 * @file MutantLinkManagerI.cpp
 * @date Sep 30, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: MutantLinkManagerI.cpp 936 2013-10-03 13:57:50Z psizun $
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

#include "get/cobo/MutantLinkManagerI.h"
#include "utl/Logging.h"
#ifdef HAVE_XGPIO
#include "xparameters.h"
#endif
#include <sstream>

namespace get
{
namespace cobo
{
//__________________________________________________________________________________________________
/**
 * Constructs a Mutant link manager (no access to XGPIO device yet).
 */
MutantLinkManagerI::MutantLinkManagerI()
{
}
//__________________________________________________________________________________________________
/**
 * Initializes the XGPIO device driver.
 */

bool MutantLinkManagerI::initializeXgpioDevice()
{
#ifdef HAVE_XGPIO
	int status = XGpio_Initialize(&GpioOutput, XPAR_CM_MODES_DEVICE_ID);
	if (status != XST_SUCCESS)
	{
		LOG_ERROR() << "XGpio_Initialize FAILED.";
		return false;
	}

	// Set the direction for all signals to be outputs
	XGpio_SetDataDirection(&GpioOutput, 1, 0x0);
#endif
	return true;
}
//__________________________________________________________________________________________________
void MutantLinkManagerI::setupMutantLink(MutantLinkMode mode, const Ice::Current & /* c */)
{
	static bool gpioInitOk = false;
	if (not gpioInitOk)
	{
		gpioInitOk = initializeXgpioDevice();
	}
#ifdef HAVE_XGPIO
	const unsigned channel = 1;
	u32 bits;
	std::string msg = "CoBo-Mutant link mode: ";
	if (mode == Normal)
	{
		bits = 0x0; // 0000
		LOG_INFO() << msg << "Entering normal mode of operation...";
	}
	else if (mode == Disconnected)
	{
		bits = 0xA; // 1010
		LOG_INFO() << msg << "CoBo-Mutant Link is disconnected...";
	}
	else if (mode == GMC)
	{
		bits = 0xC; // 1100 GMC FB
		LOG_INFO() << msg << "Entering GMC feedback mode...";
	}
	else if (mode == WSCA)
	{
		LOG_INFO() << msg << "Entering WSCA feedback mode...";
		bits = 0x3; // 0011 WSCA FB
	}
	else
	{
		LOG_ERROR() << msg << "Unknown mode " << mode;
		return;
	}
	XGpio_DiscreteWrite(&GpioOutput, channel, bits);
#else
	std::ostringstream errorMsg;
	errorMsg << "Could not set mode of CoBo-Mutant link to " << mode << ": no software access to XGPIO device!";
	LOG_ERROR() << errorMsg.str();
	throw ::mdaq::utl::CmdException(errorMsg.str());
#endif
}
//__________________________________________________________________________________________________
} /* namespace cobo */
} /* namespace get */
