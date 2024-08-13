/**
 * @file CoBoPLLManager.h
 * @date Sep 26, 2013
 * @author sizun
 * 
 * @note SVN tag: $Id: CoBoPLLManager.h 1726 2017-11-14 11:50:40Z psizun $
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

#ifndef get_rc_CoBoPLLManager_h_INCLUDED
#define get_rc_CoBoPLLManager_h_INCLUDED

#include <mdaq/hw/Control.h>
#include <get/rc/Node.h>
#include "get/rc/NodeId.h"
#include <boost/cstdint.hpp>
#include <map>

namespace get
{
namespace rc
{
//__________________________________________________________________________________________________
/**
 * Interface for configuring the M-CoBo PLL (Phase-Locked Loop) device.
 *
 * See emails from F. A-N (Sep 24th, 2013; Sep 27th, 2013; Dec 4th, 2013; Dec 5th, 2013; April 18th, 2014; March 17th, 2015)
 * See Texas Instruments LMK04806 Data Sheet
 */
class CoBoPLLManager
{
public:
	CoBoPLLManager(mdaq::hw::DevicePrx pll=0, mdaq::hw::DevicePrx ctrl=0);
	CoBoPLLManager(get::rc::Node & hNode);
	typedef std::map< Ice::Short, Ice::Int> RegisterMap;
	void loadPLLConfig(const float ckwFreq_MHz, const RegisterMap & registers = RegisterMap());
	void writePLLRegister(const uint32_t & value);
	void resetPLL();
	void setPllDevice(mdaq::hw::DevicePrx pll) { this->pll_ = pll; }
private:
	uint32_t firmwareDate();
	void setCKWFrequency(const float & ckwFreq_MHz, const RegisterMap & r);
	void setRegisters(const RegisterMap & r);
	void setMode(const uint8_t & mode);
	void setSyncPolInv(const bool isActiveHigh);
	void setStatus_CLKin0_MUX(const uint8_t & value);
	void resetWriteFIFO();
	void writeToFIFO(const uint32_t & data);
	bool isWriteFIFOFull();
	uint32_t writeFIFOVacancy();

	void resetReadFIFO();
	NodeId nodeId_;
	mdaq::hw::DevicePrx pll_; ///< Proxy to CoBo PLL device.
	mdaq::hw::DevicePrx ctrl_; ///< Proxy to CoBo control device.

	static const uint32_t RDFIFO_RESET = UINT32_C(0x0000000A); ///< read packet fifo reset
	static const uint32_t WRFIFO_RESET = UINT32_C(0x0000000A); ///< write packet fifo reset
	static const uint32_t PLL_RESET_MASK = UINT32_C(0x00020000); ///< R0[17] is the RESET bit
	static const size_t LMK04800_REG_NUM = 26;
	static const uint32_t PLL_DEFAULT_CONFIG[LMK04800_REG_NUM]; ///< Default PLL configuration.
	uint32_t pllConfig_[LMK04800_REG_NUM];
};
//__________________________________________________________________________________________________
} /* namespace rc */
} /* namespace get */
#endif /* get_rc_CoBoPLLManager_h_INCLUDED */
