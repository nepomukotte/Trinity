/*
 * $Id: RemoteRegisterProcessor.h 1780 2016-05-03 15:58:40Z anvar $
 * @file RemoteRegisterProcessor.h
 * @created May 3, 2016
 * @author anvar
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the hardware access classes for embedded systems of the Mordicus Real-Time Software
 * Development Framework.
 * ---------------------------------------------------------------------------------------------------------------------
 * FREE SOFTWARE LICENCING
 * This software is governed by the CeCILL license under French law and abiding by the rules of distribution of free
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef mdaq_hw_RemoteRegisterProcessor_h_INCLUDED
#define mdaq_hw_RemoteRegisterProcessor_h_INCLUDED

#include "mdaq/hw/Control.h"
#include "mdaq/hw/RegisterProcessor.h"
#include <memory>

namespace mdaq
{
namespace hw
{

class RemoteRegisterProcessor : public RegisterProcessor
{
public:
	class SynchCore : public RegisterProcessor::Core
	{
	public:
		SynchCore(const DevicePrx& devicePrx);
		void impl_w(const std::string& regName, const std::string& fieldName, Ice::Long value);
		void impl_r(const std::string& regName, const std::string& fieldName, Ice::Long& value);

	public:
		DevicePrx devicePrx;
	};

	class BatchCore : public RegisterProcessor::BatchCore
	{
	public:
		BatchCore(const DevicePrx& devicePrx);
		void impl_init();
		void impl_flush();

	public:
		DevicePrx devicePrx;
		RegCmdSeq regCmdBatchStatus;
	};

public:
	RemoteRegisterProcessor(const DevicePrx& devicePrx, bool isBatched_ = false);
	~RemoteRegisterProcessor();
	bool isBatched() { return dynamic_cast<BatchCore*>(&core()) != 0; }
	DevicePrx& devicePrx();
};

} /* namespace hw */
} /* namespace mdaq */

#endif /* mdaq_hw_RemoteRegisterProcessor_h_INCLUDED */
