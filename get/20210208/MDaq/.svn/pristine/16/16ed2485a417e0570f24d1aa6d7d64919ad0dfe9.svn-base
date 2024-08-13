/*
 * @file RemoteRegisterProcessor.cpp
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


#include "mdaq/hw/RemoteRegisterProcessor.h"
#include "utl/Logging.h"

namespace mdaq
{
namespace hw
{

RemoteRegisterProcessor::RemoteRegisterProcessor(const DevicePrx& devicePrx, bool isBatched_)
:	RegisterProcessor(CorePtr(isBatched_ ?
		static_cast<Core*>(new BatchCore(devicePrx)) :
		static_cast<Core*>(new SynchCore(devicePrx)))
		)
{
}

RemoteRegisterProcessor::~RemoteRegisterProcessor()
{
}

DevicePrx& RemoteRegisterProcessor::devicePrx()
{
	BatchCore* batchCore = dynamic_cast<BatchCore*>(&core());
	if (batchCore == 0)
	{
		return dynamic_cast<SynchCore*>(&core())->devicePrx;
	}
	return batchCore->devicePrx;
}

RemoteRegisterProcessor::SynchCore::SynchCore(const DevicePrx& devicePrx) : devicePrx(devicePrx)
{
}

void RemoteRegisterProcessor::SynchCore::impl_w(const std::string& regName, const std::string& fieldName, Ice::Long value)
{
	if (fieldName.empty())
	{
		devicePrx->writeRegister(regName, value);
	}
	else
	{
		devicePrx->writeField(regName, fieldName, value);
	}
}

void RemoteRegisterProcessor::SynchCore::impl_r(const std::string& regName, const std::string& fieldName, Ice::Long& value)
{
	if (fieldName.empty())
	{
		value = devicePrx->readRegister(regName);
	}
	else
	{
		value = devicePrx->readField(regName, fieldName);
	}
}

RemoteRegisterProcessor::BatchCore::BatchCore(const DevicePrx& devicePrx) : devicePrx(devicePrx)
{
}

void RemoteRegisterProcessor::BatchCore::impl_init()
{
	RegisterProcessor::BatchCore::impl_init();
	regCmdBatchStatus.clear();
}

void RemoteRegisterProcessor::BatchCore::impl_flush()
{
	// Call remote batch execution
	devicePrx->execBatch(regCmdSeq(), regCmdBatchStatus);

	// Copy results to out parameters
	CmdSeq::iterator itOut = cmdBatch.begin();
	for (RegCmdSeq::const_iterator itStatus = regCmdBatchStatus.begin(); itStatus != regCmdBatchStatus.end(); ++itStatus)
	{
		itOut->outValue() = itStatus->value; // Setting output parameter value
		++itOut;
	}
}

} /* namespace hw */
} /* namespace mdaq */
