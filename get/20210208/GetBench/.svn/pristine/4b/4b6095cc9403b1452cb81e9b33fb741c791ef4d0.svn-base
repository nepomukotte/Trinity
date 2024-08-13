/* =====================================================================================================================
 * DataRouter.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class get::daq::DataRouter
 * Created on: Jun 20, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * (c) Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the GetBench package implementing the GET-specific
 * software modules for the MDaq framework.
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
 * You can obtain this software from CEA under other licensing terms for commercial purposes. For this you will need to
 * negotiate a specific contract with a legal representative of CEA.
 * =====================================================================================================================
 */

#ifndef get_daq_DataRouter_h_INCLUDED
#define get_daq_DataRouter_h_INCLUDED

#include "mdaq/daq/RunProcessor.h"
#include "utl/net/SocketAddress.h"
#include "mdaq/daq/DataReceiver.h"
#include <memory>

// Forward declarations
namespace mdaq
{
namespace daq
{
class DataProcessorCore;
} // daq
} // mdaq
//_________________________________________________________________________________________________
namespace get
{
namespace daq
{
//_________________________________________________________________________________________________
class DataRouter : public virtual mdaq::daq::RunProcessor
{
public:
	DataRouter(const ::utl::net::SocketAddress& flowEndpoint, const std::string& flowType, const std::string& dataProcessorType);
	DataRouter(const ::utl::net::SocketAddress& flowEndpoint, const std::string& flowType, std::auto_ptr < ::mdaq::daq::DataProcessorCore > dataProcessor);
	virtual ~DataRouter();
    void runConfig(const std::string& configString, const ::Ice::Current& = ::Ice::Current());
    void runStart(const ::Ice::Current& = ::Ice::Current());
    void runStop(const ::Ice::Current& = ::Ice::Current());

private:
	DataRouter(DataRouter const & r);
	DataRouter& operator=(const DataRouter & r);
	void createDataReceiver(const ::utl::net::SocketAddress& flowEndpoint, const std::string& flowType);

protected:
	std::auto_ptr < mdaq::daq::DataReceiver > dataReceiver;
};
//_________________________________________________________________________________________________
} // namespace daq
} // namespace get

#endif /* get_daq_DataRouter_h_INCLUDED */
