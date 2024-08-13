/* =====================================================================================================================
 * EccCmd.h
 * ---------------------------------------------------------------------------------------------------------------------
 * class mdaq::EccCmd
 * class mdaq::EccEngine
 * Created on: Jan 5, 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * ---------------------------------------------------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * ---------------------------------------------------------------------------------------------------------------------
 * Contributors: Shebli Anvar (shebli.anvar@cea.fr)
 *               Joel Chavas (joel.chavas@cea.fr)
 * ---------------------------------------------------------------------------------------------------------------------
 * This software is part of the MDaq data acquisition framework @ CEA-Irfu.
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

#ifndef mdaq_EccCmd_h_INCLUDED
#define mdaq_EccCmd_h_INCLUDED

#include "mdaq/cli/ClientCmd.hpp"
#include "mdaq/Ecc.h"
#include "utl/Exception.h"
#include <boost/cstdint.hpp>
#include <iomanip>


namespace mdaq
{

template <class Proxy = EccPrx>
class EccCmd : public cli::ClientCmd<Proxy>
{
public:
	class Engine : public cli::ClientCmd<Proxy>::Engine
	{
	public:
		Engine(const std::string & name = "MDaq", const std::string & servantName = "Ecc");
		virtual typename cli::ClientCmd<Proxy>::IceProxy & servant(); // overridden
	};

public:
	EccCmd(Engine& engine, const std::string& name, const std::string& argsDescription, const std::string& comment)
		: cli::ClientCmd<Proxy>(engine, name, argsDescription, comment) {}
};

} // namespace mdaq

#include "mdaq/EccCmd.hpp"

#endif /* mdaq_EccCmd_h_INCLUDED */

