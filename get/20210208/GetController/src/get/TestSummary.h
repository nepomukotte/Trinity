/**
 * $Id: TestSummary.h 403 2012-10-05 14:47:10Z psizun $
 * @file TestSummary.h
 * @date 19 avr. 2012
 * @author sizun
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
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
 * -----------------------------------------------------------------------------
 */

#ifndef get_TstDescription_h_INCLUDED
#define get_TstDescription_h_INCLUDED

#include <mdaq/utl/net/SocketAddress.h>
#include <string>
#include <ostream>

// Forward declarations
namespace CCfg {
	class Attribute;
}

namespace get {
//_____________________________________________________________________________
class TestSummary
{
public:
	TestSummary();
	TestSummary(const CCfg::CConfig & cfg);
	virtual ~TestSummary();
	void print(std::ostream & stream);
	void fromCConfig(const CCfg::CConfig & cfg);
	void toCConfig(const CCfg::CConfig & cfg) const;

    /** @name Setters
     * See corresponding attribute description.
     */
	///@{
	void setName(std::string name);
	void setType(std::string type);
	void setHwConfigPath(std::string hwConfigPath);
	void setParamConfigPath(std::string paramConfigPath);
	void setDataPath(std::string dataPath);
	///@}

    /** @name Getters
     * See corresponding attribute description.
     */
	///@{
	std::string name() const;
	std::string type() const;
	std::string hwConfigPath() const;
	std::string paramConfigPath() const;
	std::string dataPath() const;
	const mdaq::utl::net::SocketAddress & eccEndpoint() const;
	const mdaq::utl::net::SocketAddress & targetEndpoint() const;
	mdaq::utl::net::SocketAddress & eccEndpoint();
	mdaq::utl::net::SocketAddress & targetEndpoint();
	///@}
private:
	std::string name_; ///< Name of the test
	std::string type_; ///< Type of the test
	std::string hwConfigPath_; ///< Qualified path of the hardware configuration file
	std::string conditionsConfigPath_; ///< Qualified path of the test parameter file
	std::string dataPath_; ///< Data output path
	mdaq::utl::net::SocketAddress eccEndpoint_; ///< IP address and port of the ECC server
	mdaq::utl::net::SocketAddress targetEndpoint_; ///< IP address and port of the target
};
//_____________________________________________________________________________
} /* namespace get */
#endif /* get_TstDescription_h_INCLUDED */
