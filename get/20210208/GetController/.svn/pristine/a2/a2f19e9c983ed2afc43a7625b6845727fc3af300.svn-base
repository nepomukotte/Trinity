/* =============================================================================
 * TestInfo.h
 * -----------------------------------------------------------------------------
 * class get::TestInfo
 * Created on: 24 August 2010 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: Frederic Druillole (frederic.druillole@cea.fr)
 * -----------------------------------------------------------------------------
 * This software is part of the GetController project dedicated to the GET
 * data acquisition test bench. "GET" stands for General Electronics for TPC.
 * http://www-actar-get.cea.fr
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

#ifndef get_TestInfo_h_INCLUDED
#define get_TestInfo_h_INCLUDED

#include <string>

/**
 * This class contains all the information relevant for a test.
 */
class TestInfo {

private:
	std::string testName;
	std::string testType;
	std::string hwFile;
	std::string paramFile;
	std::string dataPath;
	std::string eccIp;
	std::string daqIp;
	std::string targetIp;

public:
	TestInfo();
	virtual ~TestInfo();
	//getter and setter
	void setTestName(std::string st) {testName = st;};
	void setTestType(std::string st) {testType = st;};
	void setHwFile(std::string st) {hwFile = st;};
	void setParamFile(std::string st) {paramFile = st;};
	void setDataPath(std::string st) {dataPath = st;};
	void setEccIp(std::string st) {eccIp = st;};
	void setDaqIp(std::string st) {daqIp = st;};
	void setTargetIp(std::string st) {targetIp = st;};

	std::string getTestName() const {return testName;};
	std::string getTestType() const {return testType;};
	std::string getHwFile() const {return hwFile;};
	std::string getParamFile() const {return paramFile;};
	std::string getDataPath() const {return dataPath;};
	std::string getEccIp() const {return eccIp;};
	std::string getDaqIp() const {return daqIp;};
	std::string getTargetIp() {return targetIp;};
};

#endif /* get_TestInfo_h_INCLUDED */
