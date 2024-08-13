/*
 * @file RegisterProcessorTest.cpp
 * @created Feb 3, 2017
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

#include "mdaq/hw/RegisterProcessor.h"

#include <boost/test/auto_unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE(chained_test)
{
	// Instantiate a register processor with default core and internal trace vector.
	mdaq::hw::RegisterProcessor regProc;

	// Prepare output variables
	Ice::Long regVal = -1;
	Ice::Long fieldVal = -1;
	regProc.nop()
	       .sleep_ms(100)
		   .write("Toto", 123456789)
		   .write("Toto", "Age", 98765)
		   .read("Toto", regVal)
		   .read("Toto", "Age", fieldVal);

	bool isFlushed = regProc.isFlushed();
	BOOST_CHECK_EQUAL(regVal, 123456789);
	BOOST_CHECK_EQUAL(fieldVal, 98765);
	BOOST_CHECK_EQUAL(isFlushed, true);
}

BOOST_AUTO_TEST_CASE(cumulative_test)
{
	// Instantiate a register processor with default core and internal trace vector.
	mdaq::hw::RegisterProcessor regProc;

	// Prepare output variables
	Ice::Long regVal = -1;
	Ice::Long fieldVal = -1;
	regProc.nop()
	       .sleep_ms(100)
		   .write("Toto", 123456789)
		   .write("Toto", "Age", 98765)
		   .read("Toto", regVal)
		   .read("Toto", "Age", fieldVal);
	bool isFlushed = regProc.isFlushed();
	BOOST_CHECK_EQUAL(regVal, 123456789);
	BOOST_CHECK_EQUAL(fieldVal, 98765);
	BOOST_CHECK_EQUAL(isFlushed, false);
}

BOOST_AUTO_TEST_CASE(sleep_test)
{
	// Instantiate a register processor with default core and internal trace vector.
	mdaq::hw::RegisterProcessor regProc;

	// Prepare output variables
	Ice::Long regVal = -1;
	Ice::Long fieldVal = -1;

	// Measure timing without sleep
	IceUtil::Time startTime = IceUtil::Time::now();
	regProc.nop()
		   .write("Toto", 123456789)
		   .write("Toto", "Age", 98765)
		   .read("Toto", regVal)
		   .read("Toto", "Age", fieldVal);
	double deltaWithout = (IceUtil::Time::now()-startTime).toMilliSecondsDouble();

	// Measure timing with 100ms sleep
	startTime = IceUtil::Time::now();
	regProc.nop()
		   .write("Toto", 123456789)
		   .write("Toto", "Age", 98765)
		   .sleep_ms(100)
		   .read("Toto", regVal)
		   .read("Toto", "Age", fieldVal);
	double deltaWith = (IceUtil::Time::now()-startTime).toMilliSecondsDouble();
	BOOST_CHECK_CLOSE(deltaWith-deltaWithout, 100., 1.0);
}

