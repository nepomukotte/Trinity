/* =============================================================================
 * NetTestData.cpp
 * -----------------------------------------------------------------------------
 * class acq::NetTestData
 * Created on: 4 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the MDaq2 project.
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

#include "NetTestData.h"

NetTestData::NetTestData()
	: dataset1(50),
	  dataset2(20000),
	  dataset3(1200),
	  dataset4(255, '\0'),
	  dataset5(),
	  dataset6()
{
	for(size_t i = 0; i < dataset1.size(); ++i)
	{
		dataset1[i] = i;
	}
	for(size_t i = 0; i < dataset2.size(); ++i)
	{
		dataset2[i] = dataset2.size() - i;
	}
	for(size_t i = 0; i < dataset3.size(); ++i)
	{
		dataset3[i] = dataset3.size()/2 - i;
	}
	for(size_t i = 0; i < 255; ++i)
	{
		dataset4[i] = 127 - i;
	}
	for(size_t i = 0; i < sizeof(dataset5)/sizeof(dataset5[0]); ++i)
	{
		dataset5[i] = i;
	}
	for(size_t i = 0; i < dataset6.size(); ++i)
	{
		dataset6[i] = (i * i) / 9.87654321;
	}
}

NetTestData::~NetTestData()
{
}
