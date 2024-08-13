/* =============================================================================
 * RoundingTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 14 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
 * Author: fchateau
 * -----------------------------------------------------------------------------
 * © Commissariat a l'Energie Atomique et aux Energies Alternatives (CEA)
 * -----------------------------------------------------------------------------
 * Contributors: fchateau
 * -----------------------------------------------------------------------------
 * This software is part of the Utilities project.
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

#include <utl/Rounding.hpp>
#include <boost/test/auto_unit_test.hpp>

using namespace utl;

struct RoundingFixture
{

};

BOOST_FIXTURE_TEST_SUITE(utl_Rounding, RoundingFixture)

BOOST_AUTO_TEST_CASE(divide_down_u32)
{
	BOOST_CHECK_EQUAL(divide<RoundDown>(0u, 64u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(1024u, 64u), 16u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(1023u, 64u), 15u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(1025u, 64u), 16u);

	BOOST_CHECK_EQUAL(divide<RoundDown>(0u, 13u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(611u, 13u), 47u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(610u, 13u), 46u);

	BOOST_CHECK_EQUAL(divide<RoundDown>(0u, 1u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(120u, 1u), 120u);
	BOOST_CHECK_EQUAL(divide<RoundDown>(121u, 1u), 121u);
}

BOOST_AUTO_TEST_CASE(divide_down_i32)
{
	BOOST_CHECK_EQUAL(divide<RoundDown>(0, 32), 0);
	BOOST_CHECK_EQUAL(divide<RoundDown>(-128, 32), -4);
	BOOST_CHECK_EQUAL(divide<RoundDown>(-127, 32), -3);
	BOOST_CHECK_EQUAL(divide<RoundDown>(-129, 32), -4);

	BOOST_CHECK_EQUAL(divide<RoundDown>(-128, -32), 4);
	BOOST_CHECK_EQUAL(divide<RoundDown>(-127, -32), 3);
	BOOST_CHECK_EQUAL(divide<RoundDown>(-129, -32), 4);

	BOOST_CHECK_EQUAL(divide<RoundDown>(128, -32), -4);
	BOOST_CHECK_EQUAL(divide<RoundDown>(127, -32), -3);
	BOOST_CHECK_EQUAL(divide<RoundDown>(129, -32), -4);
}

BOOST_AUTO_TEST_CASE(divide_up_u32)
{
	BOOST_CHECK_EQUAL(divide<RoundUp>(0u, 64u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(1024u, 64u), 16u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(1023u, 64u), 16u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(1025u, 64u), 17u);

	BOOST_CHECK_EQUAL(divide<RoundUp>(0u, 13u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(611u, 13u), 47u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(610u, 13u), 47u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(612u, 13u), 48u);

	BOOST_CHECK_EQUAL(divide<RoundUp>(0u, 1u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(120u, 1u), 120u);
	BOOST_CHECK_EQUAL(divide<RoundUp>(121u, 1u), 121u);
}

BOOST_AUTO_TEST_CASE(divide_up_i32)
{
	BOOST_CHECK_EQUAL(divide<RoundUp>(0, 32), 0);

	BOOST_CHECK_EQUAL(divide<RoundUp>(611, 13), 47);
	BOOST_CHECK_EQUAL(divide<RoundUp>(610, 13), 47);
	BOOST_CHECK_EQUAL(divide<RoundUp>(612, 13), 48);

	BOOST_CHECK_EQUAL(divide<RoundUp>(-128, 32), -4);
	BOOST_CHECK_EQUAL(divide<RoundUp>(-127, 32), -4);
	BOOST_CHECK_EQUAL(divide<RoundUp>(-129, 32), -5);

	BOOST_CHECK_EQUAL(divide<RoundUp>(-128, -32), 4);
	BOOST_CHECK_EQUAL(divide<RoundUp>(-127, -32), 4);
	BOOST_CHECK_EQUAL(divide<RoundUp>(-129, -32), 5);

	BOOST_CHECK_EQUAL(divide<RoundUp>(128, -32), -4);
	BOOST_CHECK_EQUAL(divide<RoundUp>(127, -32), -4);
	BOOST_CHECK_EQUAL(divide<RoundUp>(129, -32), -5);
}

BOOST_AUTO_TEST_CASE(divide_near_u32)
{
	BOOST_CHECK_EQUAL(divide<RoundNear>(0u, 64u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(992u, 64u), 16u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(991u, 64u), 15u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(1056u, 64u), 17u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(1055u, 64u), 16u);

	BOOST_CHECK_EQUAL(divide<RoundNear>(0u, 13u), 0u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(605u, 13u), 47u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(604u, 13u), 46u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(618u, 13u), 48u);
	BOOST_CHECK_EQUAL(divide<RoundNear>(617u, 13u), 47u);
}

BOOST_AUTO_TEST_CASE(divide_near_i32)
{
	BOOST_CHECK_EQUAL(divide<RoundNear>(0, 32), 0);

	BOOST_CHECK_EQUAL(divide<RoundNear>(992, 64), 16);
	BOOST_CHECK_EQUAL(divide<RoundNear>(991, 64), 15);
	BOOST_CHECK_EQUAL(divide<RoundNear>(1056, 64), 17);
	BOOST_CHECK_EQUAL(divide<RoundNear>(1055, 64), 16);

	BOOST_CHECK_EQUAL(divide<RoundNear>(-128, 32), -4);
	BOOST_CHECK_EQUAL(divide<RoundNear>(-112, 32), -4);
	BOOST_CHECK_EQUAL(divide<RoundNear>(-111, 32), -3);
	BOOST_CHECK_EQUAL(divide<RoundNear>(-143, 32), -4);

	BOOST_CHECK_EQUAL(divide<RoundNear>(-112, -32), 4);
	BOOST_CHECK_EQUAL(divide<RoundNear>(-111, -32), 3);
	BOOST_CHECK_EQUAL(divide<RoundNear>(-143, -32), 4);

	BOOST_CHECK_EQUAL(divide<RoundNear>(112, -32), -4);
	BOOST_CHECK_EQUAL(divide<RoundNear>(111, -32), -3);
	BOOST_CHECK_EQUAL(divide<RoundNear>(143, -32), -4);
}

BOOST_AUTO_TEST_SUITE_END()
