/* =============================================================================
 * ForwardCtorTest.cpp
 * -----------------------------------------------------------------------------
 * class ::
 * Created on: 7 déc. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include <utl/ForwardCtor.hpp>
#include <boost/test/auto_unit_test.hpp>

using namespace utl;

struct B0
{
	short m_base;
	B0() : m_base(110) {}
};

struct B1
{
	int m_a0;
	B1(int a0) : m_a0(a0) {}
};

struct B2
{
	const short& m_a0;
	float m_a1;
	B2(const short& a0, float a1) : m_a0(a0), m_a1(a1) {}
};

template<class Base>
struct D : ForwardCtor<Base>
{
	int m_val;

    template<class Factory>
    D(const Factory& f)
        : ForwardCtor<Base>(f), m_val(15)
    {
    }
};

struct ForwardCtorFixture {};

BOOST_FIXTURE_TEST_SUITE(utl_ForwardCtor, ForwardCtorFixture)

BOOST_AUTO_TEST_CASE(test_0)
{
	D<B0> d(boost::in_place());
	BOOST_CHECK_EQUAL(d.m_base, 110);
	BOOST_CHECK_EQUAL(d.m_val, 15);
}

BOOST_AUTO_TEST_CASE(test_1)
{
	D<B1> d(boost::in_place(850));
	BOOST_CHECK_EQUAL(d.m_val, 15);
	BOOST_CHECK_EQUAL(d.m_a0, 850);
}


BOOST_AUTO_TEST_CASE(test_2)
{
	short s = 987;
	float f = 4.5;
	D<B2> d(boost::in_place(s, f));
	BOOST_CHECK_EQUAL(d.m_val, 15);
	BOOST_CHECK_EQUAL(d.m_a0, s);
	BOOST_CHECK_EQUAL(&d.m_a0, &s);
	BOOST_CHECK_EQUAL(d.m_a1, f);
}

BOOST_AUTO_TEST_SUITE_END()

template class D<B0>;
template class D<B1>;
template class D<B2>;
