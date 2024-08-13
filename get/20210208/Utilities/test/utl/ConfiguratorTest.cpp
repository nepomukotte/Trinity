/* =============================================================================
 * ConfiguratorTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 19 mars 2013 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include <utl/Configurator.h>
#include <boost/test/auto_unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <iomanip>

using namespace utl;

struct ConfiguratorFixture
{
	struct Config
	{
		virtual ~Config() {}
	};

	struct Config1 : Config
	{
		Config1(const char* initVal) : value(initVal) {}
		std::string value;
	};

	struct Config2 : Config
	{
		Config2(float x, int prec) : value(x), precision(prec) {}
		float value;
		int precision;
	};

	struct ConfigUnused : Config
	{
		float y;
	};

	struct Obj1
	{
		int x;
	};

	struct Obj2
	{
		std::string y;
	};

	struct ObjUnused
	{
		float z;
	};

	static void configure11(Obj1& obj, const Config1& config)
	{
		obj.x = boost::lexical_cast<int>(config.value);
	}

	static void configure21(Obj2& obj, const Config1& config)
	{
		obj.y = config.value;
		obj.y += ' ';
		obj.y += config.value;
	}

	static void configure12(Obj1& obj, const Config2& config)
	{
		obj.x = static_cast<int>(config.value * 100);
	}

	static void configure22(Obj2& obj, const Config2& config)
	{
		std::ostringstream stream;
		stream << std::fixed << std::setprecision(config.precision) << config.value;
		obj.y = stream.str();
	}

	ConfiguratorFixture()
	{
		configurator.registerFunction(&ConfiguratorFixture::configure11);
		configurator.registerFunction(&ConfiguratorFixture::configure12);
		configurator.registerFunction(&ConfiguratorFixture::configure21);
		configurator.registerFunction(&ConfiguratorFixture::configure22);
	}
	~ConfiguratorFixture()
	{

	}
	Configurator configurator;
};

BOOST_FIXTURE_TEST_SUITE(utl_Configurator, ConfiguratorFixture)

BOOST_AUTO_TEST_CASE(configurator_Obj1Config1)
{
	Config1 cfg("4852");
	Obj1 obj;
	configurator.configure(obj, cfg);

	BOOST_CHECK_EQUAL(obj.x, 4852);
}

BOOST_AUTO_TEST_CASE(configurator_Obj2Config1)
{
	Config1 cfg("Toto");
	Obj2 obj;
	configurator.configure(obj, cfg);
	BOOST_CHECK_EQUAL(obj.y, "Toto Toto");
}

BOOST_AUTO_TEST_CASE(configurator_Obj1Config2)
{
	Config2 cfg(158.523, 3);
	Obj1 obj;
	configurator.configure(obj, cfg);
	BOOST_CHECK_EQUAL(obj.x, 15852);
}

BOOST_AUTO_TEST_CASE(configurator_Obj2Config2)
{
	Config2 cfg(1.23456789, 6);
	Obj2 obj;
	configurator.configure(obj, cfg);
	BOOST_CHECK_EQUAL(obj.y, "1.234568");
}


BOOST_AUTO_TEST_CASE(configurator_polymo)
{
	std::auto_ptr<Config> cfg(new Config1("4852"));
	Obj1 obj;
	configurator.configure(obj, *cfg);
	BOOST_CHECK_EQUAL(obj.x, 4852);
}

BOOST_AUTO_TEST_CASE(configurator_missing_config)
{
	ConfigUnused cfg;
	Obj1 obj;
	BOOST_CHECK_THROW(configurator.configure(obj, cfg), Configurator::MissingEntry);
}

BOOST_AUTO_TEST_CASE(configurator_missing_obj)
{
	Config1 cfg("123.456");
	ObjUnused obj;
	BOOST_CHECK_THROW(configurator.configure(obj, cfg), Configurator::MissingEntry);
}

REGISTER_CONFIGURATOR(&ConfiguratorFixture::configure11)

BOOST_AUTO_TEST_CASE(configurator_singleton)
{
	Config1 cfg("8542");
	Obj1 obj;
	GlobalConfigurator::instance().configure(obj, cfg);
	configurator.configure(obj, cfg);
	BOOST_CHECK_EQUAL(obj.x, 8542);
}

BOOST_AUTO_TEST_SUITE_END()
