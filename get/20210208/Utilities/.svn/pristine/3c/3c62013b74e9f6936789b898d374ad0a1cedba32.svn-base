/* =============================================================================
 * FactoryTest.cpp
 * -----------------------------------------------------------------------------
 * Created on: 27 nov. 2012 at Irfu/Sedi/Lilas, CEA Saclay, F-91191, France
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

#include <utl/Factory.ipp>
#include <utl/Singleton.ipp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION >= 104300
#include <boost/functional/factory.hpp>
#endif
#include <memory>
#include <iostream>

using namespace utl;

enum AnimalType
{
	CAT,
	DOG,
	HORSE
};

struct Animal
{
	Animal(AnimalType t) : type(t) {}
	virtual ~Animal() {}
	AnimalType type;
};

struct Cat : Animal
{
	Cat() : Animal(CAT) {}
};

struct Dog : Animal
{
	Dog() : Animal(DOG) {}
};

struct Horse : Animal
{
	Horse() : Animal(HORSE) {}
};

enum MonsterType
{
	TROLL,
	GIANT,
	HYDRA
};

struct Monster
{
	Monster(int healthPoints, MonsterType t) : hp(healthPoints), type(t) {}
	virtual ~Monster() {}
	int hp;
	MonsterType type;
};

struct Troll : Monster
{
	Troll(int hp) : Monster(hp, TROLL) {}
};

struct Giant : Monster
{
	Giant(int hp) : Monster(hp, GIANT) {}
};

struct Hydra : Monster
{
	Hydra(int hp) : Monster(hp, HYDRA) {}
};

Animal* createCat() { return new Cat; }
Animal* createDog() { return new Dog; }
Animal* createHorse() { return new Horse; }

Monster* createTroll(int hp) { return new Troll(hp); }
Monster* createGiant(int hp) { return new Giant(hp); }
Monster* createHydra(int hp) { return new Hydra(hp); }

struct FactoryMapFixture
{
	FactoryMapFixture()
	{
#if BOOST_VERSION < 104300 || BOOST_VERSION >= 106100
		animalFactory
			.add("Dog", createDog)
			.add("Cat", createCat)
			.add("Horse", createHorse);
		monsterFactory
			.add("Troll", createTroll)
			.add("Giant", createGiant)
			.add("Hydra", createHydra);
#else
		animalFactory
			.add("Dog", boost::factory<Dog*>())
			.add("Cat", boost::factory<Cat*>())
			.add("Horse", boost::factory<Horse*>());
		monsterFactory
			.add("Troll", boost::factory<Troll*>())
			.add("Giant", boost::factory<Giant*>())
			.add("Hydra", boost::factory<Hydra*>());
#endif
	}
	FactoryMap<std::string, Animal*()> animalFactory;
	FactoryMap<std::string, Monster*(int)> monsterFactory;
};

BOOST_FIXTURE_TEST_SUITE(utl_FactoryMap, FactoryMapFixture)

BOOST_AUTO_TEST_CASE(arity0_test1)
{
	std::auto_ptr<Animal> animal = animalFactory["Cat"].create();
	BOOST_CHECK_EQUAL(animal->type, CAT);
}

BOOST_AUTO_TEST_CASE(arity0_test2)
{
	std::auto_ptr<Animal> animal = animalFactory["Dog"].create();
	BOOST_CHECK_EQUAL(animal->type, DOG);
}

BOOST_AUTO_TEST_CASE(arity1_test1)
{
	std::auto_ptr<Monster> monster = monsterFactory["Hydra"].create(1200);
	BOOST_CHECK_EQUAL(monster->type, HYDRA);
	BOOST_CHECK_EQUAL(monster->hp, 1200);
}

BOOST_AUTO_TEST_CASE(arity1_test2)
{
	std::auto_ptr<Monster> monster = monsterFactory["Giant"].create(150);
	BOOST_CHECK_EQUAL(monster->type, GIANT);
	BOOST_CHECK_EQUAL(monster->hp, 150);
}

BOOST_AUTO_TEST_CASE(hasFactory_test1)
{
	BOOST_CHECK(!animalFactory.hasFactory("Bird"));
}

BOOST_AUTO_TEST_CASE(hasFactory_test2)
{
	BOOST_CHECK(monsterFactory.hasFactory("Troll"));
}

bool checkException1(const MissingFactory& exc)
{
	std::ostringstream msg;
	msg << "Cannot find a factory for type \'" << typeid(Monster).name() << "\' of kind \'Dragon\'";
	return exc.factoryKey == "Dragon" and exc.getMessage() == msg.str();
}

BOOST_AUTO_TEST_CASE(exception_test1)
{
	BOOST_CHECK_EXCEPTION(monsterFactory["Dragon"].create(1500), MissingFactory, checkException1);
}

BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(utl_FactoryRegistrar)

typedef Singleton<FactoryMap<std::string, Animal*()> > AnimalFactory;

#if BOOST_VERSION < 104300
REGISTER_FACTORY(AnimalFactory, "dog", createDog)
REGISTER_FACTORY(AnimalFactory, "cat", createCat)
#else
REGISTER_FACTORY(AnimalFactory, "dog", boost::factory<Dog*>())
REGISTER_FACTORY(AnimalFactory, "cat", boost::factory<Cat*>())
#endif

BOOST_AUTO_TEST_CASE(test1)
{
	std::auto_ptr<Animal> animal = AnimalFactory::instance()["dog"].create();
	BOOST_CHECK_EQUAL(animal->type, DOG);
}

BOOST_AUTO_TEST_CASE(test2)
{
	std::auto_ptr<Animal> animal = AnimalFactory::instance()["cat"].create();
	BOOST_CHECK_EQUAL(animal->type, CAT);
}

BOOST_AUTO_TEST_CASE(test3)
{
	BOOST_CHECK_THROW(AnimalFactory::instance()["horse"].create(), MissingFactory);
}

BOOST_AUTO_TEST_SUITE_END()
