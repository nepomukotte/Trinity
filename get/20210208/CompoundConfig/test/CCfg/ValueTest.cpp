// =============================================================================
// ValueTest.cpp
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Frederic Chateau (frederic.chateau@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "CCfg/Value.h"
#include <boost/test/auto_unit_test.hpp>
#include <sstream>
#include <cstring>

namespace CCfg
{

///////////////////////////////////////////////////////////////////////////////
// test fixture definition
///////////////////////////////////////////////////////////////////////////////

struct config_value_data
{
	Value x;
	Value message;
	Value isFoo;
	Value fooCount;
	Value fooAddress;

	config_value_data()
		: x(2.5),
		message("Hello World!"),
		isFoo(true),
		fooCount(42),
		fooAddress(0xDEADBEEF)
	{
	}

	~config_value_data()
	{
	}
};

BOOST_FIXTURE_TEST_SUITE(CCfg_Value, config_value_data)

///////////////////////////////////////////////////////////////////////////////
// test methods body
///////////////////////////////////////////////////////////////////////////////


/**
 * Checks getters.
 */
BOOST_AUTO_TEST_CASE(testValueGetters)
{
	BOOST_CHECK_EQUAL(x.getReal(), 2.5);
	BOOST_CHECK_EQUAL(message.getString(), "Hello World!");
	BOOST_CHECK_EQUAL(isFoo.getBool(), true);
	BOOST_CHECK_EQUAL(fooCount.getInt(), 42);
	BOOST_CHECK_EQUAL(fooAddress.getHex(), 0xDEADBEEF);
}

/**
 * Checks conversion getters for integers.
 */
BOOST_AUTO_TEST_CASE(testIntConversion)
{
	BOOST_CHECK_MESSAGE(fooCount.getAsString() == "42", "int => string");
	BOOST_CHECK_MESSAGE(fooCount.getAsReal() == 42.0, "int => real");
	BOOST_CHECK_MESSAGE(fooCount.getAsInt() == 42, "int => int");
	BOOST_CHECK_MESSAGE(fooCount.getAsBool() == true, "int => bool");
	BOOST_CHECK_MESSAGE(fooCount.getAsHex() == 42u, "int => hex");
}

/**
 * Checks conversion getters for hexadecimals.
 */
BOOST_AUTO_TEST_CASE(testHexConversion)
{
	BOOST_CHECK_MESSAGE(fooAddress.getAsString() == "0xDEADBEEF", "hex => string");
	BOOST_CHECK_MESSAGE(fooAddress.getAsReal() == 3735928559.0, "hex => real");
	BOOST_CHECK_MESSAGE(fooAddress.getAsInt() == -559038737, "hex => int");
	BOOST_CHECK_MESSAGE(fooAddress.getAsBool() == true, "hex => bool");
	BOOST_CHECK_MESSAGE(fooAddress.getAsHex() == 0xDEADBEEF, "hex => hex");
}

/**
 * Checks conversion getters for real numbers.
 */
BOOST_AUTO_TEST_CASE(testRealConversion)
{
	BOOST_CHECK_MESSAGE(x.getAsString() == "2.5", "real => string");
	BOOST_CHECK_MESSAGE(x.getAsReal() == 2.5, "real => real");
	BOOST_CHECK_MESSAGE(x.getAsInt() == 2, "real => int");
	BOOST_CHECK_MESSAGE(x.getAsBool() == true, "real => bool");
	BOOST_CHECK_MESSAGE(x.getAsHex() == 2u, "real => hex");

	x = 5.0;
	BOOST_CHECK_MESSAGE(x.getAsString() == "5.0", "real => string");
	BOOST_CHECK_MESSAGE(x.getAsReal() == 5.0, "real => real");
	BOOST_CHECK_MESSAGE(x.getAsInt() == 5, "real => int");
	BOOST_CHECK_MESSAGE(x.getAsBool() == true, "real => bool");
	BOOST_CHECK_MESSAGE(x.getAsHex() == 5u, "real => hex");

	x = 123456789012345.0;
	BOOST_CHECK_MESSAGE(x.getAsString() == "123456789012345.0", "real => string");

	x = -1234567890123456.0;
	BOOST_CHECK_MESSAGE(x.getAsString() == "-1.23456789012346e+15", "real => string");

	x = 1.23456e-4;
	BOOST_CHECK_MESSAGE(x.getAsString() == "0.000123456", "real => string");

	x = -9.9999e-5;
	BOOST_CHECK_MESSAGE(x.getAsString() == "-9.9999e-05", "real => string");

	x = 1234567890.00001;
	BOOST_CHECK_MESSAGE(x.getAsString() == "1234567890.00001", "real => string");

	x = 1234567890.000001;
	BOOST_CHECK_MESSAGE(x.getAsString() == "1234567890.0", "real => string");

	x = 0.0;
	BOOST_CHECK_MESSAGE(x.getAsString() == "0.0", "real => string");

	x = 1.2345678912345e-17;
	BOOST_CHECK_MESSAGE(x.getAsString() == "1.2345678912345e-17", "real => string");
}

/**
 * Checks conversion getters for booleans.
 */
BOOST_AUTO_TEST_CASE(testBooleanConversion)
{
	BOOST_CHECK_MESSAGE(isFoo.getAsString() == "true", "bool => string");
	BOOST_CHECK_MESSAGE(isFoo.getAsReal() == 1.0, "bool => real");
	BOOST_CHECK_MESSAGE(isFoo.getAsInt() == 1, "bool => int");
	BOOST_CHECK_MESSAGE(isFoo.getAsBool() == true, "bool => bool");
	BOOST_CHECK_MESSAGE(isFoo.getAsHex() == 1u, "bool => hex");
}



/**
 * Checks conversion getters for integers.
 */
BOOST_AUTO_TEST_CASE(testIntConversion2)
{
	BOOST_CHECK_MESSAGE(fooCount.getAsString() == "42", "int => string");
	BOOST_CHECK_MESSAGE(fooCount.getAsReal() == 42.0, "int => real");
	BOOST_CHECK_MESSAGE(fooCount.getAsInt() == 42, "int => int");
	BOOST_CHECK_MESSAGE(fooCount.getAsBool() == true, "int => bool");
	BOOST_CHECK_MESSAGE(fooCount.getAsHex() == 42u, "int => hex");
}

/**
 * Checks conversion getters for strings.
 */
BOOST_AUTO_TEST_CASE(testStringConversion)
{
	#define CATCH_BAD_CONV(desc, expr) \
		try	\
		{	\
			expr;	\
			BOOST_FAIL(desc " should throw BadConversion");\
		} catch(Value::BadConversion& e) \
		{ \
		}

	// string as string => other
	BOOST_CHECK_MESSAGE(message.getAsString() == "Hello World!", "string as string => string");
	CATCH_BAD_CONV("string as string => real", message.getAsReal())
	CATCH_BAD_CONV("string as string => int", message.getAsInt());
	CATCH_BAD_CONV("string as string => bool", message.getAsBool());
	CATCH_BAD_CONV("string as string => hex", message.getAsHex());

	// real as string => other
	message = "2.1";
	BOOST_CHECK_MESSAGE(message.getAsString() == "2.1", "real as string => string");
	BOOST_CHECK_MESSAGE(message.getAsReal() == 2.1, "real as string => real");
	CATCH_BAD_CONV("real as string => int", message.getAsInt());
	CATCH_BAD_CONV("real as string => bool", message.getAsBool());
	CATCH_BAD_CONV("real as string => hex", message.getAsHex());

	// int as string => other
	message = "5";
	BOOST_CHECK_MESSAGE(message.getAsString() == "5", "int as string => string");
	BOOST_CHECK_MESSAGE(message.getAsReal() == 5.0, "int as string => real");
	BOOST_CHECK_MESSAGE(message.getAsInt() == 5, "int as string => int");
	CATCH_BAD_CONV("int as string => bool", message.getAsBool());
	BOOST_CHECK_MESSAGE(message.getAsHex() == 5u, "int as string => hex");

	// boolean as string => other
	message = "true";
	BOOST_CHECK_MESSAGE(message.getAsString() == "true", "bool as string => string");
	CATCH_BAD_CONV("bool as string => real", message.getAsReal());
	CATCH_BAD_CONV("bool as string => int", message.getAsInt());
	BOOST_CHECK_MESSAGE(message.getAsBool() == true, "bool as string => bool");
	CATCH_BAD_CONV("bool as string => hex", message.getAsHex());
	message = "false";
	BOOST_CHECK_MESSAGE(message.getAsString() == "false", "bool as string => string");
	CATCH_BAD_CONV("bool as string => real", message.getAsReal());
	CATCH_BAD_CONV("bool as string => int", message.getAsInt());
	BOOST_CHECK_MESSAGE(message.getAsBool() == false, "bool as string => bool");
	CATCH_BAD_CONV("bool as string => hex", message.getAsHex());

	// hex as string => other
	message = "0xBAADF00D";
	BOOST_CHECK_MESSAGE(message.getAsString() == "0xBAADF00D", "hex as string => string");
	// Converting "0xBAADF00D" into a floating-point hexadecimal is valid in C++11 which uses std::strtold instead of std::scanf.
	//CATCH_BAD_CONV("hex as string => real", message.getAsReal());
	CATCH_BAD_CONV("hex as string => int", message.getAsInt());
	CATCH_BAD_CONV("hex as string => bool", message.getAsBool());
	BOOST_CHECK_MESSAGE(message.getAsHex() == 0xBAADF00D, "hex as string => hex");
}

/**
 * Checks automatic conversion.
 */
BOOST_AUTO_TEST_CASE(testAutoConversion)
{
	double xVal = x;
	BOOST_CHECK_EQUAL(xVal, 2.5);

	const char* messageVal = message;
	BOOST_CHECK(not strcmp(messageVal, "Hello World!"));

	std::string messageVal2 = message;
	BOOST_CHECK_EQUAL(messageVal2, "Hello World!");

	bool isFooVal = isFoo;
	BOOST_CHECK_EQUAL(isFooVal, true);

	int fooCountVal = fooCount;
	BOOST_CHECK_EQUAL(fooCountVal, 42);

	unsigned fooAddressVal = fooAddress;
	BOOST_CHECK_EQUAL(fooAddressVal, 0xDEADBEEF);
}

/**
 * Checks cross-type assignment.
 */
BOOST_AUTO_TEST_CASE(testValueAssignment)
{
	x = 3.1;
	double realVal = x;
	BOOST_CHECK_EQUAL(realVal, 3.1);

	x = "Hello";
	const std::string& strVal = x;
	BOOST_CHECK_EQUAL(strVal, "Hello");

	message = true;
	bool boolVal = message;
	BOOST_CHECK_EQUAL(boolVal, true);

	isFoo = 99;
	int intVal = isFoo;
	BOOST_CHECK_EQUAL(intVal, 99);

	message = 0xBAADF00D;
	unsigned hexVal = message;
	BOOST_CHECK_EQUAL(hexVal, 0xBAADF00D);

}

/**
 * Checks type auto determination
 */
BOOST_AUTO_TEST_CASE(testIntTypeInference)
{
	x.setValue("2");
	int xVal = x;
	BOOST_CHECK_EQUAL(xVal, 2);
}

/**
 * Checks type auto determination
 */
BOOST_AUTO_TEST_CASE(testRealTypeInference)
{
	x.setValue("-5.5");
	double xVal = x;
	BOOST_CHECK_EQUAL(xVal, -5.5);
}

/**
 * Checks type auto determination
 */
BOOST_AUTO_TEST_CASE(testBoolTypeInference)
{
	x.setValue("false");
	bool xVal = x;
	BOOST_CHECK_EQUAL(xVal, false);
}

/**
 * Checks type auto determination
 */
BOOST_AUTO_TEST_CASE(testStringTypeInference)
{
	x.setValue("Hello");
	std::string xVal = x;
	BOOST_CHECK_EQUAL(xVal, "Hello");

	x.setValue("\"hi!\"");
	BOOST_CHECK_EQUAL(x.getString(), "\"hi!\"");
}

/**
 * Checks type auto determination
 */
BOOST_AUTO_TEST_CASE(testHexTypeInference)
{
	x.setValue("0xACE0FCEA");
	unsigned xVal = x;
	BOOST_CHECK_EQUAL(xVal, 0xACE0FCEA);
}

/**
 * Checks BadConversion exception.
 */
BOOST_AUTO_TEST_CASE(testBadConversion)
{
	try
	{
		message = "2.1";
		message.getAsHex();
		BOOST_FAIL("Should throw BadConversion");
	} catch(Value::BadConversion& e)
	{
		BOOST_CHECK_EQUAL(e.getMessage(),
			"Bad conversion: the string value \"2.1\" "
			"cannot be converted to hexadecimal.");
	}
}

BOOST_AUTO_TEST_CASE(testStream)
{
	std::ostringstream str;
	str << x;
	BOOST_CHECK_EQUAL(str.str(), "2.5");

	str.str(std::string());
	str << fooAddress;
	BOOST_CHECK_EQUAL(str.str(), "0xDEADBEEF");
}

BOOST_AUTO_TEST_SUITE_END()

}
