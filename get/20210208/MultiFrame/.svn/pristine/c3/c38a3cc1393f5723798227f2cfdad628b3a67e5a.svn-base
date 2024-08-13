// =============================================================================
// main.cpp
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

#include "tut.h"
#include "tut_reporter.h"

namespace tut
{
	test_runner_singleton runner;
}

using tut::runner;
using tut::reporter;

int main()
{
	// Instantiate and set the callback that displays test run reports
	reporter callback; // use std::cout by default to output report.
	runner.get().set_callback(&callback);

	// run all tests in all groups
	//runner.get().run_tests();
	runner.get().run_tests("PrimaryHeader");
	runner.get().run_tests("BlobHeader");
	runner.get().run_tests("StandardHeader");
	runner.get().run_tests("Frame");
	runner.get().run_tests("Field");
	runner.get().run_tests("BitField");
	runner.get().run_tests("FrameFormat");
	runner.get().run_tests("FrameDictionary");
	runner.get().run_tests("Signed");
	runner.get().run_tests("BlobFrame");
	runner.get().run_tests("LayeredFrame");
	runner.get().run_tests("FrameBuilder");
	//runner.get().run_test("Frame", 2);

	return 0;
}
