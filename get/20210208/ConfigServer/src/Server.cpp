// =============================================================================
// Server.cpp
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
// Developed at the CEA (Saclay, France), in the Irfu/Sedi/Lilas laboratory, by:
// Loïc HERVE (loic.herve@cea.fr)
// -----------------------------------------------------------------------------
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation.
// =============================================================================

#include "CCfg/ConfigServerI.h"
#include ICE_HEADER(ObjectAdapter.h)
#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Exception.h)
#include <iostream>

using namespace std;
using namespace CCfg;

///////////////////////////////////////////////////////////////////////////////
// Main function of the server
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char* argv[])
{
	int status = 0;
	Ice::CommunicatorPtr ic;

	try {
		ic = Ice::initialize(argc, argv);
		Ice::ObjectAdapterPtr adapter = ic->createObjectAdapterWithEndpoints("ConfigServerAdapter", "default -p 10000"); //Creation of the ConfigServer connection
		Ice::ObjectPtr object = new ConfigServerI;
		adapter->add(object, ic->stringToIdentity("ConfigServer"));
		adapter->activate();
		ic->waitForShutdown();
    } catch (const Ice::Exception& e) {
        cerr << e.what() << endl;
        status = 1;
    } catch (const char* msg) {
        cerr << msg << endl;
        status = 1;
    }
    if (ic) {
        try {
	ic->destroy();
        } catch (const Ice::Exception& e) {
            cerr << e.what() << endl;
            status = 1;
        }
    }
    return status;
}
