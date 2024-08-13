#include "ConfigServer.h"
#include ICE_HEADER(Initialize.h)
#include ICE_HEADER(Communicator.h)
#include ICE_HEADER(Exception.h)
#include <iostream>

using namespace std;
using namespace CCfg;


///////////////////////////////////////////////////////////////////////////////
// Main function of the client
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    int status = 0;

    Ice::CommunicatorPtr ic;
    try {
      ic = Ice::initialize(argc, argv);
      Ice::ObjectPrx base = ic->stringToProxy("ConfigServer:default -p 10000");
      ConfigServerPrx serv = ConfigServerPrx::checkedCast(base);
      if (!serv)
            throw "Invalid proxy";
      //HERE GO THE CODE

      //END OF CODE
	} catch (CCfg::ExceptionForIce& exIce) { //Display The CCfg Exceptions through Ice
    	cerr << "ExceptionForIce : " << exIce.reason << endl;
    	status = 1;
    }catch (const Ice::Exception& ex) {
        cerr << ex.what() << endl;
        status = 1;
   } catch (const char* msg) {
        cerr << msg << endl;
        status = 1;
    }
    if (ic)
      ic->destroy();
    return status;
}
