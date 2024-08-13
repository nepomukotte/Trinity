#include "EventBuilder.h"

int main(int argc, char* argv[]){
	
	
	std::string coboFormats = argv[1];
//	cout<<coboFormats<<endl;
	std::string filename_0 = argv[2];
	//std::string filename_1 = argv[3];
	std::string filename_2 = argv[3];
//	cout<<"H"<<endl;
	std::string outDir = argv[4];
//	cout<<"E"<<endl;
	EventBuilder eB(coboFormats);
//	cout<<"R"<<endl;
	eB.mainFlow(filename_0,  filename_2, outDir);

	return 0;
}
