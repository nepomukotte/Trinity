#include <iostream>
#include <string>
#include <fstream>
#include "PixelMap.h"
using namespace std;
 
/**
 * Function that uses a 2D array lookup table to return asadID, agetID, and channelID given sipmID.
 * @param **pixelMapArray the 2D array containing the lookup table.
 * @param sipmID the search parameter
 * @param *asadID one of the return parameters
 * @param *agetID one of the return parameters
 * @param *channelID one of the return parameters
**/
void pixelMapDecoder(int **pixelMapArray, uint16_t sipmID, uint16_t *asadID, uint16_t *agetID, uint16_t *channelID){
	int numSipm = 256;
	// Check if sipmID is out of bounds [0-511]
	if((sipmID<0)||(sipmID>numSipm)){
		cout << "sipmID is out of bounds! Must be in range [0-511]." << endl;
		*asadID = *agetID = *channelID = -1;
		return;	
	}

	// Perform simple linear search for desired sipmID
	for(int row=0; row<numSipm; row++){
		if(pixelMapArray[row][0] == sipmID){
			*asadID = pixelMapArray[row][1];
			*agetID = pixelMapArray[row][2];
			*channelID = pixelMapArray[row][3];
			return; 
		}
	}
}

/**
 * Function that uses a 2D array lookup table to return sipmID given asadID, agetID, and channelID. Completes the opposite
 * of pixelMapDecoder.
 * @param **pixelMapArray the 2D array containing the lookup table.
 * @param *sipmID the return parameter
 * @param asadID one of the search parameters
 * @param agetID one of the search parameters
 * @param channelID one of the search parameters
**/
void pixelMapEncoder(int **pixelMapArray, uint16_t *sipmID, uint16_t asadID, uint16_t agetID, uint16_t channelID){
	// Check if channelID is out of bounds. Bounds include [0-67] excluding 11,22,45,56
	if((channelID<0)||(channelID>67)||(channelID==11)||(channelID==22)||(channelID==45)||(channelID==56)){
		cout << "ChannelID is out of bounds! Must be in range 0-67 excluding 11, 22, 45, and 56." << endl;
		*sipmID = -1;
		return;
	}
	// Check if asadID is out of bounds [0-1]
	if((asadID<0)||(asadID>1)){
		cout << "AsadID is out of bounds! Must be in range 0-1." << endl;
		*sipmID = -1;
		return;
	}
	// Check if agetID is out of bounds [0-3]
	if((agetID<0)||(agetID>3)){
		cout << "AgetID is out of bounds! Must be in range 0-3." << endl;
		*sipmID = -1;
		return;
	}

	// determine row needed in from pixelMapArray. Subtractions account for certain channelIDs not being allowed.
	int row = 256*asadID + 64*agetID + channelID;
	if(channelID > 56){
		row -= 4;
	}else if(channelID > 45){
		row -= 3;
	}else if(channelID > 22){
		row -= 2;
	}else if(channelID > 11){
		row -= 1;
	}

	// set return SiPMID
	// 0 corresponds to the sipm collumn.
	*sipmID = pixelMapArray[row][0];
}

/**
 * Function which opens a csv file and allocates the data onto the heap as int values.
 * Assumming format of the csv is the following:
 * First row is a header row.
 * 512 rows and 4 columns (excluding header row).
 * Columns are in order sipmID, asadID, agetID, channelID.
 * @param pixelMapCSV string containing the file location/name of the csv file
 * @return a 2D integer array represented as a pointer to pointers to ints.
 */
int** openPixelMap(string pixelMapCSV){
	const int numRows = 256; //number of rows in csv (excluding header)
	const int numCols = 4; //number of collums in csv

	// create input file stream
	ifstream pixelMap(pixelMapCSV);

	// make sure file is open
	if(!pixelMap.is_open()) throw runtime_error("Can not open pixelMapCSV file.");

	// allocate memmory on the heap for storing pixelMap
	int **pixelMapArray = new int*[numRows];
	for(int k=0; k<numRows; ++k){
		pixelMapArray[k] = new int[numCols];
	}

	// buffer for storing line
	string line;

	// read header line
	getline(pixelMap, line);

	// assign csv values to allocated memmory
	for(int i=0; i<numRows; i++){
		//get next row
		getline(pixelMap, line);
		for(int j=0; j<numCols; j++){
			//store value
			pixelMapArray[i][j] = stoi(line.substr(0,line.find(",")));
			//get next value in row
			if(j<(numCols-1)){
				line = line.substr(line.find(",")+1);
			}
		}
	}

	return pixelMapArray;
	
}

/**
 * Function which deallocates the memmory stored on the heap created by openPixelMap.
 * @param **pixelMapArray variables whose data will be deallocated, previously return parameter from openPixelMap.
 */
void closePixelMap(int **pixelMapArray){
	const int numRows = 256; //number of rows in csv (excluding header)
	for(int i = 0; i < numRows; ++i) {
    	delete [] pixelMapArray[i];
	}
	delete [] pixelMapArray;
}

int myMain(int argc, char *argv[]){
	string fileName = argv[1];
	uint16_t sipmID, asadID, agetID, channelID;

	// // pixelMap memmory allocation test
	// int **pixelMapArray = openPixelMap(fileName);
	// for(int row=0; row<512; row++){
	// 	for(int col=0; col<4; col++){
	// 		cout << pixelMapArray[row][col] << ", ";
	// 	}
	// 	cout << endl;
	// }
	// closePixelMap(pixelMapArray);

	//pixelMapDecoder test
	int **pixelMapArray = openPixelMap(fileName);
	sipmID = 31;
	pixelMapDecoder(pixelMapArray, sipmID, &asadID, &agetID, &channelID);
	closePixelMap(pixelMapArray);

	// // pixelMapEncoder test
	// int **pixelMapArray = openPixelMap(fileName);
	// asadID = 0;
	// agetID = 0;
	// channelID = 0;
	// pixelMapEncoder(pixelMapArray, &sipmID, asadID, agetID, channelID);
	// closePixelMap(pixelMapArray);

	cout << "sipmID: " << sipmID << endl;
	cout << "asadID: " << asadID << endl;
	cout << "agetID: " << agetID << endl;
	cout << "chanelID: " << channelID << endl;

	return 0;
}
