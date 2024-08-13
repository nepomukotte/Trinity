#ifndef PIXELMPAP_H
#define PIXELMAP_H

#include<iostream>
using namespace std;

void pixelMapDecoder(int **pixelMapArray, uint16_t sipmID, uint16_t *asadID, uint16_t *agetID, uint16_t *channelID);
void pixelMapEncoder(int **pixelMapArray, uint16_t *sipmID, uint16_t asadID, uint16_t agetID, uint16_t channelID);
int** openPixelMap(string pixelMapCSV);
void closePixelMap(int **pixelMapArray);
int myMain(int argc, char *argv[]);

#endif