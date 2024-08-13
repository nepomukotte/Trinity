#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <TGraph.h>

using namespace std;

uint16_t flipfunc(uint8_t lowbyte, uint8_t highbyte)
{
  uint16_t value;
  value = ((highbyte << 8) | lowbyte);
  return value;
}

int temp_sipm () {

  vector<uint16_t> linuxtime;
  vector<vector<uint16_t>> TempofSiPM;

  int NofSIAB = 16;
  int NofFiles = 50;
  int NofSamples = 60;

  for (int i=0; i<NofFiles; i++) {
    char * buffer = NULL;
    FILE * pFile = NULL;
    long lSize;
    size_t result;
    char filename[40];
    snprintf(filename, sizeof(filename), "TVAC_SiPM/file_SiPM_%d.txt", i);
    pFile = fopen (filename,"rb");
    fseek (pFile , 0 , SEEK_END);
    lSize = ftell (pFile);
    rewind (pFile);

    buffer = (char*) malloc (sizeof(char)*lSize);
    result = fread (buffer,1,lSize,pFile);
  
    for(int j=0; j<NofSamples; j++){
      vector<uint16_t> newsample;
      linuxtime.push_back(flipfunc(*(buffer+34*j), *(buffer+34*j+1)));
      for(int k=0; k<16; k++){
        if(k == 3) {
          newsample.push_back(20);
        } else {
          newsample.push_back(flipfunc(*(buffer+34*j+2*k+2), *(buffer+34*j+2*k+3)));
        }
      }
      TempofSiPM.push_back(newsample);
    }
  
  fclose (pFile);
  free (buffer);
  }

  int TimeofTest[NofFiles*NofSamples];
  int Temperature[NofSIAB][NofFiles*NofSamples];
  
  for(int i=0; i<NofSIAB; i++){
    for(int j=0; j<NofFiles*NofSamples; j++){
      if(TempofSiPM[j][i] > 2000){
        Temperature[i][j]=(TempofSiPM[j-1][i] + TempofSiPM[j+1][i])/2;
      }
      else{
        Temperature[i][j]=TempofSiPM[j][i];
      }
      TimeofTest[j]=linuxtime[j]-28000;
    }
  }

  for(int i=0; i<NofFiles*NofSamples; i++){
    cout << (TimeofTest[i] / 60) << "\t";
    for(int j=0; j<NofSIAB; j++){
      double value = ((((Temperature[j][i] * 3.3)/1024) - 0.389) / 0.0443);
      if(value < 25) {
        cout << ((((Temperature[j][i-1] * 3.3)/1024) - 0.389) / 0.0443) << "\t";
      }
      if(value > 25) {
        cout << value << "\t";
      }
    }
    cout << endl;
  }
/*
  TGraph *g[NofSIAB];
  TMultiGraph *mg = new TMultiGraph();
  for (int i=0; i<NofSIAB; i++) {
    g[i] = new TGraph(NofFiles*NofSamples, TimeofTest, Temperature[i]);
    g[i]->SetMarkerColor(kBlue+i);
    mg->Add(g[i]);
  }
  mg->Draw("A pmc plc");
*/
  return 0;
}