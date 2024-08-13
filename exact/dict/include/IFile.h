#ifndef IFILE_H
#define IFILE_H

#include <iostream>
#include <TFile.h>
#include <TTree.h>

/*

Class used to contain TFile with events.

It allows loading the trees and using the data in methods and modules without losing reference to the TFile. 

Since ROOT Histograms attach to the TFile, creating it inside a method destroys it once the method exits. This causes SegFaults. In order to avoid this, the class IFile was created that contains all the trees in the raw data file from the CT.

Constructor:
IFile(std::string filename)
Creates a new TFile and keeps it alive, while initializing all trees from a standard DAQ from the CT.

A Tree ID was created to keep reference of the different trees

0: HLED
1; BiFocal
2: Test
3: Forced

*/

class IFile{
	public:
		TFile *ITFile;
		TTree *treeHLED;
		TTree *treeBiFocal;
		TTree *treeForced;
		TTree *treeTest;

		IFile(std::string filename, std::string mode="READ");
		~IFile();

		/*
		Get Tree ID number

		Arguments:
		std::string treeName: Not case sensitive, can be hled, bifocal, test, forced

		returns:
		int specifying the ID

		*/
		static int GetTreeID(std::string treeName);
	private:
		void GetTrees();
};
#endif