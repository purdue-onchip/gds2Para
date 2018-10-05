// Example program to test Geometry class
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "interface.h"
using namespace std;

int main(int argc, char **argv)
{
	// Analyze executable inputs
	string gdtFile;
	string designName;
	size_t indTopCell;
	unsigned long plotPixels;
	vector<double> layerBounds;
	if (argc == 1)
	{
		// Use example files
		gdtFile = "./nand2.gdt"; // string gdtFile = "./SDFFRS_X2.gdt"; // string gdtFile = "./4004.gdt";
	    designName = "nand2"; // string designName = "sdffrs"; // string designName = "4004";
		indTopCell = 3; // size_t indTopCell = 0; // size_t indTopCell = 0;
		//string topCellName = "abc2"; // string topCellName = "SDFFRS_X2"; // string topCellName = "4004";
		plotPixels = 400; // unsigned long plotPixels = 400; // unsigned long plotPixels = 8000;
		layerBounds = { -20, 180, -60, 140 }; // vector<double> layerBounds = { -2, 14, -6, 10 } // vector<double> layerBounds = { 0, 2800, 0, 2800 };
    }
	else if (argc == 2)
	{
		// Use command line input
		gdtFile = argv[1];
		designName = gdtFile.substr(0,gdtFile.find_last_of("."));
	}
	else
	{
		// Unsupported number of arguments
		cerr << "Must use GDT file name as argument or no arguments" << endl;
	}

    // Load GDT file
    Geometry geometry;
	geometry.load(gdtFile);
    geometry.print();
	//(geometry.getCell(0)).print();
	if (argc == 2)
	{
		// Assume the final cell is the top cell
		//vector<string> names = geometry.findNames(); // Analyze list of cell name otherwise
		indTopCell = (size_t) (geometry.getNumCell() - 1);
		plotPixels = 400; // Must generalize
		layerBounds = { -20, 180, -60, 140 }; // Must generalize
	}
	GeoCell topCell = geometry.getCell(indTopCell);
	string topCellName = topCell.getCellName();
    topCell.print();
    cout << "End of Testing" << endl;

	// Write file to gnuplot stored information
	cout << "Starting to Write Plot File" << endl;
	bool plotFlag = geometry.writePlotEnv(designName, topCellName, plotPixels, layerBounds);
	if (plotFlag)
	{
        cout << "Finished Writing Plot File" << endl;
	}
	else
	{
		cerr << "Could not write to location" << endl;
	}

    // Clean up (to implement later)
}
