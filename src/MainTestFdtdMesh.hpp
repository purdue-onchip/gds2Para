#ifndef GDS2PARA_MAINTESTFDTDMESH_H_
#define GDS2PARA_MAINTESTFDTDMESH_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

inline void WriteVectorIn1Line(ofstream &file_obj, const vector<double> &vec_x) {
	for (const auto& j : vec_x) {
		file_obj << j << ' ';
	}
	file_obj << endl;
}

// Write object sys to files
void WriteSysToFile(const fdtdMesh &sys) {
	ofstream file_obj;

	// Write primitive types like a bool, int, or float in one txt file
	file_obj.open("sys_PrimitiveType.txt", ios::out);
	file_obj << sys.numCdtRow << endl;
	file_obj << sys.lengthUnit << endl;
	file_obj << sys.freqUnit << endl;
	file_obj << sys.nfreq << endl;
	file_obj << sys.freqStart << endl;
	file_obj << sys.freqEnd << endl;
	file_obj << sys.freqScale << endl;
	file_obj << sys.xlim1 << endl;
	file_obj << sys.xlim2 << endl;
	file_obj << sys.ylim1 << endl;
	file_obj << sys.ylim2 << endl;
	file_obj << sys.zlim1 << endl;
	file_obj << sys.zlim2 << endl;
	file_obj << sys.numStack << endl;
	file_obj << sys.numPorts << endl;

	file_obj.close();

	// Write each vector in independent file
	file_obj.open("sys_vec_stackEps.txt", ios::out);
	for (const auto& i : sys.stackEps)
		file_obj << i << endl;
	file_obj.close();

	file_obj.open("sys_vec_stackSig.txt", ios::out);
	for (const auto& i : sys.stackSig)
		file_obj << i << endl;
	file_obj.close();

	file_obj.open("sys_vec_stackBegCoor.txt", ios::out);
	for (const auto& i : sys.stackBegCoor)
		file_obj << i << endl;
	file_obj.close();

	file_obj.open("sys_vec_stackEndCoor.txt", ios::out);
	for (const auto& i : sys.stackEndCoor)
		file_obj << i << endl;
	file_obj.close();

	file_obj.open("sys_vec_stackName.txt", ios::out);
	for (const auto& i : sys.stackName)
		file_obj << i << endl;
	file_obj.close();

	// Write sys.portCoor, only export what is used by function meshAndMark
	// output format: multiplicity ~ 1st line; x1 ~ 2nd line; x2, y1, y2, z1, z2 ~ 3-7 line.
	// Every object of class fdtdPort takes 7 lines, being one element of vector portCoor
	file_obj.open("sys_vec_portCoor.txt", ios::out);
	for (int i = 0; i < sys.portCoor.size(); i++) {
		const fdtdPort &ofdtdPort = sys.portCoor[i];

		file_obj << ofdtdPort.multiplicity << endl;      // line 1 ~ multiplicity
		WriteVectorIn1Line(file_obj, ofdtdPort.x1);      // line 2 ~ x1
		WriteVectorIn1Line(file_obj, ofdtdPort.x2);      // line 3 ~ x2
		WriteVectorIn1Line(file_obj, ofdtdPort.y1);      // line 4 ~ y1
		WriteVectorIn1Line(file_obj, ofdtdPort.y2);      // line 5 ~ y2
		WriteVectorIn1Line(file_obj, ofdtdPort.z1);      // line 6 ~ z1
		WriteVectorIn1Line(file_obj, ofdtdPort.z2);      // line 7 ~ z2
	}
	file_obj.close();

	// Write sys.conductorIn, only export what is used by function meshAndMark
	file_obj.open("sys_vec_conductorIn.txt", ios::out);
	for (int i = 0; i < sys.conductorIn.size(); i++) {
		const fdtdOneCondct &ofdtdOneCondct = sys.conductorIn[i];

		file_obj << ofdtdOneCondct.numVert << ' ' << ofdtdOneCondct.zmin << ' ' 
			<< ofdtdOneCondct.zmax << ' ' 
			<< ofdtdOneCondct.layer << endl;             // line 1 ~ numVert zmin zmax layer
		for (int j = 0; j < ofdtdOneCondct.numVert; j++) {
			file_obj << ofdtdOneCondct.x[j] << ' ';
		}
		file_obj << endl;                                // line 2 ~ x
		for (int j = 0; j < ofdtdOneCondct.numVert; j++) {
			file_obj << ofdtdOneCondct.y[j] << ' ';
		}
		file_obj << endl;                                // line 3 ~ y
	}
	file_obj.close();
	

	// Write *pointer to file 
	file_obj.open("sys_ptr_stackCdtMark.txt", ios::out);
	for (int i = 0; i < sys.numStack; i++)
		file_obj << sys.stackCdtMark[i] << endl;
	file_obj.close();


}

void Read1LineToVector(ifstream &file_obj, vector<double> *pVec) {
	double doubVal;
	string line;

	getline(file_obj, line);
	istringstream myStream(line);

	while (myStream >> doubVal) {
		(*pVec).push_back(doubVal);
	}
}

double* Read1LineToDoublePtr(ifstream &file_obj) {
	string line;
	int num_double = 0;
	double doubVal;

	getline(file_obj, line);

	istringstream myStream(line);
	
	while (myStream >> doubVal) {
		num_double++;
	}	// Get num of double numbers in the 1 text line

	myStream.clear();
	myStream.seekg(0, ios::beg);
	double *ptr = (double *)malloc(num_double * sizeof(double));
	for (int i = 0; i < num_double; i++) {
		myStream >> ptr[i];
	}

	return ptr;
}

// Read object sys from files
void ReadSysFromFile(fdtdMesh *psys) {
	double doubVecValue;
	int num_line = 0;
	string strVecValue, line;

	ifstream file_obj;

	// Read primitive types like a bool, int, or float
	file_obj.open("sys_PrimitiveType.txt", ios::in);
	if (!file_obj.is_open()) {
		perror("Error open");
		exit(EXIT_FAILURE);
	}
	file_obj >> psys->numCdtRow;
	file_obj >> psys->lengthUnit;
	file_obj >> psys->freqUnit;
	file_obj >> psys->nfreq;
	file_obj >> psys->freqStart;
	file_obj >> psys->freqEnd;
	file_obj >> psys->freqScale;
	file_obj >> psys->xlim1;
	file_obj >> psys->xlim2;
	file_obj >> psys->ylim1;
	file_obj >> psys->ylim2;
	file_obj >> psys->zlim1;
	file_obj >> psys->zlim2;
	file_obj >> psys->numStack;
	file_obj >> psys->numPorts;

	file_obj.close();

	// Read each vector from independent file
	file_obj.open("sys_vec_stackEps.txt", ios::in);
	while (file_obj >> doubVecValue) {
		psys->stackEps.push_back(doubVecValue);
	}
	file_obj.close();

	file_obj.open("sys_vec_stackSig.txt", ios::in);
	while (file_obj >> doubVecValue) {
		psys->stackSig.push_back(doubVecValue);
	}
	file_obj.close();

	file_obj.open("sys_vec_stackBegCoor.txt", ios::in);
	while (file_obj >> doubVecValue) {
		psys->stackBegCoor.push_back(doubVecValue);
	}
	file_obj.close();

	file_obj.open("sys_vec_stackEndCoor.txt", ios::in);
	while (file_obj >> doubVecValue) {
		psys->stackEndCoor.push_back(doubVecValue);
	}
	file_obj.close();

	file_obj.open("sys_vec_stackName.txt", ios::in);
	while (file_obj >> strVecValue) {
		psys->stackName.push_back(strVecValue);
	}
	file_obj.close();

	// Read sys.portCoor, only input what is used by function meshAndMark
	// Input format: multiplicity ~ 1st line; x1 ~ 2nd line; x2, y1, y2, z1, z2 ~ 3-7 line.
	// Every object of class fdtdPort takes 7 lines, being one element of vector portCoor
	file_obj.open("sys_vec_portCoor.txt", ios::in);                                     
	
	num_line = 0;
	while (getline(file_obj, line)) {
		num_line++;
	}	// Get num of lines in the file, should be num_line = 7*sys.portCoor.size()
	
	file_obj.clear();
	file_obj.seekg(0, ios::beg);
	for (int i = 0; i < num_line / 7; i++) {
		fdtdPort ofdtdPort;

		file_obj >> ofdtdPort.multiplicity;             // line 1 ~ multiplicity
		getline(file_obj, line);
		Read1LineToVector(file_obj, &ofdtdPort.x1);      // line 2 ~ x1
		Read1LineToVector(file_obj, &ofdtdPort.x2);      // line 3 ~ x2
		Read1LineToVector(file_obj, &ofdtdPort.y1);      // line 4 ~ y1
		Read1LineToVector(file_obj, &ofdtdPort.y2);      // line 5 ~ y2
		Read1LineToVector(file_obj, &ofdtdPort.z1);      // line 6 ~ z1
		Read1LineToVector(file_obj, &ofdtdPort.z2);      // line 7 ~ z2

		psys->portCoor.push_back(ofdtdPort);
	}
	file_obj.close();

	// Read sys.conductorIn, only input what is used by function meshAndMark
	file_obj.open("sys_vec_conductorIn.txt", ios::in);

	num_line = 0;
	while (getline(file_obj, line)) {
		num_line++;
	}	// Get num of lines in the file, should be num_line = 3*sys.conductorIn.size()

	file_obj.clear();
	file_obj.seekg(0, ios::beg);
	for (int i = 0; i < num_line / 3; i++) {
		fdtdOneCondct ofdtdOneCondct;

		file_obj >> ofdtdOneCondct.numVert >> ofdtdOneCondct.zmin 
			>> ofdtdOneCondct.zmax 
			>> ofdtdOneCondct.layer;                        // line 1 ~ numVert zmin zmax layer
		getline(file_obj, line);
		ofdtdOneCondct.x = Read1LineToDoublePtr(file_obj);   // line 2 ~ x
		ofdtdOneCondct.y = Read1LineToDoublePtr(file_obj);   // line 3 ~ y

		psys->conductorIn.push_back(ofdtdOneCondct);
	}
	file_obj.close();


	// Read *pointer from file 
	file_obj.open("sys_ptr_stackCdtMark.txt", ios::in);
	psys->stackCdtMark = (double *)malloc(psys->numStack * sizeof(double));
	for (int i = 0; i < psys->numStack; i++) {
		file_obj >> psys->stackCdtMark[i];
	}
	file_obj.close();

}

#ifndef SKIP_COO2CSR_MALLOC
#define SKIP_COO2CSR_MALLOC
int COO2CSR_malloc(myint *rowId, myint *ColId, double *val, myint totalnum, myint leng, myint *rowId1) {    // totalnum is the total number of entries, leng is the row number
	int i;
	int *rowId2;
	int count, start;
	int k;

	rowId2 = (int*)malloc((leng + 1) * sizeof(int));
	count = 0;
	i = 0;
	k = 0;
	rowId2[k] = 0;
	k++;
	while (i < totalnum) {
		start = rowId[i];
		while (i < totalnum && rowId[i] == start) {
			count++;
			i++;
		}
		rowId2[k] = (count);
		k++;
	}

	for (i = 0; i <= leng; i++) {
		rowId1[i] = rowId2[i];
	}

	free(rowId2); rowId2 = NULL;
	return 0;
}
#endif


#endif  // GDS2PARA_MAINTESTFDTDMESH_H_