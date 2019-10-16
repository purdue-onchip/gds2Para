#ifndef GDS2PARA_MAINTESTFDTDMESH_H_
#define GDS2PARA_MAINTESTFDTDMESH_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

inline void WriteVectorIn1Line(ofstream &file_obj, vector<double> vec_x) {
	for (const auto& j : vec_x) {
		file_obj << j << ' ';
	}
	file_obj << endl;
}

// Write object sys to file
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