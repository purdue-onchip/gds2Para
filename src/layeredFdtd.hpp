#ifndef GDS2PARA_LAYERED_FDTD_H_
#define GDS2PARA_LAYERED_FDTD_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"
#include "sysInfoIO.hpp"              // IO object "sys"
#include "pardisoSolver.hpp"          // Solve E and Z-para in pardiso

using namespace std;

int layeredFdtd(void) {

	fdtdMesh sys;

	// Read object sys (class fdtdMesh) from files
	ReadSysFromFile(&sys);

	// Write object sys to files
	WriteSysToFile(sys);

	// Mesh the domain and mark conductors
	unordered_set<double> portCoorx, portCoory;
	unordered_map<double, int> xi, yi, zi;
	clock_t t2 = clock();
	int status = meshAndMark(&sys, xi, yi, zi, &portCoorx, &portCoory);
	if (status == 0)
	{
		cout << "meshAndMark Success!" << endl;
		cout << "meshAndMark time is " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
	}
	else
	{
		cerr << "meshAndMark Fail!" << endl;
		return status;
	}

	// Set D_eps and D_sig
	clock_t t3 = clock();
	status = matrixConstruction(&sys);
	if (status == 0)
	{
		cout << "matrixConstruction Success!" << endl;
		cout << "matrixConstruction time is " << (clock() - t3) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
}
	else {
		cerr << "matrixConstruction Fail!" << endl;
		return status;
	}

	// Set port
	clock_t t4 = clock();
	status = portSet(&sys, xi, yi, zi);
	if (status == 0)
	{
		cout << "portSet Success!" << endl;
		cout << "portSet time is " << (clock() - t4) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
	}
	else
	{
		cerr << "portSet Fail!" << endl;
		return status;
	}

	// Generate Stiffness Matrix
	clock_t t5 = clock();
	status = generateStiff(&sys);
	if (status == 0)
	{
		cout << "generateStiff Success!" << endl;
		cout << "generateStiff time is " << (clock() - t5) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
	}
	else
	{
		cerr << "generateStiff Fail!" << endl;
		return status;
	}

	Solve_E_Zpara_InPardiso(&sys);

	return 0;
}

#ifndef SKIP_LAYERED_FDTD
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