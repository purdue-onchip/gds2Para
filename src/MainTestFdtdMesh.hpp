#ifndef GDS2PARA_MAINTESTFDTDMESH_H_
#define GDS2PARA_MAINTESTFDTDMESH_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

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