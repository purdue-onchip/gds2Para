#ifndef GDS2PARA_PARDISO_SOLVER_H_
#define GDS2PARA_PARDISO_SOLVER_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

// Cal all the computed freq points and store in a vector
vector<double> CalAllFreqPointsHz(const fdtdMesh &sys) {
	vector<double> vFreqHz;

	vFreqHz.push_back(sys.freqStart * sys.freqUnit);          // First frequency in sweep

	for (int id = 1; id < sys.nfreq; id++) {                  // When nfreq > 1
		if (sys.freqScale == 1) {
			vFreqHz.push_back((sys.freqStart + id * (sys.freqEnd - sys.freqStart) / (sys.nfreq - 1)) * sys.freqUnit);
		}	// Linear interpolation of frequency sweep
		else {
			vFreqHz.push_back(sys.freqStart * sys.freqUnit * pow(sys.freqEnd / sys.freqStart, (id * 1.0 / (sys.nfreq - 1))));
		}	// Logarithmic interpolation of frequency sweep
	}

	return vFreqHz;
}

// Assign source current density for a source port
void AssignSourceCurrentForSourcePort(fdtdMesh *psys, int sourcePort) {
	if (psys->J != nullptr) {
		free(psys->J);
	}	// delete previous J assignment

	psys->J = (double*)calloc(psys->N_edge, sizeof(double));
	for (int sourcePortSide = 0; sourcePortSide < psys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
		for (int indEdge = 0; indEdge < psys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
			/* Set current density for all edges within sides in port to prepare solver */
			/*cout << " port #" << sourcePort + 1 << ", side #" << sourcePortSide + 1 
				<< ", edge #" << sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge] 
				<< ": J = " << sys->portCoor[sourcePort].portDirection[sourcePortSide] << " A/m^2" << endl;*/
			psys->J[psys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]] = psys->portCoor[sourcePort].portDirection[sourcePortSide];
		}
	}
}

// Cal E field {e} by solving (-w^2*D_eps+ i wD_sig + S)x = -iwJ in Pardiso
int CalElectrFieldInPardiso(fdtdMesh *psys, complex<double> *peField, complex<double> *pRhsJ, double freqHz,
	myint *RowId, myint *ColId, double *val) {

	double omegaHz = freqHz * 2. * M_PI;

	myint size = psys->N_edge - 2 * psys->N_edge_s;
	myint *RowId1 = (myint*)malloc((size + 1) * sizeof(myint));
	myint count = 0;
	myint indi = 0;
	myint k = 0;
	complex<double> *valc;
	valc = (complex<double>*)calloc(psys->leng_S, sizeof(complex<double>));

	// Cal effective RHS current density -iwJ
	for (indi = psys->N_edge_s; indi < psys->N_edge - psys->N_edge_s; indi++) {
		pRhsJ[indi - psys->N_edge_s] = 0. + (1i) * -psys->J[indi] * omegaHz;
	}

	RowId1[k] = 0;
	k++;
	myint start;
	myint nnz = psys->leng_S;
	cout << "Start to generate CSR form for S!" << endl;
	indi = 0;
	while (indi < nnz) {
		start = RowId[indi];
		while (indi < nnz && RowId[indi] == start) {
			valc[indi] += val[indi]; // val[indi] is real
			if (RowId[indi] == ColId[indi]) {
				// valc[indi] needs omega * (-omega * epsilon  + indj * sigma) added to it
				complex<double> addedPart(-omegaHz * psys->eps[RowId[indi] + psys->N_edge_s], psys->sig[RowId[indi] + psys->N_edge_s]);
				valc[indi] += omegaHz * addedPart;
			}
			count++;
			indi++;
		}
		RowId1[k] = (count);
		k++;
	}
	cout << endl;

	myint mtype = 13;    /* Real complex unsymmetric matrix */
	myint nrhs = 1;    /* Number of right hand sides */
	void *pt[64];

	/* Pardiso control parameters */
	myint iparm[64];
	myint maxfct, mnum, phase, error, msglvl, solver;
	double dparm[64];
	int v0csin;
	myint perm;

	/* Number of processors */
	int num_procs;

	/* Auxiliary variables */
	char *var;

	msglvl = 0;    /* print statistical information */
	solver = 0;    /* use sparse direct solver */
	error = 0;
	maxfct = 1;
	mnum = 1;
	phase = 13;

	pardisoinit(pt, &mtype, iparm);
	iparm[38] = 1;
	iparm[34] = 1;    // 0-based indexing
					  //iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */

	cout << "Begin to solve (-w^2*D_eps + iwD_sig + S)x = -iwJ" << endl;
	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, valc, RowId1, ColId, &perm, &nrhs, iparm, &msglvl, pRhsJ, peField, &error);
	if (error != 0) {
		cout << endl << "Error During numerical factorization: " << error << endl;
		exit(2);
	}
	cout << "Solving (-w^2*D_eps+ i wD_sig + S)x = -iwJ is complete!" << endl;

	free(RowId1); RowId1 = NULL;
	free(valc); valc = NULL;
	return 0;
}

// Solve E field and Z-parameters in Pardiso
void Solve_E_Zpara_InPardiso(fdtdMesh *psys) {
	
	// All computed freq points
	vector<double> vFreqHz = CalAllFreqPointsHz(*psys);

	// Right-Hand-Side term -jwJ and electric field eField (V/m)
	complex<double> *pRhsJ = (complex<double>*)malloc((psys->N_edge - 2 * psys->N_edge_s) * sizeof(complex<double>));  // -jwJ
	complex<double> *peField = (complex<double>*)calloc((psys->N_edge - 2 * psys->N_edge_s), sizeof(complex<double>));

	for (int i = 0; i < vFreqHz.size(); i++) {
		for (int sourcePort = 0; sourcePort < psys->numPorts; sourcePort++) {
			AssignSourceCurrentForSourcePort(psys, sourcePort);
			CalElectrFieldInPardiso(psys, peField, pRhsJ, vFreqHz[i], psys->SRowId, psys->SColId, psys->Sval);
		}	// solve for each port
	}	// solve for each frequency


	free(pRhsJ);
	free(peField);
}


// Cal Z-parameters from E-field
void CalZParaFromField(fdtdMesh *psys, const lapack_complex_double *eField) {

	///* Build final network parameters matrix for this sourcePort by looping over response ports and adding contributions from each response port edge on first side */
	///* Z_ij = V_i / I_j = - integ[(grad V_i - part A_i / part t) * dl_i] / iinteg[(J_j) * dS_j] = - sum[(yd + yh)_respedge * leng_respedge]_oneside / sum[1 * area_source_side] */
	//for (int indPort = 0; indPort < psys->numPorts; indPort++) {
	//	int indPortSide = 0; // Only deal with first port side to get response edge line integral
	//	for (int indEdge = 0; indEdge < psys->portCoor[indPort].portEdge[indPortSide].size(); indEdge++) {
	//		myint thisEdge = psys->portCoor[indPort].portEdge[indPortSide][indEdge];
	//		if (thisEdge % (psys->N_edge_s + psys->N_edge_v) >= psys->N_edge_s) {    // This edge is along the z-axis
	//			inz = thisEdge / (psys->N_edge_s + psys->N_edge_v);
	//			leng = psys->zn[inz + 1] - psys->zn[inz];
	//		}
	//		else if (thisEdge % (psys->N_edge_s + psys->N_edge_v) >= (psys->N_cell_y) * (psys->N_cell_x + 1)) {    // This edge is along the x-axis
	//			inx = ((thisEdge % (psys->N_edge_s + psys->N_edge_v)) - (psys->N_cell_y) * (psys->N_cell_x + 1)) / (psys->N_cell_y + 1);
	//			leng = psys->xn[inx + 1] - psys->xn[inx];
	//		}
	//		else {    // This edge is along the y-axis
	//			iny = (thisEdge % (psys->N_edge_s + psys->N_edge_v)) % psys->N_cell_y;
	//			leng = psys->yn[iny + 1] - psys->yn[iny];
	//		}

	//		/*leng = pow((psys->nodepos[psys->edgelink[thisEdge * 2] * 3] - psys->nodepos[psys->edgelink[thisEdge * 2 + 1] * 3]), 2);
	//		leng += pow((psys->nodepos[psys->edgelink[thisEdge * 2] * 3 + 1] - psys->nodepos[psys->edgelink[thisEdge * 2 + 1] * 3 + 1]), 2);
	//		leng += pow((psys->nodepos[psys->edgelink[thisEdge * 2] * 3 + 2] - psys->nodepos[psys->edgelink[thisEdge * 2 + 1] * 3 + 2]), 2);
	//		leng = sqrt(leng);*/
	//		psys->x[indPort + psys->numPorts * xcol] -= yd[thisEdge] * leng; // Accumulating responses due to each response edge line integral (V)
	//	}

	//	/* Only divide matrix entry by current at end of response port calculation */
	//	//cout << "  leng = " << leng << ", first side portArea = " << psys->portCoor[sourcePort].portArea[0] << " m^2, first side portDirection = " << psys->portCoor[sourcePort].portDirection[0] << endl;
	//	double sourceCurrent = 0.; // In-phase current from unit source port edge current densities into supply point (A)
	//	for (int sourcePortSide = 0; sourcePortSide < psys->portCoor[sourcePort].multiplicity; sourcePortSide++)
	//	{
	//		sourceCurrent += psys->portCoor[sourcePort].portArea[sourcePortSide];
	//	}
	//	//cout << "  Response port voltage = " << psys->x[indPort + psys->numPorts * xcol] << " V, Source port current = " << sourceCurrent << " A" << endl;
	//	psys->x[indPort + psys->numPorts * xcol] /= sourceCurrent; // Final matrix entry (ohm)
	//}
}


#endif