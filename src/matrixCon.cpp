//#include "stdafx.h"
#include <ctime>
#include "fdtd.hpp"
#include "hypreSolver.h"
#include <mkl_service.h>
#include <mkl_dfti.h>
//#define TIME (1)   // TIME = 1 run time domain
//#define FREQ (1)   // define FREQ run frequency domain
#define SOLVECHIP (1)   // define SOLVECHIP in the frequency domain
//#define GENERATE_V0_SOLUTION
#define SKIP_VH
//#define INSIDE_OUTSIDE
void addOmegaEpsilon(fdtdMesh* sys, myint* RowId, myint* ColId, double* val, myint leng, myint N, double freq, double* val1);
void solveV0(double freq, double* rhs, complex<double>* u0, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, myint leng_v0d, sparse_matrix_t& v0ct, sparse_matrix_t& v0cat, myint leng_v0c, myint nedge, myint* AdRowId, myint* AdColId, double* Adval, myint leng_Ad, myint* AcRowId, myint* AcColId, double* Acval, myint leng_Ac, char ri);   // solve (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0) with real rhs
void V0Multiply(sparse_matrix_t& V0dt, sparse_matrix_t& V0ct, myint lengv0d, myint lengv0c, myint row, complex<double>* u1, complex<double>* u2);
void V0aMultiply(fdtdMesh* sys, sparse_matrix_t& V0dat, sparse_matrix_t& V0cat, myint lengv0d, myint lengv0c, myint nedge, complex<double>* u2, complex<double>* u1);
void generateP(fdtdMesh* sys, double freq, myint* LRowId, myint* LColId, double* Lval, myint leng_L, myint N, myint* PRowId, myint* PColId, double* Pval, myint& leng_P);
void addDiagV0bV0ba(fdtdMesh* sys, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo);
int gmres_V0dLoo(double* Doosval, sparse_matrix_t& V0dt, double* v0dn, sparse_matrix_t& V0dat, double* v0dan, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* A22dbRowId, myint* A22dbColId, double* A22dbval, myint leng_A22db, double* bm, double* res, myint N, myint N1, double freq);
void AtimesV(sparse_matrix_t& V0dt, double* v0dn, sparse_matrix_t& V0dat, double* vodan, double* Doosval, double freq, myint N1, myint N2, double* V1, double* V2, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22);
void PrecondUpperTri(double freq, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* A22dbRowId, myint* A22dbColId, double* A22dbval, myint leng_A22db, sparse_matrix_t& v0dt, double* v0dn, sparse_matrix_t& v0dat, double* v0dan, double* Doosval, double* V1, double* V2, myint N1, myint N2);
int gmres_solveA22(myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* A22dbRowId, myint* A22dbColId, double* A22dbval, myint leng_A22db, double* bm, double* x, myint N);
int combine_x_v0d_uh(double* x, fdtdMesh* sys, complex<double>* xr);
void extractColumnVectors(myint* RowId, myint* ColId, double* val, double* aval, myint num, myint lengs, myint lenge, myint* eRowId, myint* eColId, double* eval, double* eaval, myint& e_num);
void transposeMatrix(sparse_matrix_t& V, myint row, myint col, sparse_matrix_t& Vt);
int CSR2COO(myint* rowStart, myint* rowEnd, myint Rows, myint* RowId);
void gmres_V0dV0bLoo(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* bm, double* res, void* pt[64], myint iparm[64]);
void solveAooMatrix(double* Jo, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* res, void* pt[64], myint iparm[64]);
void A3b3timesV(double* Doosval, double freq, sparse_matrix_t& v0ddt, sparse_matrix_t& v0ddat, myint N1, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* v, double* res);
void Precond3b3UpperTri(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* rhs, double* x, void* pt[64], myint iparm[64]);
int combine_x_v0d_v0b_uh(sparse_matrix_t& V0ddt, sparse_matrix_t& V0bt, myint N1, myint N2, myint N3, double* res, double* x);
void findNNZ(myint* ArowStart, myint* ArowEnd, myint* AcolId, double* Aval, myint& leng_A, myint ARows, myint ACols, myint** NrowId, myint** NcolId, double** Nval);
void pardisoSol(myint* RowId1, myint* ColId, double* val, void* pt[64], myint iparm[64], double* b, double* x, myint N);
void solveMatrix(double* Jo, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint leng_v0d, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint leng_v0b, sparse_matrix_t& Soo, myint outside, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, myint* SoiRowId, myint* SoiColId, double* Soival, myint lengoi, myint* SioRowId, myint* SioColId, double* Sioval, myint lengio, myint* mapio, complex<double>* xt, myint nedge, void* pt[64], myint iparm[64]);
void solveAooMatrix_schur(double* Jo, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* res, void* pt[64], myint iparm[64]);
int solveA11A22(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]);
void solveA11A22_schur_A22(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]);
int gmres_solveA33schur(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, void* pt[64], myint iparm[64]);
void A33schurTimesV(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, void* pt[64], myint iparm[64], double* v1, double* v2);
void A11A22Precond(double freq, double* Doosval, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, myint N2, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64], double* v1, double* v2);
void A11A22timesV(double* Doosval, double freq, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, double* v, double* res);
void solveA11A22_schur(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]);
int gmres_solveA11schur(double* rhs, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]);
void A11schurtimesV(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64], double* v1, double* v2);
void DTimesV(fdtdMesh* sys, myint nedge, double freq, complex<double>* u1, complex<double>* u2);
double V2norm(complex<double>* r0, myint n);

static bool comp(pair<double, int> a, pair<double, int> b) {
	return a.first <= b.first;
};

int paraGenerator(fdtdMesh *sys, unordered_map<double, int>& xi, unordered_map<double, int>& yi, unordered_map<double, int>& zi) {
	int indi, indj, mark, k, l, n;
	int status = 0;
	int count = 0;
	int xcol = 0;
	vector<int> rowId;
	vector<int> colId;
	vector<int> temp2;
	double *temp;
	myint inx = 0, iny = 0, inz = 0;

	/* Construct V0d with row id, col id and its val */
	sys->leng_v0d1 = 0, sys->v0d1num = 0;    // store the num of v0d1 vectors, which are nodes outside the conductors
	sys->leng_v0d1a = 0, sys->v0d1anum = 0;
	myint leng_v0d2 = 0, leng_v0d2a = 0, v0d2num = 0, v0d2anum = 0;
	sys->mapd = (myint*)calloc(sys->N_node, (sizeof(myint)));
	double block1_x, block1_y, block2_x, block2_y, block3_x, block3_y;
	double sideLen = 0.; // around the conductor 10um is considered with rapid potential change
	myint node1 = 0, node2 = 0;

	block1_x = 0;// (sys->xlim2 - sys->xlim1) / 5 * sys->lengthUnit;
	block1_y = 0;// (sys->ylim2 - sys->ylim1) / 5 * sys->lengthUnit;
	block2_x = 0;// (sys->xlim2 - sys->xlim1) / 5 * sys->lengthUnit;
	block2_y = 0;// (sys->ylim2 - sys->ylim1) / 5 * sys->lengthUnit;
	block3_x = 0;
	block3_y = 0;
#ifdef PRINT_V0D_BLOCKS
	cout << "V0d's block1_x and block1_y are " << block1_x << " " << block1_y << endl;
	cout << "V0d's block2_x and block2_y are " << block2_x << " " << block2_y << endl;
	cout << "V0d's block3_x and block3_y are " << block3_x << " " << block3_y << endl;
#endif

	/* Generate V0d1 */
	sys->merge_v0d1(block1_x, block1_y, block2_x, block2_y, block3_x, block3_y, sideLen);
	
#ifdef INSIDE_OUTSIDE
	sys->merge_v0db(block1_x, block1_y, block2_x, block2_y, block3_x, block3_y, sideLen);
	myint* v0bRowId = (myint*)malloc((sys->v0dbnum - sys->v0dnum) * sizeof(myint));
	myint* v0bColId = (myint*)malloc((sys->v0dbnum - sys->v0dnum) * sizeof(myint));
	double* v0bval = (double*)malloc((sys->v0dbnum - sys->v0dnum) * sizeof(double));
	double* v0baval = (double*)malloc((sys->v0dbnum - sys->v0dnum) * sizeof(double));
	ofstream out;

	myint leng_v0b = sys->leng_v0db - sys->leng_v0d, v0bnum = sys->v0dbnum - sys->v0dnum;
	cout << "leng_v0db is " << sys->leng_v0db << " and leng_v0d is " << sys->leng_v0d << endl;
	extractColumnVectors(sys->v0dbRowId, sys->v0dbColId, sys->v0dbval, sys->v0dbaval, sys->v0dbnum, sys->leng_v0d, sys->leng_v0db - 1, v0bRowId, v0bColId, v0bval, v0baval, v0bnum);    // extract the V0b part from V0db

	/* Begin to map the inside and outside edges */
	sys->mapEdgeInsideOutside();   // generate mapio and mapioR
	/* End of mapping the inside and outside edges */

#endif

	/* Generate A11 = V0ddan'*D_eps*V0ddn */
	cout << "nnz of V0d is " << sys->v0d1num << endl;
	sys->generateAdeps_v0d1(sys->mapd, sys->v0d1num, sys->v0d1num, sys->leng_v0d1);   // Ad = V0d1a'*D_eps*V0d1
	//sys->generateAdeps_v0db(sys->mapd, sys->v0dnum, sys->v0dnum, sys->leng_v0d);   // Ad = V0da'*D_eps*V0d
	//sys->generateMd(sys->mapd, sys->v0d1num, sys->v0d1num, sys->leng_v0d1);   // Md = V0d1a'*V0d1

	//double* test = (double*)calloc(sys->leng_v0d, sizeof(double));
	//double* test1 = (double*)calloc(sys->leng_v0d, sizeof(double));
	//test[0] = 1e+12;
	//hypreSolve(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, test, sys->leng_v0d, test1, 1, 3);
	//free(test); test = NULL;
	//free(test1); test1 = NULL;


	/* Print out V0d, V0da, V0c, V0ca */
	ofstream outfile;
	//outfile.open("V0d1.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0d1num; ++indi) {
	//	outfile << sys->mapio[sys->mapEdge[sys->v0d1RowId[indi]]] + 1 << " " << sys->v0d1ColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << sys->v0d1val[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0d1a.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0d1num; ++indi) {
	//	outfile << sys->mapio[sys->mapEdge[sys->v0d1RowId[indi]]] + 1 << " " << sys->v0d1ColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << sys->v0d1aval[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0db.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0dnum; ++indi) {
	//	outfile << sys->mapio[sys->mapEdge[sys->v0dbRowId[indi]]] + 1 << " " << sys->v0dbColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << sys->v0dbval[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0dba.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0dnum; ++indi) {
	//	outfile << sys->mapio[sys->mapEdge[sys->v0dbRowId[indi]]] + 1 << " " << sys->v0dbColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << sys->v0dbaval[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0b.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < v0bnum; ++indi) {
	//	outfile << sys->mapio[sys->mapEdge[v0bRowId[indi]]] + 1 << " " << v0bColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << v0bval[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0ba.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < v0bnum; ++indi) {
	//	outfile << sys->mapio[sys->mapEdge[v0bRowId[indi]]] + 1 << " " << v0bColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << v0baval[indi] << endl;
	//}
	//outfile.close();

	//cout << "Space matrices are output!\n";
	//outfile.open("Xg.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->leng_Xg; ++indi) {
	//	outfile << sys->XgRowId[indi] + 1 << " " << sys->XgColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << sys->Xgval[indi] << endl;
	//}
	//outfile.close();

#ifdef INSIDE_OUTSIDE
	//for (int ind = 0; ind < sys->v0dnum; ++ind) {   // normalize V0db and V0dba
	//	sys->v0dbval[ind] /= sys->v0dn[sys->v0dbColId[ind]];
	//	sys->v0dbaval[ind] /= sys->v0dan[sys->v0dbColId[ind]];
	//}
	//for (int ind = sys->v0dnum; ind < sys->v0dbnum; ++ind) {   // normalize the vectors v0b and v0ba
	//	sys->v0dbval[ind] /= sys->v0dn[sys->v0dbColId[ind]];
	//	sys->v0dbaval[ind] /= sys->v0dan[sys->v0dbColId[ind]];
	//	v0bval[ind - sys->v0dnum] /= sys->v0dn[sys->v0dbColId[ind]];
	//	v0baval[ind - sys->v0dnum] /= sys->v0dan[sys->v0dbColId[ind]];
	//}
#endif

	cout << "The number of outside conductor edge is " << sys->outside << endl;
	cout << "The number of inside conductor edge is " << sys->N_edge - sys->bden - sys->outside << endl;

	sys->v0d1valo = (double*)malloc(sys->v0d1num * sizeof(double));
	for (indi = 0; indi < sys->v0d1num; indi++) {
		sys->v0d1valo[indi] = sys->v0d1val[indi];
	}
	for (indi = 0; indi < sys->v0d1num; indi++) {    // compute sqrt(D_eps)*V0d1
		sys->v0d1val[indi] *= sqrt(sys->stackEpsn[(sys->v0d1RowId[indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
	}
	sys->v0d1avalo = (double*)malloc(sys->v0d1num * sizeof(double));
	for (indi = 0; indi < sys->v0d1num; indi++) {
		sys->v0d1avalo[indi] = sys->v0d1aval[indi];
	}
	for (indi = 0; indi < sys->v0d1num; indi++) {
		sys->v0d1aval[indi] *= sqrt(sys->stackEpsn[(sys->v0d1RowId[indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
	}

	/* store the csr colId for v0d1 */
	sys->v0d1ColIdo = (myint*)malloc(sys->v0d1num * sizeof(myint));
	for (indi = 0; indi < sys->v0d1num; indi++) {
		sys->v0d1ColIdo[indi] = sys->v0d1ColId[indi];
	}
	free(sys->v0d1ColId); sys->v0d1ColId = (myint*)malloc((sys->leng_v0d1 + 1) * sizeof(myint));
	status = COO2CSR_malloc(sys->v0d1ColIdo, sys->v0d1RowId, sys->v0d1val, sys->v0d1num, sys->leng_v0d1, sys->v0d1ColId);
	if (status != 0) {
		return status;
	}

	myint V0d_row = sys->N_edge;   // sys->outside   // the row size of V0d

#ifdef INSIDE_OUTSIDE
								   /* store the csr colId for v0db */
	sys->v0dbColIdo = (myint*)malloc(sys->v0dbnum * sizeof(myint));
	for (indi = 0; indi < sys->v0dbnum; indi++) {
		sys->v0dbColIdo[indi] = sys->v0dbColId[indi];
	}
	free(sys->v0dbColId); sys->v0dbColId = (myint*)malloc((sys->leng_v0db + 1) * sizeof(myint));
	status = COO2CSR_malloc(sys->v0dbColIdo, sys->v0dbRowId, sys->v0dbval, sys->v0dbnum, sys->leng_v0db, sys->v0dbColId);
	if (status != 0) {
		return status;
	}

	/* store the csr colId for v0b */
	myint* v0bColIdo = (myint*)calloc(v0bnum, sizeof(myint));
	for (indi = 0; indi < v0bnum; indi++) {
		v0bColIdo[indi] = v0bColId[indi];
	}
	free(v0bColId); v0bColId = (myint*)malloc((leng_v0b + 1) * sizeof(myint));
	status = COO2CSR_malloc(v0bColIdo, v0bRowId, v0bval, v0bnum, leng_v0b, v0bColId);
	if (status != 0) {
		return status;
	}


	/* Set v0db to be the outside v0db */
	for (indi = 0; indi < sys->v0dbnum; ++indi) {
		sys->v0dbRowId[indi] = sys->mapio[sys->mapEdge[sys->v0dbRowId[indi]]];
	}
	/* Set V0b to be outside V0b */
	for (indi = 0; indi < v0bnum; ++indi) {
		v0bRowId[indi] = sys->mapio[sys->mapEdge[v0bRowId[indi]]];
	}

	V0d_row = sys->outside;
	double* v0dbavalmu = (double*)malloc(sys->v0dbnum * sizeof(double));
	for (indi = 0; indi < sys->v0dbnum; ++indi) {
		v0dbavalmu[indi] = sys->v0dbaval[indi] / MU;
	}

	/* Set v0d to be the outside v0d */
	for (indi = 0; indi < sys->v0d1num; ++indi) {
		//if (sys->markEdge[sys->v0d1RowId[indi]])
		//	cout << "Something wrong with V0d!" << endl;
		sys->v0d1RowId[indi] = sys->mapio[sys->mapEdge[sys->v0d1RowId[indi]]];
	}


#endif
	double* v0d1avalmu = (double*)malloc(sys->v0d1num * sizeof(double));
	for (indi = 0; indi < sys->v0d1num; ++indi) {
		v0d1avalmu[indi] = sys->v0d1avalo[indi] / MU;
	}

	//cout << "Number of NNZ in V0d1 is " << sys->v0d1num << endl;

	sparse_status_t s;

	/* V0d^T's csr form handle for MKL : to construct Loos */
	sparse_matrix_t V0dt;
	s = mkl_sparse_d_create_csr(&V0dt, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, V0d_row, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1valo);

	/* V0da^T's csr form handle for MKL */
	sparse_matrix_t V0dat;
	s = mkl_sparse_d_create_csr(&V0dat, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, V0d_row, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1avalo);


	/* V0da^T/MU's csr form handle for MKL */
	sparse_matrix_t V0damut;
	s = mkl_sparse_d_create_csr(&V0damut, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, V0d_row, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, v0d1avalmu);

#ifdef INSIDE_OUTSIDE
	/* V0b^T's csr form handle for MKL (V0d on the dielectric nodes and around conductors) */
	sparse_matrix_t V0bt;   // V0b normalized
	s = mkl_sparse_d_create_csr(&V0bt, SPARSE_INDEX_BASE_ZERO, leng_v0b, V0d_row, &v0bColId[0], &v0bColId[1], v0bRowId, v0bval);

	sparse_matrix_t V0b;
	transposeMatrix(V0bt, leng_v0b, V0d_row, V0b);   // get the transpose of V0bt

	myint v0bRows, v0bCols;
	MKL_INT *v0browStart, *v0browEnd, *v0bcolId;
	MKL_INT *v0b2colId;
	double* v0bval1;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	sparse_status_t s0;

	s0 = mkl_sparse_d_export_csr(V0b, &indexing, &v0bRows, &v0bCols, &v0browStart, &v0browEnd, &v0bcolId, &v0bval1);
	myint ii = 0, rowi = 0;
	//out.open("v0b.txt", ofstream::out | ofstream::trunc);
	//while (rowi < v0bRows) {
	//	while (ii < v0browEnd[rowi]) {
	//		out << rowi + 1 << " " << v0bcolId[ii] + 1 << " ";
	//		out << setprecision(15) << v0bval1[ii] << endl;
	//		ii++;
	//	}
	//	rowi++;
	//}
	////for (ii = 0; ii < v0bnum; ++ii) {
	////	out << v0bRowId[ii] + 1 << " " << v0bColIdo[ii] + 1 << " ";
	////	out << setprecision(15) << v0bval[ii] << endl;
	////}
	//out.close();

	/* V0ba^T's csr form handle for MKL (V0da on the dielectric nodes and around conductors) */
	sparse_matrix_t V0bat;   // V0ba normalized
	s = mkl_sparse_d_create_csr(&V0bat, SPARSE_INDEX_BASE_ZERO, leng_v0b, V0d_row, &v0bColId[0], &v0bColId[1], v0bRowId, v0baval);

	/* V0dd^T csr form handle for MKL */
	sparse_matrix_t V0ddt;   // V0dd normalized
	s = mkl_sparse_d_create_csr(&V0ddt, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d, V0d_row, &sys->v0dbColId[0], &sys->v0dbColId[1], sys->v0dbRowId, sys->v0dbval);

	/* V0dba^T csr form handle for MKL */
	sparse_matrix_t V0ddat;   // V0dda normalized
	s = mkl_sparse_d_create_csr(&V0ddat, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d, V0d_row, &sys->v0dbColId[0], &sys->v0dbColId[1], sys->v0dbRowId, sys->v0dbaval);
#endif

	int *argc;
	char ***argv;
	MPI_Init(argc, argv);

	//status = setHYPREMatrix(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_v0d1, ad, parcsr_ad);
	/* End */

	//sys->AdRowId1 = (int*)malloc((sys->leng_v0d1 + 1) * sizeof(int));
	//status = COO2CSR_malloc(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, sys->leng_v0d1, sys->AdRowId1);
	//if (status != 0)

	//    return status;

	//cout << "The number of nonzeros in Ad is " << sys->leng_Ad << endl;

	/* Construct V0c with row id, col id and its val */
	sys->leng_v0c = 0, sys->v0cnum = 0;
	sys->leng_v0ca = 0, sys->v0canum = 0;

	int numPortCdt = 0;
	count = 0;
	sys->cindex.push_back(-1);    // the last index in the sparse form for each conductor in V0c, indi denotes ith conductor (starting from 1)
	sys->acu_cnno.push_back(0);    // how many v0c are in each conductor (accumulative), indi denotes the ith conductor (starting from 1)


	count = 0;
	block1_x = 0;// (sys->xlim2 - sys->xlim1) / 3 * sys->lengthUnit;// (sys->xn[sys->nx - 1] - sys->xn[0]) / 10;
	block1_y = 0;// (sys->ylim2 - sys->ylim1) / 3 * sys->lengthUnit;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 10;
	block2_x = 0;// (sys->xlim2 - sys->xlim1) / 5 * sys->lengthUnit;
	block2_y = 0;// (sys->ylim2 - sys->ylim1) / 5 * sys->lengthUnit;

	clock_t t1 = clock();

	sys->mapc = (myint*)calloc(sys->N_node, sizeof(myint));
	sys->merge_v0c(block1_x, block1_y, block2_x, block2_y);

	cout << "Length of V0c is " << sys->leng_v0c << " number of non-zeros in V0c is " << sys->v0cnum << endl;
	cout << "Length of V0ca is " << sys->leng_v0ca << " number of non-zeros in V0ca is " << sys->v0canum << endl;
	cout << "V0c is generated!" << endl;

	indj = 0;

	/*for (indi = 0; indi < sys->numCdt + 1; indi++) {
	cout << sys->acu_cnno[indi] << " ";
	}
	cout << endl;*/
	t1 = clock();
	sys->generateAc(sys->mapc, sys->v0cnum, sys->v0canum, sys->leng_v0c);

#ifdef PRINT_VERBOSE_TIMING
	//cout << "Time to generate Ac is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
	//cout << "Number of non-zeros in Ac is " << sys->leng_Ac << endl;
#endif

	for (indi = 0; indi < sys->numCdt; indi++) {
		free(sys->conductor[indi].node); sys->conductor[indi].node = NULL;
	}
	free(sys->conductor); sys->conductor = NULL;
	/*  trial of first set HYPRE matrix Ac */
	//HYPRE_IJMatrix ac;
	//HYPRE_ParCSRMatrix parcsr_ac;
	//status = setHYPREMatrix(sys->AcRowId, sys->AcColId, sys->Acval, sys->leng_v0c, ac, parcsr_ac);
	/* End */


	sys->v0cvalo = (double*)malloc(sys->v0cnum * sizeof(double));
	for (indi = 0; indi < sys->v0cnum; indi++)
		sys->v0cvalo[indi] = sys->v0cval[indi];    // v0cvalo is the v0c values without D_sig
	for (indi = 0; indi < sys->v0cnum; indi++) {
		if (sys->markEdge[sys->v0cRowId[indi]] != 0) {
			sys->v0cval[indi] *= sqrt(SIGMA);       // Compute the sparse form of D_sig*V0c
		}
	}
	sys->v0cavalo = (double*)malloc(sys->v0canum * sizeof(double));
	for (indi = 0; indi < sys->v0canum; indi++)
		sys->v0cavalo[indi] = sys->v0caval[indi];
	for (indi = 0; indi < sys->v0canum; indi++) {
		if (sys->markEdge[sys->v0cRowId[indi]] != 0) {
			sys->v0caval[indi] *= sqrt(SIGMA);
		}
	}

	sys->v0cColIdo = (myint*)malloc(sys->v0cnum * sizeof(myint));
	for (indi = 0; indi < sys->v0cnum; indi++) {
		sys->v0cColIdo[indi] = sys->v0cColId[indi];
	}
	free(sys->v0cColId); sys->v0cColId = (myint*)malloc((sys->leng_v0c + 1) * sizeof(myint));
	status = COO2CSR_malloc(sys->v0cColIdo, sys->v0cRowId, sys->v0cval, sys->v0cnum, sys->leng_v0c, sys->v0cColId);
	if (status != 0) {
		return status;
	}


	//outfile.open("Ad.txt", ofstream::trunc | ofstream::out);
	//for (indi = 0; indi < sys->leng_Ad; ++indi) {
	//	outfile << sys->AdRowId[indi] + 1 << " " << sys->AdColId[indi] + 1 << " ";
	//	outfile << setprecision(15) << sys->Adval[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0c.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0cnum; ++indi) {
	//	outfile << sys->v0cRowId[indi] + 1 << " " << sys->v0cColIdo[indi] + 1 << " " << sys->v0cvalo[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0ca.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0cnum; ++indi) {
	//	outfile << sys->v0cRowId[indi] + 1 << " " << sys->v0cColIdo[indi] + 1 << " " << sys->v0cavalo[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("inside.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->N_edge; ++indi) {
	//	if (sys->markEdge[indi]) {
	//		outfile << indi << endl;
	//	}
	//}
	//outfile.close();
	/*for (indi = 0; indi < sys->numCdt + 1; indi++) {
	cout << sys->acu_cnno[indi] << " ";
	}
	cout << endl;
	for (indi = 0; indi < sys->numCdt + 1; indi++) {
	cout << sys->cindex[indi] << " ";
	}
	cout << endl;*/


	double *a;
	int *ia, *ja;


	/* Pick up a y0c2 cooresponding to one source port */
	complex<double> Zresult;
	double *bd1, *bd2;
	double *bdc1, *bdc2;
	double mdone;
	int ione;
	lapack_int *ipiv;
	int info;
	double *workspace;
	double *xd2;
	double *temp1;
	int startCol;
	double freq; // Frequency

	startCol = 0;
	sys->Y = (complex<double>*)calloc(sys->numPorts * sys->numPorts * sys->nfreq, sizeof(complex<double>));
	sys->x.assign(sys->numPorts * sys->numPorts * sys->nfreq, complex<double>(0., 0.)); // Use complex double constructor to assign initial output matrix for single-frequency solve

	double *b, *xc;    // the array of the right hand side
	xcol = 0;

	int port = 0, sourcePort = 0;    // show which port it is using
	int node;
	double *crhs;

	char transa;
	int m;
	double alpha = 1, beta = 0;
	char matdescra[6];
	matdescra[0] = 'G'; matdescra[3] = 'C';    // general matrix multi, 0-based indexing
	cout << endl << "Begin to solve for network parameters!" << endl;

	double *ydcp;
	double *y0c, *y0cs;
	double *yc, *yca;
	double *yccp;
	double *dRhs, *dRhs2, *crhss;
	complex<double> *yd;
	double *yd1, *yd2, *yd2a;
	double *v0caJ, *v0daJ, *y0d, *ydt, *y0d2, *ydat;
	lapack_complex_double *u0, *u0a;
	double *u0d, *u0c;
	double nn, nna;
	complex<double>* xr;
	struct matrix_descr descr;


	/* V0ca^T's csr form handle for MKL */
	sparse_matrix_t V0cat;
	s = mkl_sparse_d_create_csr(&V0cat, SPARSE_INDEX_BASE_ZERO, sys->leng_v0c, sys->N_edge, &sys->v0cColId[0], &sys->v0cColId[1], sys->v0cRowId, sys->v0cavalo);

	/* V0c^T's csr form handle for MKL */
	sparse_matrix_t V0ct;
	s = mkl_sparse_d_create_csr(&V0ct, SPARSE_INDEX_BASE_ZERO, sys->leng_v0c, sys->N_edge, &sys->v0cColId[0], &sys->v0cColId[1], sys->v0cRowId, sys->v0cvalo);


	lapack_complex_double *tmp;
	lapack_complex_double *m_h, *m_hc;
	lapack_complex_double *rhs_h, *rhs_h0;
	lapack_complex_double *J;
	lapack_complex_double *y_h;
	lapack_int info1;
	lapack_int iter;
	complex<double> *final_x;
	lapack_complex_double *J_h;
	double *ferr, *berr;


	/* Begin to generate Soo, Soi, Sio */
#ifdef INSIDE_OUTSIDE
	myint lengoo = 0, lengoi = 0, lengio = 0;   // nnz for each submatrix
	matrixInsideOutside_count(sys->mapio, sys->outside, sys->inside, sys->SRowId, sys->SColId, sys->leng_S, lengoo, lengoi, lengio);   // count nnz in each submatrix
	myint* SooRowId = (myint*)malloc(lengoo * sizeof(myint));
	myint* SooColId = (myint*)malloc(lengoo * sizeof(myint));
	double* Sooval = (double*)malloc(lengoo * sizeof(double));
	myint* SoiRowId = (myint*)malloc(lengoi * sizeof(myint));
	myint* SoiColId = (myint*)malloc(lengoi * sizeof(myint));
	double* Soival = (double*)malloc(lengoi * sizeof(double));
	myint* SioRowId = (myint*)malloc(lengio * sizeof(myint));
	myint* SioColId = (myint*)malloc(lengio * sizeof(myint));
	double* Sioval = (double*)malloc(lengio * sizeof(double));
	matrixInsideOutside(sys->mapio, sys->outside, sys->inside, sys->SRowId, sys->SColId, sys->Sval, sys->leng_S, SooRowId, SooColId, Sooval, SoiRowId, SoiColId, Soival, SioRowId, SioColId, Sioval);   // assign each matrix's rowId, colId, val
																																																		/* Transfer Soo to CSR format with SooRowId1 */
	myint* SooRowId1 = (myint*)malloc((sys->outside + 1) * sizeof(myint));
	status = COO2CSR_malloc(SooRowId, SooColId, Sooval, lengoo, sys->outside, SooRowId1);
	/* End of transferring Soo to CSR format with SooRowId1 */
	double *Sooval_old; //dj
	Sooval_old = (double*)malloc(lengoo * sizeof(double));
	for (int ind = 0; ind < lengoo; ++ind) {
		Sooval_old[ind] = Sooval[ind];
	}

	double* Doosval;
	Doosval = (double*)calloc(sys->outside, sizeof(double));
	//out.open("Doos.txt", ofstream::trunc | ofstream::out);
	for (int ind = 0; ind < sys->outside; ++ind) {
		Doosval[ind] = sys->getEps(sys->mapEdgeR[sys->mapioR[ind]]);
		//out << setprecision(15) << Doosval[ind] << endl;
	}
	//out.close();


	/* calculate D_epsoo+dt^2*Soo */
	//outfile.open("Soo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < lengoo; ++ind) {
	//Sooval[ind] = Sooval[ind] * pow(DT, 2);
	//if (SooRowId[ind] == SooColId[ind]) {
	//Sooval[ind] += sys->getEps(sys->mapEdgeR[sys->mapioR[SooRowId[ind]]]);
	//}
	//outfile << SooRowId[ind] + 1 << " " << SooColId[ind] + 1 << " ";
	//outfile << setprecision(15) << Sooval[ind] << endl;
	//}
	//outfile.close();
	/* End of generating Soo, Soi, Sio */

	/* Begin to generate Loo */
	sparse_matrix_t A, Adb;

	s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, V0dt, V0damut, &A);   // A = V0d*V0da'/MU + V0b*V0ba'/MU
																		   //s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, V0dbt, V0dbat, &Adb);   // Adb = v0d*v0da'/MU + V0b*V0ba'/MU
																		   //sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
																		   //myint ARows, ACols, leng_Aoo = 0;
																		   //MKL_INT *ArowStart, *ArowEnd, *AcolId;
																		   //double *Aval;
																		   //s0 = mkl_sparse_d_export_csr(A, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
																		   ////cout << "V0d*V0da' row number is " << ARows << " and column number is " << ACols << endl;
																		   //leng_Aoo = ArowEnd[ARows - 1];    // how many non-zeros are in A
																		   //myint* AoorowId1 = (myint*)malloc((sys->outside + 1) * sizeof(myint));
																		   ////out.open("AoorowId1.txt", std::ofstream::out | std::ofstream::trunc);
																		   //for (myint ind = 0; ind < sys->outside; ++ind) {
																		   //	AoorowId1[ind] = ArowStart[ind];
																		   //	//out << AoorowId1[ind] << endl;
																		   //}
																		   //AoorowId1[sys->outside] = ArowEnd[sys->outside - 1];
																		   ////out << AoorowId1[sys->outside] << endl;
																		   ////out.close();
																		   //for (myint ind = 0; ind < leng_Aoo; ++ind) {
																		   //	Aval[ind] /= MU;
																		   //}


																		   //if (status != 0) {
																		   //	return status;
																		   //}
	myint* LooRowId = NULL, *LoodbRowId = NULL;    // Loos = Soos+V0d*V0da'/mu+V0b*V0ba'/mu
												   //myint* LooRowId1 = NULL;
	myint* LooColId = NULL, *LoodbColId = NULL;
	double* Looval = NULL, *Loodbval = NULL;
	myint leng_Loo, leng_Loodb;
	sparseMatrixSum(sys, A, SooRowId1, SooColId, Sooval, sys->outside, &LooRowId, &LooColId, &Looval, leng_Loo);   // Soos + V0d*V0da'/MU
																												   //sparseMatrixSum(sys, Adb, SooRowId1, SooColId, Sooval, sys->outside, &LoodbRowId, &LoodbColId, &Loodbval, leng_Loodb);   // Soos + V0d*V0da'/MU + V0b*V0ba'/MU

																												   /* output the Loo */

																												   //addDiagV0bV0ba(sys, LoodbRowId, LoodbColId, Loodbval, leng_Loodb);    // generate Soos+V0d*V0da'/mu+V0b*V0ba'/mu-diag(V0b*V0ba'/mu)

																												   //out.open("Loo.txt", ofstream::out | ofstream::trunc);
																												   //for (int ind = 0; ind < leng_Loo; ++ind) {
																												   //	out << LooRowId[ind] + 1 << " " << LooColId[ind] + 1 << " ";
																												   //	out << setprecision(15) << Looval[ind] << endl;
																												   //}
																												   //out.close();
																												   //cout << "Loo is generated!\n";
																												   //out.open("Loodb.txt", ofstream::out | ofstream::trunc);
																												   //for (int ind = 0; ind < leng_Loodb; ++ind) {
																												   //	out << LoodbRowId[ind] + 1 << " " << LoodbColId[ind] + 1 << " ";
																												   //	out << setprecision(15) << Loodbval[ind] << endl;
																												   //}
																												   //out.close();

	myint* LooRowId1 = (myint*)calloc(sys->outside + 1, sizeof(myint));
	COO2CSR_malloc(LooRowId, LooColId, Looval, leng_Loo, sys->outside, LooRowId1);
	double* Looval_old = (double*)calloc(leng_Loo, sizeof(double));
	for (int ind = 0; ind < leng_Loo; ind++) {
		Looval_old[ind] = Looval[ind];
	}
	free(Looval); Looval = NULL;   // for future use Looval

	//myint* LoodbRowId1 = (myint*)calloc((sys->outside + 1), sizeof(myint));
	//status = COO2CSR_malloc(LoodbRowId, LoodbColId, Loodbval, leng_Loodb, sys->outside, LoodbRowId1);
	//double* Loodbval_old = (double*)calloc(leng_Loodb, sizeof(double));
	//for (int ind = 0; ind < leng_Loodb; ind++) {
	//	Loodbval_old[ind] = Loodbval[ind];
	//}
	//free(Loodbval); Loodbval = NULL;   // for future use Loodbval
	////sparseMatrixSum1(sys, AoorowId1, AcolId, Aval, SooRowId1, SooColId, Sooval, sys->outside);



#endif
	myint nedge = sys->N_edge - sys->bden;

	//outfile.open("epsoo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < sys->outside; ++ind) {
	//	outfile << sys->getEps(sys->mapEdgeR[sys->mapioR[ind]]) << endl;
	//}
	//outfile.close();
	/* End of generating Doo */

	///* free the space of V0d*V0da' */
	////free(AoorowId1); AoorowId1 = NULL;
	//mkl_sparse_destroy(A);
	/* end of free of the space V0d*V0da' */



	/* generate Laplacian matrix for this mesh */
	/* Generate the Laplacian matrix */
	myint leng = 0, leng_L;
	myint *LrowId, *LcolId;
	double *Lval;
	leng = generateLaplacian_count(sys);   // count how many nnz in the matrix
	LrowId = (myint*)malloc(leng * sizeof(myint));
	LcolId = (myint*)malloc(leng * sizeof(myint));
	Lval = (double*)malloc(leng * sizeof(double));
	status = generateLaplacian(sys, LrowId, LcolId, Lval);
	leng_L = leng;

	//outfile.open("L.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < leng; ++ind) {
	//	outfile << LrowId[ind] + 1 << " " << LcolId[ind] + 1 << " ";
	//	outfile << setprecision(15) << Lval[ind] << endl;
	//} 
	//outfile.close();
	//cout << "L is generated!\n";
	/* End of generating the Laplacian matrix */

	/* debug to test L to be solved with HYPRE */
	//double* Lr = (double*)calloc(sys->N_edge - sys->bden, sizeof(double));
	//Lr[0] = 1;
	//double* Lx = (double*)calloc(sys->N_edge - sys->bden, sizeof(double));
	//status = hypreSolve(LrowId, LcolId, Lval, leng, Lr, sys->N_edge - sys->bden, Lx, 1, 3);
	//free(Lr); Lr = NULL;
	//free(Lx); Lx = NULL;
	/* end of debugging to test L to be solved with HYPRE */


#ifdef TIME
	/* Time domain method to simulate */
	/* Generating left hand matrix (D_eps + dt * D_sig + dt ^ 2 * S) */
	double* sval = (double*)malloc(sys->leng_S * sizeof(double));   // sys->SRowId, sys->SColId, val, COO format of this matrix
	myint* SRowId1 = (myint*)malloc((nedge + 1) * sizeof(myint));    // matrix's CSR form's rowId, with the start index in each row
	sys->backDiffLeftMatrix(sys->SRowId, sys->SColId, sys->Sval, sys->leng_S, sval, DT);
	status = COO2CSR_malloc(sys->SRowId, sys->SColId, sval, sys->leng_S, nedge, SRowId1);
	/* End of generating left hand matrix (D_eps + dt * D_sig + dt ^ 2 * S) */


	double* val = (double*)malloc(leng * sizeof(double));
	sys->backDiffLeftMatrix(LrowId, LcolId, Lval, leng, val, DT);     // left hand matrix (D_eps + dt * D_sig + dt ^ 2 * L), where L is Laplacian matrix, result is put in val
																	  //cout << "Begin to generate Ll matrix!\n";
	sys->generateLaplacianLeft(LrowId, LcolId, Lval, leng, DT);    // Compute the matrix [V0a'*D*V0, V0a'*D; D*V0, D+L] and put in sys->LlrowId, sys->LlcolId, sys->Llval
																   //cout << "Already generate L matrix!\n";
																   //outfile.open("Ll.txt", std::ofstream::out | std::ofstream::trunc);
																   //for (int ind = 0; ind < sys->leng_Ll; ++ind) {
																   //   	outfile << sys->LlRowId[ind] + 1 << " " << sys->LlColId[ind] + 1 << " " << std::setprecision(15) << sys->Llval[ind] << endl;
																   //}
																   //outfile.close();
																   //cout << "Ll is generated!\n";

																   /* transfer Ll's rowId to csr format */
	sys->LlRowIdo = (myint*)malloc(sys->leng_Ll * sizeof(myint));
	for (int indi = 0; indi < sys->leng_Ll; indi++) {
		sys->LlRowIdo[indi] = sys->LlRowId[indi];
	}
	free(sys->LlRowId); sys->LlRowId = (myint*)malloc((sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden + 1) * sizeof(myint));
	status = COO2CSR_malloc(sys->LlRowIdo, sys->LlColId, sys->Llval, sys->leng_Ll, sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden, sys->LlRowId);
	if (status != 0) {
		return status;
	}

	sparse_matrix_t Ll;
	s = mkl_sparse_d_create_csr(&Ll, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden, sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden, &sys->LlRowId[0], &sys->LlRowId[1], sys->LlColId, sys->Llval);

	/* A11 rowId, colId, val and A22 rowId, colId, val */
	myint* A12RowId, *A12ColId, *A21RowId, *A21ColId, *A22RowId, *A22ColId;
	double* A12val, *A21val, *A22val;

	myint leng_A12 = 0, leng_A21 = 0, leng_A22 = 0;
	get1122Block_count(sys->leng_v0d1 + sys->leng_v0c, sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden, sys, leng_A12, leng_A21, leng_A22);
	A12RowId = (myint*)calloc(leng_A12, sizeof(myint));
	A12ColId = (myint*)calloc(leng_A12, sizeof(myint));
	A12val = (double*)calloc(leng_A12, sizeof(double));
	A21RowId = (myint*)calloc(leng_A21, sizeof(myint));
	A21ColId = (myint*)calloc(leng_A21, sizeof(myint));
	A21val = (double*)calloc(leng_A21, sizeof(double));
	A22RowId = (myint*)calloc(leng_A22, sizeof(myint));
	A22ColId = (myint*)calloc(leng_A22, sizeof(myint));
	A22val = (double*)calloc(leng_A22, sizeof(double));
	get1122Block(sys->leng_v0d1 + sys->leng_v0c, sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden, sys, A12RowId, A12ColId, A12val, A21RowId, A21ColId, A21val, A22RowId, A22ColId, A22val);
	//outfile.open("A22.txt", std::ofstream::trunc | std::ofstream::out);
	//for (int ind = 0; ind < leng_A22; ++ind) {
	//	outfile << A22RowId[ind] + 1 << " " << A22ColId[ind] + 1 << " " << A22val[ind] << endl;
	//}
	//outfile.close();

	for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {

		/* backward difference */
		find_Vh_back(sys, sourcePort, V0ct, V0cat, V0dt, V0dat, A12RowId, A12ColId, A12val, leng_A12, A21RowId, A21ColId, A21val, leng_A21, A22RowId, A22ColId, A22val, leng_A22, SRowId1, sys->SColId, sval, Ll);


		/* central difference */
		//find_Vh_central(sys, sourcePort);

	}
	free(SRowId1); SRowId1 = NULL;
	free(sval); sval = NULL;
	free(A12RowId); A12RowId = NULL;
	free(A12ColId); A12ColId = NULL;
	free(A12val); A12val = NULL;
	free(A21RowId); A21RowId = NULL;
	free(A21ColId); A21ColId = NULL;
	free(A21val); A21val = NULL;
	free(A22RowId); A22RowId = NULL;
	free(A22ColId); A22ColId = NULL;
	free(A22val); A22val = NULL;
	/* End of time domain method to simulate */
#endif

	/* Print out current for Matlab debug */
	outfile.open("J.txt", std::ofstream::out | std::ofstream::trunc);
	for (int sourcePort = sys->numPorts - 1; sourcePort < sys->numPorts; ++sourcePort) {
		double* cur = (double*)calloc((sys->N_edge), sizeof(double));   // current
		for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
			for (int indEdge = 0; indEdge < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
				/* Set current density for all edges within sides in port to prepare solver */
				cur[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]] = sys->portCoor[sourcePort].portDirection[sourcePortSide];
				//cout << "current left markEdge is " << sys->markEdge[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge] - sys->N_cell_y] << endl;
				//myint node1, node2;
				//int indx, indy, indz;
				//sys->compute_edgelink(sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge], node1, node2);
				//sys->compute_node_index(node1, indx, indy, indz);
				//cout << sys->markEdge[(indz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1) * (sys->N_cell_y + 1) + indy] << endl;
				//cout << sys->markEdge[(indz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1) * (sys->N_cell_y + 1) + indy + 1] << endl;
			}
		}
		for (int ind = 0; ind < sys->N_edge; ++ind) {
			outfile << cur[ind] << endl;
		}
		free(cur); cur = NULL;
	}
	outfile.close();
	cout << "current is generated!\n";

	/* Frequency domain to solve inside outside part */
	for (int ind = 0; ind < sys->nfreq; ++ind) {
#ifdef FREQ
		double freq = sys->freqNo2freq(ind);

		/* generate Loodb-omega^2*D_epsoo */
		double* Looval = (double*)malloc(leng_Loo * sizeof(double));
		addOmegaEpsilon(sys, LooRowId, LooColId, Looval_old, leng_Loo, sys->outside, freq, Looval);
		/* end of generating Loo-omega^2*D_epsoo */
		//outfile.open("A33.txt", std::ofstream::out | std::ofstream::trunc);
		//for (int indi = 0; indi < leng_Loo; ++indi) {
		//	outfile << LooRowId[indi] + 1 << " " << LooColId[indi] + 1 << " ";
		//	outfile << setprecision(15) << Looval[indi] << endl;
		//}
		//outfile.close();
		//cout << "A33 is generated!\n";


		//out.open("A11.txt", ofstream::out | ofstream::trunc);
		//for (int indi = 0; indi < sys->leng_Ad; ++indi) {
		//	out << sys->AdRowId[indi] + 1 << " " << sys->AdColId[indi] + 1 << " ";
		//	out << setprecision(15) << sys->Adval[indi] * (-1) * pow(freq * 2 * M_PI, 2) << endl;
		//}
		//out.close();

		//outfile.open("Loodb.txt", ofstream::out | ofstream::trunc);
		//for (int ind = 0; ind < leng_Loodb; ++ind) {
		//	outfile << LoodbRowId[ind] + 1 << " " << LoodbColId[ind] + 1 << " ";
		//	outfile << setprecision(15) << Loodbval[ind] << endl;
		//}
		//outfile.close();

		/* Frequency domain to solve the inside and outside part
		[-omega^2*D_epsoo+Soo, Soi;
		Sio,                  i*omega*D_sigii+Sii] */
		/* Begin to output (-omega^2*D_epsoo+Soo) */
		addOmegaEpsilon(sys, SooRowId, SooColId, Sooval_old, lengoo, sys->outside, freq, Sooval);

		sparse_matrix_t Soo;
		s0 = mkl_sparse_d_create_csr(&Soo, SPARSE_INDEX_BASE_ZERO, sys->outside, sys->outside, &SooRowId1[0], &SooRowId1[1], SooColId, Sooval);   // Soo denotes -omega^2*D_epsoo+Soo
		sparse_matrix_t sol, a22;
		s0 = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, Soo, V0b, &sol);    // sol=Soo*V0bn
		s0 = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, V0bat, sol, &a22);    // a22=V0ban'*Soo*V0bn


		myint A22Rows, A22Cols, leng_A22;
		MKL_INT *A22rowStart, *A22rowEnd;
		MKL_INT *A22colId;
		double* A22val;

		s0 = mkl_sparse_d_export_csr(a22, &indexing, &A22Rows, &A22Cols, &A22rowStart, &A22rowEnd, &A22colId, &A22val);
		myint* A22RowId, *A22ColId;
		double* A22Val;
		findNNZ(A22rowStart, A22rowEnd, A22colId, A22val, leng_A22, A22Rows, A22Cols, &A22RowId, &A22ColId, &A22Val);
		cout << "A22Rows and A22Cols " << A22Rows << " " << A22Cols << endl;

		myint* A22RowId1 = (myint*)calloc((A22Rows + 1), sizeof(myint));
		status = COO2CSR_malloc(A22RowId, A22ColId, A22Val, leng_A22, leng_v0b, A22RowId1);

		/* print out A */
		//out.open("a22.txt", ofstream::out | ofstream::trunc);
		//for (int indi = 0; indi < sys->leng_v0d + leng_v0b + sys->outside; ++indi) {
		//	double* v = (double*)calloc(sys->leng_v0d + leng_v0b + sys->outside, sizeof(double));
		//	double* A = (double*)calloc(sys->leng_v0d + leng_v0b + sys->outside, sizeof(double));
		//	v[indi] = 1;
		//	Precond3b3UpperTri(freq, Doosval, V0ddt, V0ddat, sys->leng_v0d, V0bt, V0bat, leng_v0b, Soo, sys->outside, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, v, A);
		//	for (int indj = 0; indj < sys->leng_v0d + leng_v0b + sys->outside; ++indj) {
		//		if (A[indj] != 0) {
		//			out << indj + 1 << " " << indi + 1 << " ";
		//			out << setprecision(15) << A[indj] << endl;
		//		}
		//	}
		//	free(A); A = NULL;
		//	free(v); v = NULL;
		//}
		//out.close();
		//cout << "P matrix inverse is generated!\n";

		//out.open("A22.txt", ofstream::trunc | ofstream::out);
		//for (myint indi = 0; indi < leng_A22; ++indi) {
		//	out << A22RowId[indi] + 1 << " " << A22ColId[indi] + 1 << " ";
		//	out << setprecision(15) << A22Val[indi] << endl;
		//}
		//out.close();
		/* End of outputting (-omega^2*D_epsoo+Soo) */

		/* begin to factorize A22 */

		myint mtype = 11;    /* Real unsymmetric matrix */
		myint nrhs = 1;    /* Number of right hand sides */
		void *pt[64];
		myint N2 = leng_v0b;

		/* Pardiso control parameters */
		myint iparm[64];
		myint maxfct, mnum, phase, error, msglvl, solver;
		double dparm[64];
		int v0csin;

		msglvl = 0;    /* print statistical information */
		solver = 0;    /* use sparse direct solver */
		error = 0;
		maxfct = 1;
		mnum = 1;
		phase = 12;    /* Analysis, numerical factorization */

		pardisoinit(pt, &mtype, iparm);
		iparm[38] = 1;
		iparm[34] = 1;    // 0-based indexing
		iparm[3] = 2;    // number of processors
						 //iparm[59] = 2;    // out of core version to solve very large problem
						 //iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */

		double * ddum;
		myint idum;    // integer dummy

					   /* debug testing to see the performance of different solvers */
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &N2, A22Val, A22RowId1, A22ColId, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		cout << "A22 factorization is done!\n";
		//complex<double>* xsol = (complex<double>*)calloc(nedge * sys->numPorts, sizeof(complex<double>));
		double* xsol = (double*)calloc(sys->outside, sizeof(double));
		for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
			double* Jo = (double*)calloc(sys->outside, sizeof(double));
			for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
				for (int indEdge = 0; indEdge < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
					/* Set current density for all edges within sides in port to prepare solver */
					Jo[sys->mapio[sys->mapEdge[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]]]] = sys->portCoor[sourcePort].portDirection[sourcePortSide] * (-1) * 2 * M_PI * freq;
				}
			}


			//status = hypreSolve(LooRowId, LooColId, Looval, leng_Loo, Jo, sys->outside, xsol, 1, 3);
			//mkl_gmres_A(Jo, xsol, LooRowId, LooColId, Looval_old, leng_Loo, sys->outside);

			/* solve the o-i system */
			//complex<double>* xt = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//solveMatrix(Jo, freq, Doosval, V0ddt, V0ddat, sys->leng_v0d, sys->v0dn, sys->v0dan, V0bt, V0bat, leng_v0b, Soo, sys->outside, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, SoiRowId, SoiColId, Soival, lengoi, SioRowId, SioColId, Sioval, lengoi, sys->mapio, xt, nedge, pt, iparm);   // solve the inside outside part
			//sys->Construct_Z_V0_Vh(xt, ind, sourcePort);
			//cout << "Port " << sourcePort << " and frequency is " << freq << endl;
			//for (int indi = 0; indi < nedge; ++indi) {
			//	xsol[sourcePort * nedge + indi] = xt[indi];
			//}
			//free(xt); xt = NULL;
			//out.open("antenna_time.txt", ofstream::out | ofstream::app);
			//cout << "Solve 3*3 time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
			//out.close();
			/* end of solving the o-i system */

			/* solve (-omega^2*D_epsoo+Soo)\(Jo) */
			double* res = (double*)calloc(sys->outside + sys->leng_v0d + leng_v0b, sizeof(double));
			//solveAooMatrix(Jo, freq, Doosval, V0ddt, V0ddat, sys->leng_v0d, sys->v0dn, sys->v0dan, V0bt, V0bat, leng_v0b, Soo, sys->outside, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, res, pt, iparm);
			solveAooMatrix_schur(Jo, freq, Doosval, V0ddt, V0ddat, sys->leng_v0d, sys->v0dn, sys->v0dan, V0bt, V0bat, leng_v0b, Soo, sys->outside, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, res, pt, iparm);
			//complex<double>* x = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			////	status = combine_x_v0d_uh(res, sys, x);
			//status = combine_x_v0d_v0b_uh(V0ddt, V0bt, sys->leng_v0d, leng_v0b, sys->outside, res, xsol);
			////out.open("plasma_pkg_time_47.5GHz.txt", ofstream::out | ofstream::app);
			////out << freq << " " << sourcePort << " " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
			////out.close();

			//for (int indi = 0; indi < sys->outside; ++indi) {
			//	x[sys->mapioR[indi]] = 0 + 1i * xsol[indi];
			//}
			//sys->Construct_Z_V0_Vh(x, ind, sourcePort);
			//cout << "Port " << sourcePort << " and frequency is " << freq << endl;
			//free(x); x = NULL;
			free(res); res = NULL;
			free(Jo); Jo = NULL;

		}
		phase = -1;     // Release internal memory
		pardiso(pt, &maxfct, &mnum, &mtype, &phase, &N2, A22Val, A22RowId1, A22ColId, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
		/* End of debugging to see the performance of different solvers */


		/* Only solve oo part with pardiso */
		reference_oo(sys, ind, SooRowId1, SooColId, Sooval, xsol);
		/* End of only solving oo part with pardiso */

		//t1 = clock();
		//reference(sys, ind, sys->SRowId, sys->SColId, sys->Sval, xsol);
		//cout << "Pardiso solve time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;

		/* start to solve with (-omega^2*D_epsoo+Soo) */
		//complex<double>* y;
		//t1 = clock();
		//y = (complex<double>*)calloc(nedge * sys->numPorts, sizeof(complex<double>));
		//solveFreqIO(sys, ind, y, SooRowId, SooRowId1, SooColId, Sooval, lengoo, SioRowId, SioColId, Sioval, lengio, LooRowId1, LooColId, Looval, leng_Loo);
		//for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
		//	sys->Construct_Z_V0_Vh(&y[sourcePort * nedge], ind, sourcePort);
		//}
		//free(y); y = NULL;
		//cout << "Direct solving (-omega^2*D_epsoo+Soo) time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
		/* end of solving with inside and outside */



		/* Solve [V0a'*D*V0, V0a'*D;
		D*V0,      D+L] with preconditioner which is the L part */
		//solveV0L(sys, ind, LrowId, LcolId, Lval, leng, (sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden) * 2, V0dt, V0dat, V0ct, V0cat);
		/* End of solving */

		/* Start to generate matrix [-omega^2*D_eps+S, -omega*D_sig;
		omega*D_sig,      -omega^2*D_eps+S] */
		//  myint leng_S1 = (sys->leng_S + sys->inside - sys->outside) * 2;   // nnz
		//  myint* S1RowId = (myint*)malloc(leng_S1 * sizeof(myint));
		//  myint* S1ColId = (myint*)malloc(leng_S1 * sizeof(myint));
		//  double* S1val = (double*)malloc(leng_S1 * sizeof(double));
		//  generateL1(sys, freq, sys->SRowId, sys->SColId, sys->Sval, sys->leng_S, nedge, S1RowId, S1ColId, S1val, leng_S1);
		//  /* End of generating matrix [-omega^2*D_eps+S, -omega*D_sig;
		//  omega*D_sig,      -omega^2*D_eps+S] */
		//for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
		//	double* J = (double*)calloc(nedge, sizeof(double));
		//	sys->setCurrent(sourcePort, J);
		//	for (int ind = 0; ind < nedge; ++ind) {
		//		J[ind] *= -freq * 2 * M_PI;   // -omega*J   the whole J
		//	}
		//	

		//	/* Solve (-omega^2*D_eps+1i*omega*D_sig+S) by gmres */
		//	double* xrri = (double*)calloc(2 * nedge, sizeof(double));
		//	double* bri = (double*)calloc(2 * nedge, sizeof(double));
		//	for (int ind = 0; ind < nedge; ++ind) {
		//		bri[nedge + ind] = J[ind];
		//	}
		//	mkl_gmres_A(bri, xrri, S1RowId, S1ColId, S1val, leng_S1, 2 * nedge);
		//	complex<double>* xr = (complex<double>*)calloc(nedge, sizeof(complex<double>));
		//	for (int ind = 0; ind < nedge; ++ind) {
		//		xr[ind] = xrri[ind] + 1i * xrri[ind + nedge];
		//	}
		//	sys->Construct_Z_V0_Vh(xr, ind, sourcePort);

		//	free(J); J = NULL;
		//	free(xrri); xrri = NULL;
		//	free(bri); bri = NULL;
		//	free(xr); xr = NULL;
		//}
		//free(S1RowId); S1RowId = NULL;
		//free(S1ColId); S1ColId = NULL;
		//free(S1val); S1val = NULL;
		/* End of solving -omega^2*D_eps+1i*omega*D_sig+S by gmres */


		free(Looval); Looval = NULL;
		free(Sooval); Sooval = NULL;
		//free(Loodbval); Loodbval = NULL;
		mkl_sparse_destroy(Soo);
		mkl_sparse_destroy(sol);
		mkl_sparse_destroy(a22);
		free(A22RowId1); A22RowId1 = NULL;
		free(A22RowId); A22RowId = NULL;
		free(A22ColId); A22ColId = NULL;
		free(A22Val); A22Val = NULL;
		free(xsol); xsol = NULL;

#endif

#ifdef SOLVECHIP
		double freq = sys->freqNo2freq(ind);
		cout << "Frequenc SOLVECHIP is " << freq << endl;
		/* Start to generate matrix [-omega^2*D_eps+L, -omega*D_sig;
		omega*D_sig,      -omega^2*D_eps+L] */
		//myint leng_L1 = (leng + sys->inside - sys->outside) * 2;   // sys->inside - sys->outside is the number of edges inside
		//myint* L1RowId = (myint*)malloc(leng_L1 * sizeof(myint));
		//myint* L1ColId = (myint*)malloc(leng_L1 * sizeof(myint));
		//double* L1val = (double*)malloc(leng_L1 * sizeof(double));
		//generateL1(sys, freq, LrowId, LcolId, Lval, leng, nedge, L1RowId, L1ColId, L1val, leng_L1);
		//outfile.open("L1.txt", std::ofstream::trunc | std::ofstream::out);
		//for (int ind = 0; ind < leng_L1; ++ind) {
		//	outfile << L1RowId[ind] + 1 << " " << L1ColId[ind] + 1 << " ";
		//	outfile << setprecision(15) << L1val[ind] << endl;
		//}
		//outfile.close();
		/* End of generating matrix [-omega^2*D_eps+L, -omega*D_sig;
		omega*D_sig,      -omega^2*D_eps+L] */


		/* Start to generate matrix [-omega^2*D_eps+S, -omega*D_sig;
		omega*D_sig,      -omega^2*D_eps+S] */
		//myint leng_S1 = (sys->leng_S + sys->inside - sys->outside) * 2;   // nnz
		//myint* S1RowId = (myint*)malloc(leng_S1 * sizeof(myint));
		//myint* S1ColId = (myint*)malloc(leng_S1 * sizeof(myint));
		//double* S1val = (double*)malloc(leng_S1 * sizeof(double));
		//generateL1(sys, freq, sys->SRowId, sys->SColId, sys->Sval, sys->leng_S, nedge, S1RowId, S1ColId, S1val, leng_S1);
		/* End of generating matrix [-omega^2*D_eps+S, -omega*D_sig;
		omega*D_sig,      -omega^2*D_eps+S] */

		for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
			double* J = (double*)calloc(nedge, sizeof(double));
			sys->setCurrent(sourcePort, J);
			for (int ind = 0; ind < nedge; ++ind) {
				J[ind] *= -freq * 2 * M_PI;   // -omega*J   the whole J
			}

			/* Solve -omega^2*D_eps+1i*omega*D_sig+S by gmres */
			//double* xrri = (double*)calloc(2 * nedge, sizeof(double));
			//double* bri = (double*)calloc(2 * nedge, sizeof(double));
			//for (int ind = 0; ind < nedge; ++ind) {
			//	bri[nedge + ind] = J[ind];
			//}
			//mkl_gmres_A(bri, xrri, L1RowId, L1ColId, L1val, leng_L1, 2 * nedge);
			////hypreSolve(L1RowId, L1ColId, L1val, leng_L1, bri, 2 * nedge, xrri, 1, 3);

			////mkl_gmres_A(bri, xrri, S1RowId, S1ColId, S1val, leng_S1, 2 * nedge);
			////complex<double>* xr = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			////for (int ind = 0; ind < nedge; ++ind) {
			////	xr[ind] = xrri[ind] + 1i * xrri[ind + nedge];
			////}
			////sys->Construct_Z_V0_Vh(xr, ind, sourcePort);
			////free(xrri); xrri = NULL;
			////free(bri); bri = NULL;
			////free(xr); xr = NULL;
			/* End of solving -omega^2*D_eps+1i*omega*D_sig+S by gmres */


			/* Solve [V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0, V0a'*(-omega^2*D_eps+1i*omega*D_sig);
			(-omega^2*D_eps+1i*omega*D_sig)*V0,      -omega^2*D_eps+1i*omega*D_sig+L] */
			
			double* V0aJ = (double*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(double));   // v0aJ = V0a'*(-omega*J)   question? can whole array put in multiply?
			alpha = 1;
			beta = 0;
			descr.type = SPARSE_MATRIX_TYPE_GENERAL;
			s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, J, beta, V0aJ);
			s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, J, beta, &V0aJ[sys->leng_v0d1]);
			//for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
			//	V0aJ[ind] /= sys->v0dan[ind];
			//	//if (sys->v0dan[ind] == 0)
			//	//	cout << "Wrong!\n";
			//}
			//for (int ind = 0; ind < sys->leng_v0c; ++ind) {
			//	V0aJ[sys->leng_v0d1 + ind] /= sys->v0can[ind];   // v0aJ= = V0a'*(-omega*J)
			//}
			complex<double>* u1 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			solveV0(freq, V0aJ, u1, V0dt, V0dat, sys->leng_v0d1, V0ct, V0cat, sys->leng_v0c, nedge, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, sys->AcRowId, sys->AcColId, sys->Acval, sys->leng_Ac, 'i');   // u1 = (V0a'*D*V0)\(V0a'*(-1i*omega*J))
			
			complex<double>* u2 = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			V0Multiply(V0dt, V0ct, sys->leng_v0d1, sys->leng_v0c, nedge, u1, u2);   // u2 = V0*(V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\(V0a'*1i*v0aJ)
			DTimesV(sys, nedge, freq, u2, u2);   // u2 = D*u2
			//double* u2ri = (double*)calloc(nedge * 2, sizeof(double));
			for (int ind = 0; ind < nedge; ++ind) {
				complex<double> temp(0, J[ind]);
				u2[ind] = temp - u2[ind];   // u2 = -1i*omega*J-(-omega^2*D_eps+1i*omega*D_sig)*V0*(V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\(V0a'*(-1i*omega*J))
				//u2ri[ind] = u2[ind].real();   // u2.real()
				//u2ri[nedge + ind] = u2[ind].imag();   // u2.imag()
			}

			complex<double>* uh = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			solve_L_schur(u2, uh, nedge, freq, Doosval, V0dt, V0dat, sys->leng_v0d1, V0ct, V0cat, sys->leng_v0c, LrowId, LcolId, Lval, leng_L);   // uh = ((D+L)-D*V0*(V0a'*D*V0)^(-1)*V0a'*D)\(-1i*omega*J-D*V0*(V0a'*D*V0)\(V0a'*(-1i*omega*J)))
			
			///* solve [-omega^2*D_eps+L, -omega*D_sig;
			//omega*D_sig,      -omega^2*D_eps+L]*[xr; xi] = [u2r; u2i] */
			//double* xhri = (double*)calloc(nedge * 2, sizeof(double));
			////hypreSolve(L1RowId, L1ColId, L1val, leng_L1, u2ri, nedge * 2, xhri, 1, 3);
			//mkl_gmres_A(u2ri, xhri, L1RowId, L1ColId, L1val, leng_L1, 2 * nedge);
			////hypreSolve(LrowId, LcolId, Lval, leng, u2ri, nedge, xhri, 1, 3);
			//complex<double>* xh = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//for (int ind = 0; ind < nedge; ++ind) {
			//	xh[ind] = xhri[ind] + 1i * xhri[nedge + ind];   // xh is got
			//}
			//

			///* (-omega^2*D_eps+1i*omega*D_sig)*xh */
			//for (int ind = 0; ind < nedge; ++ind) {
			//	complex<double> epsig;
			//	if (sys->markEdge[sys->mapEdgeR[ind]])
			//		epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * freq * 2 * M_PI * SIGMA;
			//	else
			//		epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * 0;
			//	u2[ind] = xh[ind] * epsig;   // u2 = (-omega^2*D_eps+1i*omega*D_sig)*xh
			//	u2[ind] = 1i * (J[ind] - u2[ind].imag()) - u2[ind].real();   // u2 = -1i*omega*J-(-omega^2*D_eps+1i*omega*D_sig)*xh
			//}
			//V0aMultiply(sys, V0dat, V0cat, sys->leng_v0d1, sys->leng_v0c, nedge, u2, u1);   // u1 = V0a'*(-1i*omega*J-(-omega^2*D_eps+1i*omega*D_sig)*xh)
			//double* tempr = (double*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(double));
			//double* tempi = (double*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(double));
			//complex<double>* u01 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			//complex<double>* u02 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			//for (int ind = 0; ind < sys->leng_v0d1 + sys->leng_v0c; ++ind) {
			//	tempr[ind] = u1[ind].real();
			//	tempi[ind] = u1[ind].imag();
			//}
			//solveV0(sys, freq, tempr, u01, V0dt, V0dat, V0ct, V0cat, 'r');   // u01 = (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\real(V0a'*(1i*v0aJ-(-omega^2*D_eps+1i*omega*D_sig)*xh))
			//solveV0(sys, freq, tempi, u02, V0dt, V0dat, V0ct, V0cat, 'i');   // u01 = (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\imag(V0a'*(1i*v0aJ-(-omega^2*D_eps+1i*omega*D_sig)*xh))
			//complex<double>* u0 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			//for (int ind = 0; ind < sys->leng_v0d1 + sys->leng_v0c; ++ind) {
			//	u0[ind] = u01[ind].real() + u02[ind].real() + 1i * (u01[ind].imag() + u02[ind].imag());
			//}
			//complex<double>* x0 = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//complex<double>* x = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//V0Multiply(V0dt, V0ct, sys->leng_v0d1, sys->leng_v0c, nedge, u0, x0);   // x0 = V0*u0
			//for (int ind = 0; ind < nedge; ++ind) {
			//	x[ind] = x0[ind].real() + xh[ind].real() + 1i * (x0[ind].imag() + xh[ind].imag());   // x = V0*u0+xh
			//}
			////outfile.open("xx.txt", ofstream::out | ofstream::trunc);
			////for (int ind = 0; ind < nedge; ++ind) {
			////	outfile << x[ind].real() << " " << x[ind].imag() << endl;
			////}
			////outfile.close();


			//sys->Construct_Z_V0_Vh(x, ind, sourcePort);

			/////* Calculate the error norm */
			//xr = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//sys->reference1(ind, sourcePort, xr);
			//double error = 0, b = 0;
			//for (int ind = 0; ind < nedge; ++ind) {
			//	error += pow(xr[ind].real() - x[ind].real(), 2) + pow(xr[ind].imag() - x[ind].imag(), 2);
			//	b += pow(xr[ind].real(), 2) + pow(xr[ind].imag(), 2);
			//}
			//cout << "Frequency " << freq << " error is " << error / b << endl;

			//free(xr); xr = NULL;
			/* End of calculating the error norm */

			free(J); J = NULL;
			//free(V0aJ); V0aJ = NULL;
			//free(u1); u1 = NULL;
			//free(u2); u2 = NULL;
			//free(u2ri); u2ri = NULL;
			//free(xhri); xhri = NULL;
			//free(xh); xh = NULL;
			//free(tempr); tempr = NULL;
			//free(tempi); tempi = NULL;
			//free(u01); u01 = NULL;
			//free(u02); u02 = NULL;
			//free(u0); u0 = NULL;
			//free(x0); x0 = NULL;
			//free(x); x = NULL;
			/* End of solving [V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0, V0a'*(-omega^2*D_eps+1i*omega*D_sig);
			(-omega^2*D_eps+1i*omega*D_sig)*V0,      -omega^2*D_eps+1i*omega*D_sig+L] */
		}
		//free(L1RowId); L1RowId = NULL;
		//free(L1ColId); L1ColId = NULL;
		//free(L1val); L1val = NULL;

#endif
	}
#ifdef GENERATE_V0_SOLUTION
	for (int sourcePort = 0; sourcePort < sys->numPorts; sourcePort++) {
		sys->J = (double*)calloc(sys->N_edge, sizeof(double));
		for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
			//double sourceArea = 0;
			//sourceArea += sys->portCoor[sourcePort].portArea[sourcePortSide];
			for (int indEdge = 0; indEdge < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
				/* Set current density for all edges within sides in port to prepare solver */
				sys->J[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]] = sys->portCoor[sourcePort].portDirection[sourcePortSide];
			}
			//cout << "SourcePort " << sourcePort << " area is " << sourceArea << endl;
		}


		v0daJ = (double*)calloc(sys->leng_v0d1, sizeof(double));
		y0d = (double*)calloc(sys->leng_v0d1, sizeof(double));

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;

		s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, sys->J, beta, v0daJ);
		for (indi = 0; indi < sys->leng_v0d1; indi++) {
			v0daJ[indi] *= -1.0;
		}
		/* solve V0d system */

		//status = hypreSolve(ad, parcsr_ad, sys->leng_Ad, v0daJ, sys->leng_v0d1, y0d);
		status = hypreSolve(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, v0daJ, sys->leng_v0d1, y0d, 1, 3);
		/* End of solving */


		for (indi = 0; indi < sys->leng_v0d1; indi++) {
			y0d[indi] /= (2 * M_PI * sys->freqStart * sys->freqUnit);    // y0d is imaginary
		}
		ydt = (double*)calloc(sys->N_edge, sizeof(double));
		ydat = (double*)calloc(sys->N_edge, sizeof(double));
		yd1 = (double*)malloc(sys->N_edge * sizeof(double));

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, y0d, beta, ydt);    // -V0d*(D_eps0\(V0da'*rsc))
		s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dat, descr, y0d, beta, ydat);    // -V0da*(D_eps0\(V0da'*rsc))

		u0 = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * 2, sizeof(lapack_complex_double));
		u0a = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * 2, sizeof(lapack_complex_double));
		nn = 0;
		nna = 0;
		for (indi = 0; indi < sys->N_edge; indi++) {
			nn += ydt[indi] * ydt[indi];
			nna += ydat[indi] * ydat[indi];
		}
		nn = sqrt(nn);
		nna = sqrt(nna);
		for (indi = 0; indi < sys->N_edge - sys->bden; indi++) {
			u0[indi].real = ydt[sys->mapEdgeR[indi]] / nn;    // u0d is one vector in V0
			u0a[indi].real = ydat[sys->mapEdgeR[indi]] / nna;    // u0da
		}

		/* Compute C right hand side */
		y0c = (double*)calloc(sys->leng_v0c, sizeof(double));
		v0caJ = (double*)calloc(sys->leng_v0c, sizeof(double));

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, sys->J, beta, v0caJ);
		// myint count_non = 0; // This is not used later in the function but was used for printing earlier
		for (indi = 0; indi < sys->leng_v0c; indi++) {
			v0caJ[indi] *= -1.0;
		}

		crhs = (double*)calloc(sys->leng_v0c, sizeof(double));
		for (indi = 0; indi < sys->N_edge; indi++) {
			yd1[indi] = ydt[indi];
			ydt[indi] *= -1.0 * (2 * M_PI*sys->freqStart * sys->freqUnit) * sys->stackEpsn[(indi + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
		}

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, ydt, beta, crhs);
		free(ydt); ydt = NULL;
		free(ydat); ydat = NULL;

		double v0caJn, crhsn;
		v0caJn = 0;
		crhsn = 0;
		for (indi = 0; indi < sys->leng_v0c; indi++) {
			v0caJ[indi] += crhs[indi];
			v0caJn += v0caJ[indi] * v0caJ[indi];
			crhsn += crhs[indi] * crhs[indi];
		}
		v0caJn = sqrt(v0caJn);
		crhsn = sqrt(crhsn);
		free(crhs); crhs = NULL;


		//cout << "Time between the first and the second HYPRE is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
		/*solve c system block by block*/

		for (indi = 1; indi <= 1; indi++) {

			//a = (double*)malloc((sys->cindex[indi] - sys->cindex[indi - 1]) * sizeof(double));
			//ia = (int*)malloc((sys->cindex[indi] - sys->cindex[indi - 1]) * sizeof(int));
			//ja = (int*)malloc((sys->cindex[indi] - sys->cindex[indi - 1]) * sizeof(int));
			//crhss = (double*)malloc((sys->acu_cnno[indi] - sys->acu_cnno[indi - 1]) * sizeof(double));
			//y0cs = (double*)calloc((sys->acu_cnno[indi] - sys->acu_cnno[indi - 1]), sizeof(double));

			//for (indj = sys->acu_cnno[indi - 1]; indj <= sys->acu_cnno[indi] - 1; indj++) {
			//    crhss[indj - sys->acu_cnno[indi - 1]] = v0caJ[indj];
			//}
			//for (indj = sys->cindex[indi - 1] + 1; indj <= sys->cindex[indi]; indj++) {
			//    a[indj - sys->cindex[indi - 1] - 1] = sys->Acval[indj];
			//    ia[indj - sys->cindex[indi - 1] - 1] = sys->AcRowId[indj] - sys->AcRowId[sys->cindex[indi - 1] + 1];
			//    ja[indj - sys->cindex[indi - 1] - 1] = sys->AcColId[indj] - sys->AcColId[sys->cindex[indi - 1] + 1];
			//    //cout << a[indj - sys->cindex[indi - 1] - 1] << " " << ia[indj - sys->cindex[indi - 1] - 1] << " " << ja[indj - sys->cindex[indi - 1] - 1] << endl;
			//}

			//status = hypreSolve(ia, ja, a, (sys->cindex[indi] - sys->cindex[indi - 1]), crhss, (sys->acu_cnno[indi] - sys->acu_cnno[indi - 1]), y0cs);

			//for (indj = sys->acu_cnno[indi - 1]; indj <= sys->acu_cnno[indi] - 1; indj++) {
			//    y0c[indj] = y0cs[indj - sys->acu_cnno[indi - 1]];
			//}
			//if (sys->acu_cnno[indi] - sys->acu_cnno[indi - 1] > 100000) {
			//    outfile1.open("Ac.txt", std::ofstream::out | std::ofstream::trunc);
			//    for (indj = 0; indj < sys->cindex[indi] - sys->cindex[indi - 1]; indj++) {
			//        outfile1 << sys->AcRowId[indi] + 1 << " " << sys->AcColId[indi] + 1 << " " << sys->Acval[indi] << endl;
			//    }
			//    outfile1.close();
			//}
			//free(a); a = NULL;
			//free(ia); ia = NULL;
			//free(ja); ja = NULL;
			//free(crhss); crhss = NULL;
			//free(y0cs); y0cs = NULL;

			//status = hypreSolve(ac, parcsr_ac, sys->leng_Ac, v0caJ, sys->leng_v0c, y0c);
			status = hypreSolve(sys->AcRowId, sys->AcColId, sys->Acval, sys->leng_Ac, v0caJ, sys->leng_v0c, y0c, 1, 3);

		}

		free(v0caJ); v0caJ = NULL;


		/* V0cy0c */
		yc = (double*)calloc(sys->N_edge, sizeof(double));
		yca = (double*)calloc(sys->N_edge, sizeof(double));
		yccp = (double*)malloc(sys->N_edge * sizeof(double));
		dRhs2 = (double*)calloc(sys->leng_v0d1, sizeof(double));
		y0d2 = (double*)calloc(sys->leng_v0d1, sizeof(double));

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ct, descr, y0c, beta, yc);
		s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0cat, descr, y0c, beta, yca);

		free(y0c); y0c = NULL;

		for (indi = 0; indi < sys->N_edge; indi++) {
			yccp[indi] = -yc[indi] * sys->stackEpsn[(indi + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
		}

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, yccp, beta, dRhs2);
		free(yccp); yccp = NULL;


		status = hypreSolve(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, dRhs2, sys->leng_v0d1, y0d2, 1, 3);
		free(dRhs2); dRhs2 = NULL;

		yd2 = (double*)calloc(sys->N_edge, sizeof(double));
		yd2a = (double*)calloc(sys->N_edge, sizeof(double));

		alpha = 1;
		beta = 0;
		descr.type = SPARSE_MATRIX_TYPE_GENERAL;
		s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, y0d2, beta, yd2);
		s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dat, descr, y0d2, beta, yd2a);
		free(y0d2); y0d2 = NULL;
		nn = 0;
		nna = 0;
		for (indi = 0; indi < sys->N_edge; indi++) {
			nn += (yd2[indi] + yc[indi]) * (yd2[indi] + yc[indi]);
			nna += (yd2a[indi] + yca[indi]) * (yd2a[indi] + yca[indi]);
		}
		nn = sqrt(nn);
		nna = sqrt(nna);

		for (indi = 0; indi < sys->N_edge - sys->bden; indi++) {
			u0[sys->N_edge - sys->bden + indi].real = (yd2[sys->mapEdgeR[indi]] + yc[sys->mapEdgeR[indi]]) / nn;    // u0c is the other vector in u0
			u0a[sys->N_edge - sys->bden + indi].real = (yd2a[sys->mapEdgeR[indi]] + yca[sys->mapEdgeR[indi]]) / nna;    // u0ca
		}

		yd = (complex<double>*)malloc(sys->N_edge * sizeof(complex<double>));
		for (int id = 0; id < sys->N_edge; id++) {
			yd[id] = yd2[id] - (1i)*(yd1[id]);
		}

		free(yd2); yd2 = NULL;
		free(yd1); yd1 = NULL;
		free(yd2a); yd2a = NULL;
		free(yca); yca = NULL;

		//sys->y = (complex<double>*)malloc(sys->N_edge*sizeof(complex<double>));
		for (indi = 0; indi < sys->N_edge; indi++) {
			yd[indi] += yc[indi];
		}
		free(yc); yc = NULL;


		/*out.open("U0.txt", std::ofstream::out | std::ofstream::trunc);
		for (indi = 0; indi < sys->N_edge - sys->bden; indi++) {
		for (int inde = 0; inde < 2; inde++) {
		out << u0[inde * (sys->N_edge - sys->bden) + indi].real << " ";
		}
		out << endl;
		}
		out.close();
		out.open("U0a.txt", std::ofstream::out | std::ofstream::trunc);
		for (indi = 0; indi < sys->N_edge - sys->bden; indi++) {
		for (int inde = 0; inde < 2; inde++) {
		out << u0a[inde * (sys->N_edge - sys->bden) + indi].real << " ";
		}
		out << endl;
		}
		out.close();*/
		sys->Construct_Z_V0(yd, sourcePort);
	}
#endif

	/* Calculate the Vh part */
#ifndef SKIP_VH
	cout << "Begin to solve Vh!\n";
	// u0a = V0a*(V0a'*V0a)*V0a'*u0
	/*for (indi = 0; indi < 2; indi++) {
	}*/

	// find the Vh eigenmodes
	//cout << "Begin to find Vh!\n";
	//status = find_Vh_central(sys, u0, u0a, sourcePort);
	//status = find_Vh_back(sys, sourcePort);
	//cout << "Finish finding Vh!\n";

	for (indi = 0; indi < sys->nfreq; indi++) {
		// this point's frequency

		freq = sys->freqNo2freq(indi);


		// Vh = Vh - u0*((u0a'*A*u0)\(u0a'*A*Vh))
		lapack_complex_double* V_re2 = (lapack_complex_double*)malloc((sys->N_edge - sys->bden) * sys->leng_Vh * sizeof(lapack_complex_double));
		for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++) {    // A*Vh
			for (myint inde2 = 0; inde2 < sys->leng_Vh; inde2++) {
				if (sys->markEdge[sys->mapEdgeR[inde]] != 0) {    // if this edge is inside the conductor
					V_re2[inde2 * (sys->N_edge - sys->bden) + inde].real = sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde])
						- 2 * M_PI * freq * sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].imag * SIGMA;
					V_re2[inde2 * (sys->N_edge - sys->bden) + inde].imag = sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].real * (freq * 2 * M_PI) * SIGMA
						- sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde]);
				}
				else {
					V_re2[inde2 * (sys->N_edge - sys->bden) + inde].real = sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde]);
					V_re2[inde2 * (sys->N_edge - sys->bden) + inde].imag = -sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde]);
				}
			}
		}


		lapack_complex_double *y_re = (lapack_complex_double*)calloc(2 * sys->leng_Vh, sizeof(lapack_complex_double));    // u0a'*A*Vh
		status = matrix_multi('T', u0a, (sys->N_edge - sys->bden), 2, V_re2, (sys->N_edge - sys->bden), sys->leng_Vh, y_re);
		free(V_re2); V_re2 = NULL;

		lapack_complex_double *tmp3 = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * 2, sizeof(lapack_complex_double));
		for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++) {    // A*u0
			for (myint inde2 = 0; inde2 < 2; inde2++) {
				if (sys->markEdge[sys->mapEdgeR[inde]] != 0) {
					tmp3[inde2 * (sys->N_edge - sys->bden) + inde].real = u0[inde2 * (sys->N_edge - sys->bden) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde])
						- 2 * M_PI * freq * u0[inde2 * (sys->N_edge - sys->bden) + inde].imag * SIGMA;
					tmp3[inde2 * (sys->N_edge - sys->bden) + inde].imag = u0[inde2 * (sys->N_edge - sys->bden) + inde].real * (freq * 2 * M_PI) * SIGMA
						- u0[inde2 * (sys->N_edge - sys->bden) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde]);
				}
				else {
					tmp3[inde2 * (sys->N_edge - sys->bden) + inde].real = u0[inde2 * (sys->N_edge - sys->bden) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde]);
					tmp3[inde2 * (sys->N_edge - sys->bden) + inde].imag = -u0[inde2 * (sys->N_edge - sys->bden) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->getEps(sys->mapEdgeR[inde]);
				}
			}
		}

		lapack_complex_double *tmp4 = (lapack_complex_double*)calloc(2 * 2, sizeof(lapack_complex_double));    // u0a'*A*u0
		status = matrix_multi('T', u0a, (sys->N_edge - sys->bden), 2, tmp3, (sys->N_edge - sys->bden), 2, tmp4);    // u0a'*A*u0
		ipiv = (lapack_int*)malloc(2 * sizeof(lapack_int));
		lapack_complex_double *y_new = (lapack_complex_double*)calloc(2 * sys->leng_Vh, sizeof(lapack_complex_double));
		/*outfile.open("ma.txt", std::ofstream::out | std::ofstream::trunc);
		for (myint inde = 0; inde < 2; inde++){
		for (myint inde1 = 0; inde1 < 2; inde1++){
		outfile << tmp4[inde1 * 2 + inde].real << " " << tmp4[inde1 * 2 + inde].imag << " ";
		}
		outfile << endl;
		}*/
		info = LAPACKE_zcgesv(LAPACK_COL_MAJOR, 2, sys->leng_Vh, tmp4, 2, ipiv, y_re, 2, y_new, 2, &iter);    // (u0a'*A*u0)\(u0a'*A*V_re)
		free(ipiv); ipiv = NULL;
		free(y_re); y_re = NULL;
		free(tmp3); tmp3 = NULL;
		free(tmp4); tmp4 = NULL;

		lapack_complex_double *m_new = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * sys->leng_Vh, sizeof(lapack_complex_double));
		status = matrix_multi('N', u0, (sys->N_edge - sys->bden), 2, y_new, 2, sys->leng_Vh, m_new);    // u0*((u0a'*A*u0)\(u0a'*A*V_re))
		free(y_new); y_new = NULL;
		lapack_complex_double* Vh = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * (sys->leng_Vh), sizeof(lapack_complex_double));
		for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++) {
			for (myint inde2 = 0; inde2 < sys->leng_Vh; inde2++) {
				Vh[inde2 * (sys->N_edge - sys->bden) + inde].real = sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].real - m_new[inde2 * (sys->N_edge - sys->bden) + inde].real;
				Vh[inde2 * (sys->N_edge - sys->bden) + inde].imag = sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].imag - m_new[inde2 * (sys->N_edge - sys->bden) + inde].imag;
			}
		}
		free(m_new); m_new = NULL;


		// Vh'*(A+C)*Vh
		int inde;
		myint start;
		tmp = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * sys->leng_Vh, sizeof(lapack_complex_double));
		for (myint j = 0; j < sys->leng_Vh; j++) {    // calculate (A+C)*V_re1
			inde = 0;
			while (inde < sys->leng_S) {
				start = sys->SRowId[inde];
				while (inde < sys->leng_S && sys->SRowId[inde] == start) {

					tmp[j * (sys->N_edge - sys->bden) + sys->SRowId[inde]].real += sys->Sval[inde] * Vh[j * (sys->N_edge - sys->bden) + sys->SColId[inde]].real;
					tmp[j * (sys->N_edge - sys->bden) + sys->SRowId[inde]].imag += sys->Sval[inde] * Vh[j * (sys->N_edge - sys->bden) + sys->SColId[inde]].imag;


					inde++;
				}
			}

			for (inde = 0; inde < sys->N_edge - sys->bden; inde++) {
				if (sys->markEdge[sys->mapEdgeR[inde]] != 0) {
					tmp[j * (sys->N_edge - sys->bden) + inde].real += -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * Vh[j * (sys->N_edge - sys->bden) + inde].real
						- freq * 2 * M_PI * SIGMA * Vh[j * (sys->N_edge - sys->bden) + inde].imag;
					tmp[j * (sys->N_edge - sys->bden) + inde].imag += -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * Vh[j * (sys->N_edge - sys->bden) + inde].imag
						+ freq * 2 * M_PI * SIGMA * Vh[j * (sys->N_edge - sys->bden) + inde].real;
				}
				else {
					tmp[j * (sys->N_edge - sys->bden) + inde].real += -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * Vh[j * (sys->N_edge - sys->bden) + inde].real;
					tmp[j * (sys->N_edge - sys->bden) + inde].imag += -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * Vh[j * (sys->N_edge - sys->bden) + inde].imag;
				}
			}
		}


		m_h = (lapack_complex_double*)calloc(sys->leng_Vh * sys->leng_Vh, sizeof(lapack_complex_double));
		status = matrix_multi('T', Vh, (sys->N_edge - sys->bden), sys->leng_Vh, tmp, (sys->N_edge - sys->bden), sys->leng_Vh, m_h);    // V_re1'*(A+C)*V_re1

		rhs_h = (lapack_complex_double*)calloc(sys->leng_Vh * 1, sizeof(lapack_complex_double));
		J = (lapack_complex_double*)calloc(sys->N_edge - sys->bden, sizeof(lapack_complex_double));
		for (inde = 0; inde < sys->N_edge; inde++) {
			if (sys->J[inde] != 0) {
				J[sys->mapEdge[inde]].imag = -1 * freq * 2 * M_PI;
			}
		}

		status = matrix_multi('T', Vh, (sys->N_edge - sys->bden), sys->leng_Vh, J, (sys->N_edge - sys->bden), 1, rhs_h);    // -1i*omega*V_re1'*J
		free(J); J = NULL;


		/* V_re1'*A*u */
		free(tmp);
		tmp = (lapack_complex_double*)calloc((sys->N_edge - sys->bden), sizeof(lapack_complex_double));
		for (inde = 0; inde < sys->N_edge - sys->bden; inde++) {
			if (sys->markEdge[sys->mapEdgeR[inde]] != 0) {
				tmp[inde].real = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * yd[sys->mapEdgeR[inde]].real() - freq * 2 * M_PI * SIGMA * yd[sys->mapEdgeR[inde]].imag() * sys->freqStart * sys->freqUnit / freq;
				tmp[inde].imag = freq * 2 * M_PI * SIGMA * yd[sys->mapEdgeR[inde]].real() - pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * yd[sys->mapEdgeR[inde]].imag() * sys->freqStart * sys->freqUnit / freq;
			}
			else {
				tmp[inde].real = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * yd[sys->mapEdgeR[inde]].real();
				tmp[inde].imag = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[inde]) * yd[sys->mapEdgeR[inde]].imag() * sys->freqStart * sys->freqUnit / freq;
			}
		}
		rhs_h0 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
		status = matrix_multi('T', Vh, sys->N_edge - sys->bden, sys->leng_Vh, tmp, sys->N_edge - sys->bden, 1, rhs_h0);    // V_re1'*A*u
		for (inde = 0; inde < sys->leng_Vh; inde++) {
			rhs_h[inde].real = rhs_h[inde].real - rhs_h0[inde].real;
			rhs_h[inde].imag = rhs_h[inde].imag - rhs_h0[inde].imag;
		}
		free(rhs_h0); rhs_h0 = NULL;
		free(tmp); tmp = NULL;


		ipiv = (lapack_int*)malloc(sys->leng_Vh * sizeof(lapack_int));
		info1 = LAPACKE_zgesv(LAPACK_COL_MAJOR, sys->leng_Vh, 1, m_h, sys->leng_Vh, ipiv, rhs_h, sys->leng_Vh);// , y_h, sys->leng_Vh, &iter);    // yh is generated
		free(ipiv); ipiv = NULL;
		free(m_h); m_h = NULL;

		y_h = (lapack_complex_double*)calloc((sys->N_edge - sys->bden), sizeof(lapack_complex_double));
		status = matrix_multi('N', Vh, (sys->N_edge - sys->bden), sys->leng_Vh, rhs_h, sys->leng_Vh, 1, y_h);

		final_x = (complex<double>*)malloc((sys->N_edge - sys->bden) * sizeof(complex<double>));
		for (inde = 0; inde < sys->N_edge - sys->bden; inde++) {
			final_x[inde] = yd[sys->mapEdgeR[inde]].real() + y_h[inde].real + 1i * (yd[sys->mapEdgeR[inde]].imag() *sys->freqStart * sys->freqUnit / freq + y_h[inde].imag);
		}

		free(y_h); y_h = NULL;
		free(rhs_h); rhs_h = NULL;

		xr = (complex<double>*)calloc(sys->N_edge - sys->bden, sizeof(complex<double>));

		sys->reference1(indi, sourcePort, xr);

		// Construct Z parameters
		sys->Construct_Z_V0_Vh(final_x, indi, sourcePort);


		/* Compare with the solution from (-w^2*D_eps+iw*D_sig+S)xr=-iwJ */
		double err = 0, total_norm = 0, err0 = 0;
		for (inde = 0; inde < sys->N_edge - sys->bden; inde++) {
			err += (xr[inde].real() - final_x[inde].real()) * (xr[inde].real() - final_x[inde].real()) + (xr[inde].imag() - final_x[inde].imag()) * (xr[inde].imag() - final_x[inde].imag());
			err0 += (xr[inde].real() - yd[sys->mapEdgeR[inde]].real()) * (xr[inde].real() - yd[sys->mapEdgeR[inde]].real()) + (xr[inde].imag() - yd[sys->mapEdgeR[inde]].imag() *sys->freqStart * sys->freqUnit / freq) * (xr[inde].imag() - yd[sys->mapEdgeR[inde]].imag() *sys->freqStart * sys->freqUnit / freq);
			total_norm += (xr[inde].real()) * (xr[inde].real()) + (xr[inde].imag()) * (xr[inde].imag());
		}
		err = sqrt(err);
		err0 = sqrt(err0);
		total_norm = sqrt(total_norm);

		cout << "Freq " << freq << " the y0 total error is " << err0 / total_norm << endl;
		cout << "Freq " << freq << " the total error is " << err / total_norm << endl;

		free(final_x); final_x = NULL;
		free(Vh); Vh = NULL;
		free(xr); xr = NULL;
	}
	free(sys->Vh); sys->Vh = NULL;
	free(u0); u0 = NULL;
	free(u0a); u0a = NULL;
#endif

#ifdef GENERATE_V0_SOLUTION
	free(yd); yd = NULL;
	free(sys->J); sys->J = NULL;
	/*free(u0); u0 = NULL;*/
	//free(sys->y); sys->y = NULL;
	xcol++;

#endif


	//free(Pval); Pval = NULL;

	MPI_Finalize();

#ifndef SKIP_STIFF_REFERENCE
	/*  Generate the reference results and S parameters in .citi file for different frequencies with multiple right hand side */


	for (indi = 0; indi < sys->nfreq; indi++) {

		if (sys->nfreq == 1) {    // to avoid (sys->nfreq - 1)
			freq = sys->freqStart * sys->freqUnit;
		}
		else {
			if (sys->freqScale == 1) {
				freq = (sys->freqStart + indi * (sys->freqEnd - sys->freqStart) / (sys->nfreq - 1)) * sys->freqUnit;
			}
			else {
				freq = sys->freqStart * sys->freqUnit * pow(sys->freqEnd / sys->freqStart, (indi * 1.0 / (sys->nfreq - 1)));
			}
		}
		cout << "Frequency " << freq << "'s z parameter matrix is shown below as" << endl;

		status = reference(sys, indi, sys->SRowId, sys->SColId, sys->Sval);


		/*for (int indj = 0; indj < sys->numPorts; indj++){
		for (int indk = 0; indk < sys->numPorts; indk++){
		cout << sys->x[indi * (sys->numPorts * sys->numPorts) + indj * sys->numPorts + indk] << " ";
		}
		cout << endl;
		}*/
	}

#endif

	/* Report the Z-parameters and Prepare to Export Them */
#ifdef PRINT_V0_Z_PARAM
	sys->print_z();
#endif

#ifdef PRINT_V0_Vh_Z_PARAM
	sys->print_z();
#endif
	sys->cindex.clear();
	sys->acu_cnno.clear();
	sys->stackEpsn.clear();


	mkl_sparse_destroy(V0dt);
	mkl_sparse_destroy(V0dat);
	mkl_sparse_destroy(V0damut);
	mkl_sparse_destroy(V0ct);
	mkl_sparse_destroy(V0cat);


	//mkl_sparse_destroy(V0dbt);
	//mkl_sparse_destroy(V0dbat);


#ifdef FREQ
	/* free Soo matrix */
	free(SooRowId); SooRowId = NULL;
	free(SooRowId1); SooRowId1 = NULL;
	free(SooColId); SooColId = NULL;
	free(Sooval); Sooval = NULL;
	//free(LooRowId); LooRowId = NULL;
	//free(LooColId); LooColId = NULL;
	//free(Looval_old); Looval_old = NULL;
	/* end of freeing Soo matrix */

	/* free Soi matrix */
	free(SoiRowId); SoiRowId = NULL;
	free(SoiColId); SoiColId = NULL;
	free(Soival); Soival = NULL;
	/* end of freeing Soi matrix */

	/* free Sio matrix */
	free(SioRowId); SioRowId = NULL;
	free(SioColId); SioColId = NULL;
	free(Sioval); Sioval = NULL;
	/* end of freeing Sio matrix */

	/* free Doos matrix */
	free(Doosval); Doosval = NULL;
	free(v0dbavalmu); v0dbavalmu = NULL;

	mkl_sparse_destroy(V0bt);
	mkl_sparse_destroy(V0b);
	mkl_sparse_destroy(V0bat);
	mkl_sparse_destroy(V0ddt);
	mkl_sparse_destroy(V0ddat);
#endif
	free(v0d1avalmu); v0d1avalmu = NULL;

	//HYPRE_IJMatrixDestroy(ad);
	//HYPRE_IJMatrixDestroy(ac);

	return 0;
}

void V0V0aTV0V0aTv(sparse_matrix_t& V0t, sparse_matrix_t& V0at, myint leng_v0, myint* MRowId, myint* MColId, double* Mval, myint leng_M, double* v1, double* v2) {
	/* apply V0d*(V0da'*V0d)*(V0da') to v1 and generate v2 */

	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	int status;

	double* tmp = (double*)calloc(leng_v0, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0at, descr, v1, beta, tmp);   // v0a'*v1
	double* sol = (double*)calloc(leng_v0, sizeof(double));
	status = hypreSolve(MRowId, MColId, Mval, leng_M, tmp, leng_v0, sol, 1, 3);   // (V0da'*V0d)\v0da'*v1

	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0t, descr, sol, beta, v2);   // v2=V0d*(V0da'*V0d)\(v0da'*v1)

	free(tmp); tmp = NULL;
	free(sol); sol = NULL;
}



int COO2CSR(vector<int> &rowId, vector<int> &ColId, vector<double> &val) {
	myint i;
	vector<int> rowId2;
	int count, start;

	rowId2.push_back(0);
	count = 0;
	i = 0;
	while (i < rowId.size()) {
		start = rowId[i];
		while (i < rowId.size() && rowId[i] == start) {
			count++;
			i++;
		}
		rowId2.push_back(count);
	}

	rowId.clear();
	rowId = rowId2;
	return 0;
}

int COO2CSR_malloc(myint *rowId, myint *ColId, double *val, myint totalnum, myint leng, myint *rowId1) {
	/* Transfer the COO format to CSR format
	The COO format is rowwise
	totalnum is the total number of entries
	leng is the row number
	*/
	myint i, start, k;

	i = 0;
	k = 0;
	rowId1[k] = 0;
	k++;
	while (i < totalnum) {
		start = rowId[i];   // start + 1 denotes the row # (starting from 1)
		while (k < start + 1) {
			rowId1[k] = i;
			k++;
		}
		while (i < totalnum && rowId[i] == start) {
			i++;
		}
		rowId1[start + 1] = i;
		k = start + 2;
	}

	return 0;
}



int setsideLen(int node, double sideLen, int *markLayerNode, int *markProSide, fdtdMesh *sys) {
	queue<int> q;
	int *visited;
	int indx, indy;
	int mark;
	double startx, starty;

	if (sideLen == 0)
		return 0;

	visited = (int*)calloc(sys->N_node_s, sizeof(int));
	q.push(node);
	visited[node] = 1;
	startx = sys->xn[node / (sys->N_cell_y + 1)];
	starty = sys->yn[node % (sys->N_cell_y + 1)];
	// bfs
	while (!q.empty()) {
		indx = q.front() / (sys->N_cell_y + 1);
		indy = q.front() % (sys->N_cell_y + 1);

		if (indx != sys->nx - 1) {    // it must have a right x edge, thus right x node

			if (visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0) {    // this node is in dielectric and this node is not visited
				if (sqrt(pow((sys->xn[indx + 1] - startx), 2) + pow((sys->yn[indy] - starty), 2)) <= sideLen) {    // this node is within the block area
					q.push((indx + 1)*(sys->N_cell_y + 1) + indy);
					visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
					markProSide[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
				}
			}
		}

		if (indx != 0) {    // it must have a left x edge, thus left x node
			if (visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0) {    // this node is in dielectric and this node is not visited
				if (sqrt(pow((sys->xn[indx - 1] - startx), 2) + pow((sys->yn[indy] - starty), 2)) <= sideLen) {    // this node is within the block area
					q.push((indx - 1)*(sys->N_cell_y + 1) + indy);
					visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
					markProSide[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
				}
			}
		}

		if (indy != sys->ny - 1) {    // it must have a farther y edge, thus farther y node
			if (visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 0) {    // this node is in dielectric and this node is not visited
				if (sqrt(pow((sys->xn[indx] - startx), 2) + pow((sys->yn[indy + 1] - starty), 2)) <= sideLen) {    // this node is within the block area
					q.push((indx)*(sys->N_cell_y + 1) + indy + 1);
					visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
					markProSide[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
				}
			}
		}

		if (indy != 0) {    // it must have a closer y edge, thus closer y node
			if (visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0) {    // this node is in dielectric and this node is not visited
				if (sqrt(pow((sys->xn[indx] - startx), 2) + pow((sys->yn[indy - 1] - starty), 2)) <= sideLen) {    // this node is within the block area
					q.push((indx)*(sys->N_cell_y + 1) + indy - 1);
					visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
					markProSide[(indx)* (sys->N_cell_y + 1) + indy - 1] = 1;
				}
			}
		}

		q.pop();
	}

	free(visited); visited = NULL;

	return 0;
}

void matrixInsideOutside_count(myint* mapio, myint outside, myint inside, myint* rowId, myint* colId, myint leng, myint& lengoo, myint& lengoi, myint& lengio) {
	/* count each submatrix nnz number
	sys : used to provide inside outside edge numbers
	rowId : matrix rowId
	colId : matrix colId
	leng: : nnz in this matrix
	lengoo : oo nnz
	lengoi : oi nnz
	lengio : io nnz */

	for (myint ind = 0; ind < leng; ++ind) {
		if (mapio[rowId[ind]] >= 0 && mapio[rowId[ind]] < outside && mapio[colId[ind]] >= 0 && mapio[colId[ind]] < outside) {
			lengoo++;
		}
		else if (mapio[rowId[ind]] >= 0 && mapio[rowId[ind]] < outside && mapio[colId[ind]] >= outside && mapio[colId[ind]] < inside) {
			lengoi++;
		}
		else if (mapio[rowId[ind]] >= outside && mapio[rowId[ind]] < inside && mapio[colId[ind]] >= 0 && mapio[colId[ind]] < outside) {
			lengio++;
		}
	}
	return;
}

void matrixInsideOutside(myint* mapio, myint outside, myint inside, myint* rowId, myint* colId, double* val, myint leng, myint* ooRowId, myint* ooColId, double* ooval, myint* oiRowId, myint* oiColId, double* oival, myint* ioRowId, myint* ioColId, double* ioval) {
	/* assign the matrix that is oo, oi, io
	rowId : matrix rowId
	colId : matrix colId
	val : matrix value
	leng : matrix nnz
	lengoo : count the nnz in oo
	lengoi : count the nnz in oi
	lengio : count the nnz in io*/
	myint lengoo = 0, lengoi = 0, lengio = 0;
	for (myint ind = 0; ind < leng; ++ind) {
		if (mapio[rowId[ind]] >= 0 && mapio[rowId[ind]] < outside && mapio[colId[ind]] >= 0 && mapio[colId[ind]] < outside) {
			ooRowId[lengoo] = mapio[rowId[ind]];
			ooColId[lengoo] = mapio[colId[ind]];
			ooval[lengoo] = val[ind];
			lengoo++;
		}
		else if (mapio[rowId[ind]] >= 0 && mapio[rowId[ind]] < outside && mapio[colId[ind]] >= outside && mapio[colId[ind]] < inside) {
			oiRowId[lengoi] = mapio[rowId[ind]];
			oiColId[lengoi] = mapio[colId[ind]] - outside;
			oival[lengoi] = val[ind];
			lengoi++;
		}
		else if (mapio[rowId[ind]] >= outside && mapio[rowId[ind]] < inside && mapio[colId[ind]] >= 0 && mapio[colId[ind]] < outside) {
			ioRowId[lengio] = mapio[rowId[ind]] - outside;
			ioColId[lengio] = mapio[colId[ind]];
			ioval[lengio] = val[ind];
			lengio++;
		}
	}
	return;
}

void matrixOutside_count(fdtdMesh* sys, myint* ArowStart, myint* ArowEnd, myint* AcolId, double* Aval, myint ARows, myint& leng_Aoo) {
	/* Count how many nnz in oo part of the matrix
	ArowStart : the start index in each row
	ArowEnd : the start index of the next row
	AcolId : colId of this matrix
	Aval : value of this matrix
	ARows : how many rows are there in the matrix
	leng_Aoo : the nnz of this matrix as output */
	for (myint ind = 0; ind < ARows; ++ind) {
		for (myint indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
			leng_Aoo++;
		}
	}
}

void matrixOutside(fdtdMesh* sys, myint* ArowStart, myint* ArowEnd, myint* AcolId, double* Aval, myint ARows, myint* AoorowId, myint* AoocolId, double* Aooval, double scale) {
	/* assign the matrix nnz for the ourside part
	ArowStart :  the start index in each row
	ArowEnd : the start index in the next row
	AcolId : colId of this matrix
	Aval : value of this matrix
	ARows : row number in this matrix
	AoorowId : output, the outside part rowId
	AoocolId : output, the outside part colId
	Aooval : output, the outside part of value
	scale : scalar scale of the value */
	myint leng_Aoo = 0;
	for (myint ind = 0; ind < ARows; ++ind) {
		for (myint indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
			AoorowId[leng_Aoo] = ind;
			AoocolId[leng_Aoo] = AcolId[indi];
			Aooval[leng_Aoo] = Aval[indi] * scale;
			leng_Aoo++;
		}
	}
}

void sparseMatrixSum(fdtdMesh* sys, sparse_matrix_t& A, myint* browId1, myint* bcolId, double* bval, myint Rows, myint** rRowId, myint** rColId, double** rval, myint& lengr) {
	/* Add two sparse matrices
	arowId1 : csr form of rowId of matrix a
	acolId : a matrix colId
	aval : a matrix value
	browId : csr form of rowId of matrix b
	bcolId : b matrix colId
	bval : b matrix value
	Rows : row # of matrices a or b
	rrowId : resultant matrix rowId
	rcolId : resultant matrix colId
	rval : resultant matrix value
	rleng : resultant matrix nnz */

	sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
	sparse_matrix_t a, b, r;
	sparse_status_t s0;
	double alpha = 1.0;

	//s0 = mkl_sparse_d_create_csr(&a, SPARSE_INDEX_BASE_ZERO, Rows, Rows, &arowId1[0], &arowId1[1], acolId, aval);
	////if (s0 == SPARSE_STATUS_SUCCESS)
	////	cout << "SPARSE_STATUS_SUCCESS\n";

	s0 = mkl_sparse_d_create_csr(&b, SPARSE_INDEX_BASE_ZERO, Rows, Rows, &browId1[0], &browId1[1], bcolId, bval);
	//if (s0 == SPARSE_STATUS_SUCCESS)
	//	cout << "SPARSE_STATUS_SUCCESS";

	s0 = mkl_sparse_d_add(operation, A, alpha, b, &r);   // add the two sparse matrices

	myint ARows, ACols;
	MKL_INT *ArowStart, *ArowEnd;
	MKL_INT *AcolId;
	double *Aval;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	s0 = mkl_sparse_d_export_csr(r, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
	lengr = ArowEnd[ARows - 1];

	checkCSRRepeats(ArowStart, ArowEnd, AcolId, Aval, lengr, sys->outside);

	*rRowId = (myint*)calloc(lengr, sizeof(myint));
	*rColId = (myint*)calloc(lengr, sizeof(myint));
	*rval = (double*)calloc(lengr, sizeof(double));

	myint count = 0;
	for (myint ind = 0; ind < ARows; ++ind) {
		for (myint indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
			(*rRowId)[count] = ind;
			(*rColId)[count] = AcolId[indi];
			(*rval)[count] = Aval[indi];
			count++;
		}
	}


	//mkl_sparse_destroy(a);
	mkl_sparse_destroy(b);
	mkl_sparse_destroy(r);
}

void sparseMatrixSum1(fdtdMesh* sys, myint* arowId1, myint* acolId, double* aval, myint* browId1, myint* bcolId, double* bval, myint Rows) {
	/* Add two sparse matrices
	arowId1 : csr form of rowId of matrix a
	acolId : a matrix colId
	aval : a matrix value
	browId : csr form of rowId of matrix b
	bcolId : b matrix colId
	bval : b matrix value
	Rows : row # of matrices a or b
	rrowId : resultant matrix rowId
	rcolId : resultant matrix colId
	rval : resultant matrix value
	rleng : resultant matrix nnz */

	sparse_operation_t operation = SPARSE_OPERATION_NON_TRANSPOSE;
	sparse_matrix_t a, b, r;
	sparse_status_t s0;
	double alpha = 1.0;

	s0 = mkl_sparse_d_create_csr(&a, SPARSE_INDEX_BASE_ZERO, Rows, Rows, &arowId1[0], &arowId1[1], acolId, aval);
	//if (s0 == SPARSE_STATUS_SUCCESS)
	//	cout << "SPARSE_STATUS_SUCCESS\n";

	s0 = mkl_sparse_d_create_csr(&b, SPARSE_INDEX_BASE_ZERO, Rows, Rows, &browId1[0], &browId1[1], bcolId, bval);
	//if (s0 == SPARSE_STATUS_SUCCESS)
	//	cout << "SPARSE_STATUS_SUCCESS";

	s0 = mkl_sparse_d_add(operation, a, alpha, b, &r);   // add the two sparse matrices

	myint ARows, ACols;
	MKL_INT *ArowStart, *ArowEnd;
	MKL_INT *AcolId;
	double *Aval;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	s0 = mkl_sparse_d_export_csr(r, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
	sys->leng_Loo = ArowEnd[ARows - 1];

	checkCSRRepeats(ArowStart, ArowEnd, AcolId, Aval, sys->leng_Loo, sys->outside);

	sys->LooRowId = (myint*)malloc(sys->leng_Loo * sizeof(myint));
	sys->LooColId = (myint*)malloc(sys->leng_Loo * sizeof(myint));
	sys->Looval = (double*)malloc(sys->leng_Loo * sizeof(double));

	myint count = 0;
	for (myint ind = 0; ind < ARows; ++ind) {
		for (myint indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
			sys->LooRowId[count] = ind;
			sys->LooColId[count] = AcolId[indi];
			sys->Looval[count] = Aval[indi];
			count++;
		}
	}

	mkl_sparse_destroy(a);
	mkl_sparse_destroy(b);
	mkl_sparse_destroy(r);
}

int checkCSRRepeats(myint* RowIdStart, myint* RowIdEnd, myint* ColId, double* val, myint& leng, myint N) {
	/* check the CSR matrix format repeat elements and remove the 0 elements
	This matrix should be rowwise
	RowIdStart : matrix rowId start
	RowIdEnd : matrix rowId end, RowIdStart and RowIdEnd are actually the same array
	ColId : matrix colId
	val : matrix value
	leng : matrix nnz
	N : size of the matrix */
	myint i = 0, j = 0;   // i denotes the old array checking index, j denotes the new array index
	myint row = 0;

	while (row < N) {
		map<myint, double> col_val;   // the column number and its value and sort the column number
		myint start = j;   // new row start
		myint end;   // new row end
		while (i < RowIdEnd[row]) {
			if (col_val.find(ColId[i]) == col_val.end() && abs(val[i]) != 0) {   // this column hasn't shown up and this element's value is not zero
				col_val[ColId[i]] = val[i];
			}
			else if (col_val.find(ColId[i]) != col_val.end() && abs(val[i]) != 0) {   // this column has shown up
				col_val[ColId[i]] += val[i];
				if (col_val[ColId[i]] == 0)
					col_val.erase(ColId[i]);
			}
			i++;
		}
		end = start + col_val.size();
		i = RowIdEnd[row];    // the next row's start position
		RowIdStart[row] = start;
		RowIdEnd[row] = end;
		for (auto ci : col_val) {
			ColId[j] = ci.first;
			val[j] = ci.second;
			j++;
		}
		row++;
		col_val.clear();
	}
	leng = j;

	return 0;
}

void addOmegaEpsilon(fdtdMesh* sys, myint* RowId, myint* ColId, double* val, myint leng, myint N, double freq, double* val1) {
	/* generate the matrix M-omega^2*D_epsoo
	sys : used to provide D_eps
	RowId : matrix rowId
	ColId : matrix colId
	val : matrix value
	leng : nnz
	N : size of the matrix
	freq : frequency considered */
	for (int ind = 0; ind < leng; ++ind) {
		val1[ind] = val[ind];
		if (RowId[ind] == ColId[ind]) {
			val1[ind] -= pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[sys->mapioR[RowId[ind]]]);
		}
	}
}

void solveV0(double freq, double* rhs, complex<double>* u0, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, myint leng_v0d, sparse_matrix_t& v0ct, sparse_matrix_t& v0cat, myint leng_v0c, myint nedge, myint* AdRowId, myint* AdColId, double* Adval, myint leng_Ad, myint* AcRowId, myint* AcColId, double* Acval, myint leng_Ac, char ri) {
	/* solve (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0) with real right hand side
	sys : provide V0d, V0da, V0c, V0ca
	freq : frequency
	rhs : the real right hand side
	u0 : the solution */
	int status;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;

	double* temp = (double*)calloc(leng_v0d, sizeof(double));
	double* temp1 = (double*)calloc(nedge, sizeof(double));
	double* temp2 = (double*)calloc(leng_v0c, sizeof(double));
	double* u0c = (double*)calloc(leng_v0c, sizeof(double));
	double* u0di = (double*)calloc(leng_v0d, sizeof(double));
	double* u0dr = (double*)calloc(leng_v0d, sizeof(double));
	double* rhs1 = (double*)calloc(leng_v0d, sizeof(double));
	for (int indi = 0; indi < leng_v0d; ++indi) {   // Ad is not normalized with V0da and V0d
		rhs1[indi] = rhs[indi]; // sys->v0dan[indi] *
	}
	status = hypreSolve(AdRowId, AdColId, Adval, leng_Ad, rhs1, leng_v0d, temp, 1, 3);   // temp = (V0da'*D_eps*V0d)\bd
	//for (int ind = 0; ind < leng_v0d; ++ind) {
	//	temp[ind] *= sys->v0dn[ind];
	//}
	//
	//for (int ind = 0; ind < leng_v0d; ++ind) {
	//	temp[ind] /= sys->v0dn[ind];
	//}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, v0dt, descr, temp, beta, temp1);   // temp1 = V0d*(V0da'*D_eps*V0d)\bd
	for (int ind = 0; ind < nedge; ++ind) {
		temp1[ind] *= Doosval[ind];   // temp1 = D_eps*V0d*(V0da'*D_eps*V0d)\bd
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, v0cat, descr, temp1, beta, temp2);   // temp2 = V0ca'*D_eps*V0d*(V0da'*D_eps*V0d)\bd
	for (int ind = 0; ind < leng_v0c; ++ind) {
		//temp2[ind] /= sys->v0can[ind];
		temp2[ind] = rhs[leng_v0d + ind] - temp2[ind];   // temp2 = bc-V0ca'*D_eps*V0d*(V0da'*D_eps*V0d)\bd
	}

	//for (int indi = 0; indi < leng_v0c; ++indi) {   // Ac is not normalized with V0ca and V0c
	//	temp2[indi] *= sys->v0can[indi];
	//}
	status = hypreSolve(AcRowId, AcColId, Acval, leng_Ac, temp2, leng_v0c, u0c, 1, 3);
	for (int indi = 0; indi < leng_v0c; ++indi) {
		u0c[indi] /= -((freq * 2 * M_PI));// / sys->v0cn[indi]);    // u0c = (V0ca'*(omega*D_sig)*V0c)(bc-V0ca'*D_eps*V0d*(V0da'*D_eps*V0d)\bd), imaginary
		temp2[indi] = u0c[indi];
	}

	//for (int ind = 0; ind < leng_v0c; ++ind) {
	//	temp2[ind] /= sys->v0cn[ind];
	//}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, v0ct, descr, temp2, beta, temp1);   // temp1 = V0c*u0c, imaginary
	for (int ind = 0; ind < nedge; ++ind) {
		temp1[ind] *= Doosval[ind];   // temp1 = D_eps*V0c*u0c
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, v0dat, descr, temp1, beta, temp);
	//for (int ind = 0; ind < leng_v0d; ++ind) {
	//	temp[ind] /= sys->v0dan[ind];
	//}
	status = hypreSolve(AdRowId, AdColId, Adval, leng_Ad, temp, leng_v0d, u0di, 1, 3);   // u0di = (V0da'*D_eps*V0d)\(V0da'*D_eps*V0c*u0c), imaginary
	for (int ind = 0; ind < leng_v0d; ++ind) {
		u0di[ind] *= -1;// *sys->v0dn[ind];   // u0di = -(V0da'*D_eps*V0d)\(V0da'*D_eps*V0c*u0c), imaginary
	}
	status = hypreSolve(AdRowId, AdColId, Adval, leng_Ad, rhs1, leng_v0d, u0dr, 1, 3);   // u0dr = (V0da'*D_eps*V0d)\(bd), real
	for (int ind = 0; ind < leng_v0d; ++ind) {
		u0dr[ind] /= ((-pow(freq * 2 * M_PI, 2)));// / sys->v0dn[ind]);   // u0dr = (V0da'*(-omega^2*D_eps)*V0d)\(bd), real
	}

	/* Assign to the final solution */
	if (ri == 'r') {
		for (int ind = 0; ind < leng_v0d; ++ind) {
			u0[ind] = u0dr[ind] + 1i * u0di[ind];
		}
		for (int ind = 0; ind < leng_v0c; ++ind) {
			u0[leng_v0d + ind] = 0 + 1i * u0c[ind];
		}
	}
	else if (ri == 'i') {
		for (int ind = 0; ind < leng_v0d; ++ind) {
			u0[ind] = 1i * u0dr[ind] - u0di[ind];
		}
		for (int ind = 0; ind < leng_v0c; ++ind) {
			u0[leng_v0d + ind] = -u0c[ind];
		}
	}

	free(temp); temp = NULL;
	free(temp1); temp1 = NULL;
	free(temp2); temp2 = NULL;
	free(u0c); u0c = NULL;
	free(u0di); u0di = NULL;
	free(u0dr); u0dr = NULL;
	free(rhs1); rhs1 = NULL;
}

void V0Multiply(sparse_matrix_t& V0dt, sparse_matrix_t& V0ct, myint lengv0d, myint lengv0c, myint row, complex<double>* u1, complex<double>* u2) {
	/* Do V0*u1 and put the result in u2 */
	double alpha = 1;
	double beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* tempr = (double*)calloc(lengv0d + lengv0c, sizeof(double));
	double* tempi = (double*)calloc(lengv0d + lengv0c, sizeof(double));
	double* tempdr = (double*)calloc(row, sizeof(double));
	double* tempdi = (double*)calloc(row, sizeof(double));
	for (int ind = 0; ind < lengv0d; ++ind) {
		tempr[ind] = u1[ind].real();// / sys->v0dn[ind];
		tempi[ind] = u1[ind].imag();// / sys->v0dn[ind];
	}
	for (int ind = lengv0d; ind < lengv0d + lengv0c; ++ind) {
		tempr[ind] = u1[ind].real();// / sys->v0cn[ind - lengv0d];
		tempi[ind] = u1[ind].imag();// / sys->v0cn[ind - lengv0d];
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, tempr, beta, tempdr);   // tempdr = V0d*u1d.real
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, tempi, beta, tempdi);   // tempdi = V0d*u1d.imag

	double* tempcr = (double*)calloc(row, sizeof(double));
	double* tempci = (double*)calloc(row, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ct, descr, &tempr[lengv0d], beta, tempcr);   // tempcr = V0c*u1c.real
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ct, descr, &tempi[lengv0d], beta, tempci);   // tempci = V0c*u1c.imag

	for (int ind = 0; ind < row; ++ind) {
		u2[ind] = tempdr[ind] + tempcr[ind] + 1i * (tempdi[ind] + tempci[ind]);
	}

	free(tempr); tempr = NULL;
	free(tempi); tempi = NULL;
	free(tempdr); tempdr = NULL;
	free(tempdi); tempdi = NULL;
	free(tempcr); tempcr = NULL;
	free(tempci); tempci = NULL;
}

void generateL1(fdtdMesh* sys, double freq, myint* LRowId, myint* LColId, double* Lval, myint leng_L, myint N, myint* L1RowId, myint* L1ColId, double* L1val, myint& leng_L1) {
	myint i = 0, start, ind = 0;

	while (i < leng_L) {
		start = LRowId[i];

		while (i < leng_L && LRowId[i] == start) {
			L1RowId[ind] = start;
			L1ColId[ind] = LColId[i];
			L1val[ind] = Lval[i];
			if (LRowId[i] == LColId[i]) {
				L1val[ind] -= pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[LRowId[i]]);
			}
			ind++;
			i++;
		}
		if (sys->markEdge[sys->mapEdgeR[start]]) {
			L1RowId[ind] = start;
			L1ColId[ind] = N + start;
			L1val[ind] = -freq * 2 * M_PI * SIGMA;
			ind++;
		}
	}
	i = 0;
	while (i < leng_L) {
		start = LRowId[i];
		if (sys->markEdge[sys->mapEdgeR[start]]) {
			L1RowId[ind] = N + start;
			L1ColId[ind] = start;
			L1val[ind] = freq * 2 * M_PI * SIGMA;
			ind++;
		}
		while (i < leng_L && LRowId[i] == start) {
			L1RowId[ind] = N + start;
			L1ColId[ind] = N + LColId[i];
			L1val[ind] = Lval[i];
			if (LRowId[i] == LColId[i]) {
				L1val[ind] -= pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[LRowId[i]]);
			}
			i++;
			ind++;
		}

	}
	leng_L1 = ind;
}

void V0aMultiply(fdtdMesh* sys, sparse_matrix_t& V0dat, sparse_matrix_t& V0cat, myint lengv0d, myint lengv0c, myint nedge, complex<double>* u2, complex<double>* u1) {
	double alpha = 1;
	double beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* tempr = (double*)calloc(nedge, sizeof(double));
	double* tempi = (double*)calloc(nedge, sizeof(double));
	double* tempdr = (double*)calloc(lengv0d, sizeof(double));
	double* tempdi = (double*)calloc(lengv0d, sizeof(double));
	double* tempcr = (double*)calloc(lengv0c, sizeof(double));
	double* tempci = (double*)calloc(lengv0c, sizeof(double));
	for (int ind = 0; ind < nedge; ++ind) {
		tempr[ind] = u2[ind].real();
		tempi[ind] = u2[ind].imag();
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, tempr, beta, tempdr);
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, tempi, beta, tempdi);
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, tempr, beta, tempcr);
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, tempi, beta, tempci);
	for (int ind = 0; ind < lengv0d; ++ind) {
		tempdr[ind] /= sys->v0dan[ind];
		tempdi[ind] /= sys->v0dan[ind];
		u1[ind] = tempdr[ind] + 1i * tempdi[ind];
	}
	for (int ind = 0; ind < lengv0c; ++ind) {
		tempcr[ind] /= sys->v0can[ind];
		tempci[ind] /= sys->v0can[ind];
		u1[ind + lengv0d] = tempcr[ind] + 1i * tempci[ind];
	}
	free(tempr); tempr = NULL;
	free(tempi); tempi = NULL;
	free(tempdr); tempdr = NULL;
	free(tempdi); tempdi = NULL;
	free(tempcr); tempcr = NULL;
	free(tempci); tempci = NULL;
}

void generateP(fdtdMesh* sys, double freq, myint* LRowId, myint* LColId, double* Lval, myint leng_L, myint N, myint* PRowId, myint* PColId, double* Pval, myint& leng_P) {
	/* Generate matrix [-omega^2*D_eps+L,0;
	0, -omega^2*D_eps+L]
	PRowId : the generated matrix P's rowId
	PColId : the generated matrix P's colId
	Pval : the generated matrix P's value */

	myint i = 0, ind = 0;

	while (i < leng_L) {
		PRowId[ind] = LRowId[i];
		PColId[ind] = LColId[i];
		Pval[ind] = Lval[i];
		if (LRowId[i] == LColId[i]) {
			Pval[ind] -= pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[LRowId[i]]);
		}
		i++;
		ind++;
	}
	i = 0;
	while (i < leng_L) {
		PRowId[ind] = N + LRowId[i];
		PColId[ind] = N + LColId[i];
		Pval[ind] = Lval[i];
		if (LRowId[i] == LColId[i]) {
			Pval[ind] -= pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[LRowId[i]]);
		}
		i++;
		ind++;
	}
	leng_P = ind;


}

void addDiagV0bV0ba(fdtdMesh* sys, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo) {
	/* return Loo + diag(V0b*V0ba'/mu)
	LooRowId : Loo row Id
	LooColId : Loo col Id
	Looval : Loo value
	leng_Loo : nnz in Loo */
	myint node1, node2;
	double l, la;

	for (int ind = 0; ind < leng_Loo; ++ind) {
		if (LooRowId[ind] == LooColId[ind]) {   // the diagonal
			sys->compute_edgelink(sys->mapEdgeR[sys->mapioR[LooRowId[ind]]], node1, node2);
			//if (sys->markNode[node1] != 0 && sys->markNode[node2] != 0)   // this edge is the shared edge
			//	continue;
			if (sys->markNode[node1] != 0) {   // node1 is on the conductor boundary
				sys->nodeLength(node1, node2, 1, l, la);
				Looval[ind] -= (1 / l) * (1 / la) / MU;
			}
			if (sys->markNode[node2] != 0) {   // node2 is on the conductor boundary
				sys->nodeLength(node1, node2, 2, l, la);
				Looval[ind] -= (1 / l) * (1 / la) / MU;
			}
		}
	}
}

int gmres_V0dLoo(double* Doosval, sparse_matrix_t& V0dt, double* v0dn, sparse_matrix_t& V0dat, double* v0dan, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* A22dbRowId, myint* A22dbColId, double* A22dbval, myint leng_A22db, double* bm, double* res, myint N, myint N1, double freq) {
	/* solve [V0dan'*Doos*V0dn, V0dan'*Doos;
	Doos*V0dn, Doos+Soos+V0d*V0da'/mu] with preconditioner
	[V0dan'*Doos*V0dn, V0dan'*Doos;
	0, Doos+Soos+V0d*V0da'/mu]
	Doosval : the diagonal value of Doos without D_epsoo
	V0dt : V0d's handler
	v0dn : norm of each V0d vector
	V0dat : V0da's handler
	v0dan : norm of each V0da vector
	A11RowId : V0dan'*Deps*V0dn's row ID
	A11ColId : V0dan'*Deps*V0dn's col ID
	A11val : V0dan'*Deps*V0dn's value normalized
	leng_A11 : nnz in V0da'*Deps*V0d
	A22RowId : A22's row ID
	A22ColId : A22's col ID
	A22val : A22's value
	leng_A22 : nnz in A22
	A22dbRowId : A22 preconditioner row ID
	A22dbColId : A22 preconditioner col ID
	A22dbval : A22 preconditioner value
	leng_A22db : nn in A22 preconditioner
	bm : right hand side
	res : the solution after solving this
	N : the matrix size */

	ofstream out;
	myint size = 128;
	myint N2 = N - N1;
	/*------------------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/*------------------------------------------------------------------------------------*/
	MKL_INT ipar[size];
	ipar[14] = 500;    // restart number
	cout << "The size of the problem is " << N << endl;
	//double b[N];
	//double expected_solution[N];
	//double computed_solution[N];
	//double residual[N];
	double* b = (double*)malloc(N * sizeof(double));
	double* expected_solution = (double*)malloc(N * sizeof(double));
	double* computed_solution = (double*)malloc(N * sizeof(double));
	double* residual = (double*)malloc(N * sizeof(double));
	double dpar[size];
	double* tmp = (double*)calloc(N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1, sizeof(double));
	/*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	ivar = N;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double dvar;
	int status, maxit = 500;
	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector b for future use
	/*---------------------------------------------------------------------------*/
	i = 0;
	double bmn = 0;
	for (i = 0; i < N; ++i) {
		b[i] = bm[i];
	}
	/*--------------------------------------------------------------------------
	/* Initialize the initial guess
	/*--------------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		computed_solution[i] = 0;// bm[i];
	}


	/*--------------------------------------------------------------------------
	/* Initialize the solver
	/*--------------------------------------------------------------------------*/
	dfgmres_init(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
	ipar[10] = 1;   // use preconditioner
	ipar[7] = 0;
	dpar[0] = 1.0E-3;

	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
ONE: dfgmres(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	//bmn = 0;
	//for (int ind = 0; ind < ivar; ++ind) {
	//	bmn += computed_solution[ind] * computed_solution[ind];
	//}
	//cout << "Right hand side norm is " << bmn << endl;
	if (RCI_request == 0) {
		goto COMPLETE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------*/
	if (RCI_request == 1) {
		/* A11*v1 = V0dan'*Doos*V0dn*v1 and A21*v1 */
		AtimesV(V0dt, v0dn, V0dat, v0dan, Doosval, freq, N1, N2, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1], A22RowId, A22ColId, A22val, leng_A22);

		goto ONE;
	}

	/* do the user-defined stopping test */
	if (RCI_request == 2) {
		ipar[12] = 1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */

		for (myint ind = 0; ind < N; ++ind) {
			residual[ind] = 0;   // before using sparseMatrixVecMul, the resultant vector should be first initialized
		}
		AtimesV(V0dt, v0dn, V0dat, v0dan, Doosval, freq, N1, N2, b, residual, A22RowId, A22ColId, A22val, leng_A22);
		dvar = -1.0E0;
		i = 1;
		daxpy(&ivar, &dvar, bm, &i, residual, &i);
		dvar = cblas_dnrm2(ivar, residual, i) / cblas_dnrm2(ivar, bm, i);    // relative residual
		cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
		if (dvar < 0.00001 || itercount > maxit) goto COMPLETE;
		else goto ONE;
	}

	/* apply the preconditioner on the vector tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1] */
	if (RCI_request == 3) {
		/* Apply [A11, A12;
		0, A22] as preconditioner */
		for (int ind = 0; ind < N; ++ind) {
			tmp[ipar[22] - 1 + ind] = 0;
		}
		PrecondUpperTri(freq, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22ColId, A22val, leng_A22, A22dbRowId, A22dbColId, A22dbval, leng_A22db, V0dt, v0dn, V0dat, v0dan, Doosval, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1], N1, N2);
		goto ONE;
	}

	/* check if the norm of the next generated vector is not zero up to rounding and computational errors. */
	if (RCI_request == 4) {
		if (dpar[6] < 1.0E-12) goto COMPLETE;
		else goto ONE;
		//goto ONE;
	}

	else {
		goto FAILED;
	}

COMPLETE: ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
	for (i = 0; i < ivar; ++i) {
		res[i] = computed_solution[i];
	}
	goto SUCCEDED;
FAILED: cout << "The solver has returned the ERROR code " << RCI_request << endl;

SUCCEDED: free(b); b = NULL;
	free(expected_solution); expected_solution = NULL;
	free(computed_solution); computed_solution = NULL;
	free(residual); residual = NULL;
	return 0;

}

void gmres_V0dV0bLoo(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* bm, double* res, void* pt[64], myint iparm[64]) {
	/* gmres to solve 3*3 system
	[v0dan'*Doos*v0dn, v0dan'*Doos*v0bn, v0dan'*Doos;
	v0ban'*Doos*v0dn, v0ban'*(Doos+Soos)*v0bn, v0ban'*(Doos+Soos);
	Doos*V0dn, (Doos+Soos)*v0dn, Doos+Loos] with upper triangular matrix
	frefrequency
	Doosval : D_epsoo
	V0ddt : V0dnq :
	V0ddat : V0dan
	N1 : column number of V0d
	V0bt : V0bn
	V0bat : V0ban
	N2 : column number of V0b
	Soo : Soo+Doos
	N3 : matrix col or row # of Soo
	LooRowId : row ID of Loos+Doos
	LooColId : col ID of Loos+Doos
	Looval : value of Loos+Doos
	leng_Loo : nnz of Loos+Doos
	bm : righ hand side
	res : solution result
	*/

	ofstream out;
	out.open("iteration_error.txt", ofstream::out | ofstream::trunc);
	myint size = 128;
	myint N = N1 + N2 + N3;
	/*------------------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/*------------------------------------------------------------------------------------*/
	MKL_INT ipar[size];
	ipar[14] = 500;    // restart number
	cout << "The size of the problem is " << N << endl;
	//double b[N];
	//double expected_solution[N];
	//double computed_solution[N];
	//double residual[N];
	double* b = (double*)malloc(N * sizeof(double));
	double* expected_solution = (double*)malloc(N * sizeof(double));
	double* computed_solution = (double*)malloc(N * sizeof(double));
	double* residual = (double*)malloc(N * sizeof(double));
	double dpar[size];
	double* tmp = (double*)malloc((N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1) * sizeof(double));
	cout << "The memory of tmp is " << N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1 << endl;
	/*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	ivar = N;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double dvar;
	int status, maxit = 500;
	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector b for future use
	/*---------------------------------------------------------------------------*/
	i = 0;
	double bmn = 0;
	for (i = 0; i < N; ++i) {
		b[i] = bm[i];
	}
	/*--------------------------------------------------------------------------
	/* Initialize the initial guess
	/*--------------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		computed_solution[i] = 0;// bm[i];
	}

	/*--------------------------------------------------------------------------
	/* Initialize the solver
	/*--------------------------------------------------------------------------*/

	dfgmres_init(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
	dpar[0] = 0.01;
	ipar[10] = 1;   // use preconditioner
	ipar[7] = 0;
	ipar[14] = 500;   // restart number

					  /*---------------------------------------------------------------------------
					  /* Check the correctness and consistency of the newly set parameters
					  /*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
ONE: dfgmres(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	//bmn = 0;
	//for (int ind = 0; ind < ivar; ++ind) {
	//	bmn += computed_solution[ind] * computed_solution[ind];
	//}
	//cout << "Right hand side norm is " << bmn << endl;
	if (RCI_request == 0) {
		goto COMPLETE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------*/
	if (RCI_request == 1) {
		/* [v0dan'*Doos*v0dn, v0dan'*Doos*v0bn, v0dan'*Doos;
		v0ban'*Doos*v0dn, v0ban'*(Doos+Soos)*v0bn, v0ban'*(Doos+Soos);
		Doos*V0dn, (Doos+Soos)*v0dn, Doos+Loos]*v */
		//out.open("x1.txt", ofstream::out | ofstream::trunc);
		//for (int ind = 0; ind < N; ++ind) {
		//	out << setprecision(15) << tmp[ipar[21] - 1 + ind] << endl;
		//}
		//out.close();

		A3b3timesV(Doosval, freq, V0ddt, V0ddat, N1, V0bt, V0bat, N2, Soo, N3, LooRowId, LooColId, Looval, leng_Loo, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
		//out.open("x2.txt", ofstream::out | ofstream::trunc);
		//for (int ind = 0; ind < N; ++ind) {
		//	out << setprecision(15) << tmp[ipar[22] - 1 + ind] << endl;
		//}
		//out.close();

		goto ONE;
	}

	/* do the user-defined stopping test */
	if (RCI_request == 2) {
		ipar[12] = 1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */
		//out.open("x1.txt", ofstream::out | ofstream::trunc);
		//for (int ind = 0; ind < N; ++ind) {
		//	out << setprecision(15) << b[ind] << endl;
		//}
		//out.close();

		//for (myint ind = 0; ind < N; ++ind) {
		//	residual[ind] = 0;
		//}

		A3b3timesV(Doosval, freq, V0ddt, V0ddat, N1, V0bt, V0bat, N2, Soo, N3, LooRowId, LooColId, Looval, leng_Loo, b, residual);
		//out.open("x2.txt", ofstream::out | ofstream::trunc);
		//for (int ind = 0; ind < N; ++ind) {
		//	out << setprecision(15) << residual[ind] << endl;
		//}
		//out.close();

		dvar = -1.0E0;
		i = 1;
		daxpy(&ivar, &dvar, bm, &i, residual, &i);
		dvar = cblas_dnrm2(ivar, residual, i) / cblas_dnrm2(ivar, bm, i);    // relative residual
		out << itercount << " " << dvar << endl;
		if (dvar < dpar[0] || itercount > maxit) goto COMPLETE;
		else goto ONE;
	}

	/* apply the preconditioner on the vector tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1] */
	if (RCI_request == 3) {
		/* Apply [A11, A12;
		0, A22] as preconditioner */

		//for (int ind = 0; ind < N; ++ind) {
		//	tmp[ipar[22] - 1 + ind] = 0;
		//}
		//out.open("x1.txt", ofstream::out | ofstream::trunc);
		//for (int ind = 0; ind < N; ++ind) {
		//	out << setprecision(15) << tmp[ipar[21] - 1 + ind] << endl;
		//}
		//out.close();

		Precond3b3UpperTri(freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1], pt, iparm);   // apply the upper triangular matrix as the preconditioner

		goto ONE;
	}

	/* check if the norm of the next generated vector is not zero up to rounding and computational errors. */
	if (RCI_request == 4) {
		if (dpar[6] < 1.0E-12) goto COMPLETE;
		else goto ONE;
		//goto ONE;
	}

	else {
		goto FAILED;
	}

COMPLETE: ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
	for (i = 0; i < ivar; ++i) {
		res[i] = computed_solution[i];
	}
	MKL_Free_Buffers();
	goto SUCCEDED;
FAILED: cout << "The solver has returned the ERROR code " << RCI_request << endl;
	MKL_Free_Buffers();

SUCCEDED: free(b); b = NULL;
	free(tmp); tmp = NULL;
	free(expected_solution); expected_solution = NULL;
	free(computed_solution); computed_solution = NULL;
	free(residual); residual = NULL;
	out.close();
}

void AtimesV(sparse_matrix_t& V0dt, double* v0dn, sparse_matrix_t& V0dat, double* v0dan, double* Doosval, double freq, myint N1, myint N2, double* V1, double* V2, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22) {
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* A21V1 = (double*)calloc(N2, sizeof(double));
	double* tmp_cp = (double*)calloc(N1, sizeof(double));
	for (int ind = 0; ind < N1; ++ind) {
		tmp_cp[ind] = V1[ind] / v0dn[ind];
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, tmp_cp, beta, A21V1);   // A21V1 = V0dn*V1
	for (int ind = 0; ind < N2; ++ind) {
		A21V1[ind] *= Doosval[ind] * pow(freq * 2 * M_PI, 2) * (-1);   // A21V1 = Doos*V0dn*V1
	}
	double* A11V1 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, A21V1, beta, A11V1);   // A11V1 = A11*V1
	for (int ind = 0; ind < N1; ++ind) {
		A11V1[ind] /= v0dan[ind];
	}

	/* A12*v2 = V0dan'*Doos*v2 */
	double* DoosV2 = (double*)calloc(N2, sizeof(double));
	for (int ind = 0; ind < N2; ++ind) {
		DoosV2[ind] = V1[N1 + ind] * Doosval[ind] * pow(freq * 2 * M_PI, 2) * (-1);   // DoosV2 = Doos*V2
	}
	double* A12V2 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, DoosV2, beta, A12V2);   // A12V2 = A12*V2
	for (int ind = 0; ind < N1; ++ind) {
		A12V2[ind] /= v0dan[ind];
	}

	/* A22*v2 */
	double* A22V2 = (double*)calloc(N2, sizeof(double));
	sparseMatrixVecMul(A22RowId, A22ColId, A22val, leng_A22, &V1[N1], A22V2);

	for (int ind = 0; ind < N1; ++ind) {
		V2[ind] = A11V1[ind] + A12V2[ind];
	}
	for (int ind = 0; ind < N2; ++ind) {
		V2[N1 + ind] = A21V1[ind] + A22V2[ind];
	}

	free(A21V1); A21V1 = NULL;
	free(tmp_cp); tmp_cp = NULL;
	free(A11V1); A11V1 = NULL;
	free(DoosV2); DoosV2 = NULL;
	free(A12V2); A12V2 = NULL;
	free(A22V2); A22V2 = NULL;

}

void A3b3timesV(double* Doosval, double freq, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* v, double* res) {
	/* Do [v0dan'*Doos*v0dn, v0dan'*Doos*v0bn, v0dan'*Doos;
	v0ban'*Doos*v0dn, v0ban'*(Doos+Soos)*v0bn, v0ban'*(Doos+Soos);
	Doos*V0dn, (Doos+Soos)*v0dn, Doos+Loos] multiply a vector */
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	ofstream out;


	/* first line */
	double* v0dnV1 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, v, beta, v0dnV1);   // v0dnV1 = V0dn*v1
	double* v0bnV2 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &v[N1], beta, v0bnV2);   // v0bnV2 = V0bn*v2

	double* temp = (double*)malloc(N3 * sizeof(double));
	for (int ind = 0; ind < N3; ++ind) {
		temp[ind] = v0dnV1[ind] + v0bnV2[ind] + v[N1 + N2 + ind];   // temp = V0dn*v1 + V0bn*v2 + v3
		temp[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // temp = (-omega^2*D_epsoo)*(V0dn*v1 + V0bn*v2 + v3)
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, temp, beta, res);   // res1 = V0dan'*(-omega^2*D_epsoo)*(V0dn*v1 + V0bn*v2 + v3)

																								  /* second line */
	double* DoosV0dnV1 = (double*)malloc(N3 * sizeof(double));
	for (int ind = 0; ind < N3; ++ind) {
		DoosV0dnV1[ind] = (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind] * v0dnV1[ind];    // DoosV0dnV1 = -omega^2*D_epsoo*V0dn*v1
	}
	double* SoosV0bnV2 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, v0bnV2, beta, SoosV0bnV2);   // SoosV0bnV2 = (-omega^2*D_epsoo+Soos)*V0bn*V2
	double* SoosV3 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, &v[N1 + N2], beta, SoosV3);   // SoosV3 = (-omega^2*D_epsoo+Soos)*V3
	for (int ind = 0; ind < N3; ++ind) {
		temp[ind] = DoosV0dnV1[ind] + SoosV0bnV2[ind] + SoosV3[ind];    // temp = (-omega^2*D_epsoo)*V0dn*V1+(-omega^2*D_epsoo+Soos)*(V0bn*V2+V3)
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, temp, beta, &res[N1]);    // res2 = V0ban'*((-omega^2*D_epsoo)*V0dn*V1+(-omega^2*D_epsoo+Soos)*(V0bn*V2+V3))

																									   /* third line */
	double* LoosV3 = (double*)calloc(N3, sizeof(double));   // use sparseMatrixVecMul the initial vector should be zero
	sparseMatrixVecMul(LooRowId, LooColId, Looval, leng_Loo, &v[N1 + N2], LoosV3);    // LoosV3 = (-omega^2*D_epsoo+Loos)*V3
	for (int ind = 0; ind < N3; ++ind) {
		res[N1 + N2 + ind] = DoosV0dnV1[ind] + SoosV0bnV2[ind] + LoosV3[ind];
	}

	free(v0dnV1); v0dnV1 = NULL;
	free(v0bnV2); v0bnV2 = NULL;
	free(temp); temp = NULL;
	free(DoosV0dnV1); DoosV0dnV1 = NULL;
	free(SoosV0bnV2); SoosV0bnV2 = NULL;
	free(SoosV3); SoosV3 = NULL;
	free(LoosV3); LoosV3 = NULL;

}

void PrecondUpperTri(double freq, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* A22dbRowId, myint* A22dbColId, double* A22dbval, myint leng_A22db, sparse_matrix_t& V0dt, double* v0dn, sparse_matrix_t& V0dat, double* v0dan, double* Doosval, double* b, double* x, myint N1, myint N2) {
	/* Apply the upper triangular matrix as the preconditioner
	[A11, A12;
	0, A22]
	freq : frequency considered
	A11RowId : V0dan'*Deps*V0dn's row ID, without omega
	A11ColId : V0dan'*Deps*V0dn's col ID, without omega
	A11val : V0dan'*Deps*V0dn's value, without omega
	leng_A11 : nnz in A11
	A22RowId : A22 row ID
	A22ColId : A22 col ID
	A22val : A22 value
	leng_A22 : nnz in A22
	A22dbRowId : A22 preconditioner's row ID
	A22dbColId : A22 preconditioner's col ID
	A22dbval : A22 preconditioner's value
	leng_A22db : nnz in A22 preconditioner
	Doosval : Doos diagonal matrix, with no omega D_epsoo
	b : rhs
	x : results */

	int status;
	ofstream out;

	/* x2 = A22^(-1)b2 */
	//gmres_solveA22(A22RowId, A22ColId, A22val, leng_A22, A22dbRowId, A22dbColId, A22dbval, leng_A22db, &b[N1], &x[N1], N2);
	//status = hypreSolve(A22RowId, A22ColId, A22val, leng_A22, &b[N1], N2, &x[N1], 1, 3);   // x = A22^(-1)b2
	status = hypreSolve(A22dbRowId, A22dbColId, A22dbval, leng_A22db, &b[N1], N2, &x[N1], 0, 3);   // x = A22^(-1)b2
																								   //myint* A22dbRowId1 = (myint*)calloc((N2 + 1), sizeof(myint));
																								   //status = COO2CSR_malloc(A22RowId, A22ColId, A22val, leng_A22, N2, A22dbRowId1);
																								   //status = pardisoSolve(A22dbRowId1, A22ColId, A22val, &b[N1], &x[N1], N2);

																								   /* x1 = A11^(-1)(b1-A12*x2) */
	double* Doosx2 = (double*)calloc(N2, sizeof(double));
	for (int ind = 0; ind < N2; ++ind) {
		Doosx2[ind] = Doosval[ind] * x[N1 + ind] * (-1) * pow(freq * 2 * M_PI, 2);   // Doosx2 = Doos*x2
	}
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* A12x2 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, Doosx2, beta, A12x2);   // A12x2 = V0da'*Doos*x2
	for (int ind = 0; ind < N1; ++ind) {
		A12x2[ind] /= v0dan[ind];   // A12x2 = A12*x2
		A12x2[ind] = b[ind] - A12x2[ind];   // A12x2 = b1 - A12*x2
	}
	status = hypreSolve(A11RowId, A11ColId, A11val, leng_A11, A12x2, N1, x, 0, 3);   // x = A11^(-1)(b1-A12*x2)
	for (int ind = 0; ind < N1; ++ind) {
		x[ind] /= (-1) * pow(2 * M_PI * freq, 2);   // x = A11^(-1)(b1-A12*x2)
	}

	free(Doosx2); Doosx2 = NULL;
	free(A12x2); A12x2 = NULL;
	//free(A22dbRowId1); A22dbRowId1 = NULL;
}

void Precond3b3UpperTri(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* rhs, double* x, void* pt[64], myint iparm[64]) {
	/* Apply the 3*3 upper triangular matrix as the preconditioner
	[v0dan'*Doos*v0dn, v0dan'*Doos*v0bn, v0dan'*Doos;
	0, v0ban'*(Doos+Soos)*v0bn, v0ban'*(Doos+Soos);
	0, 0, Doos+Loos]
	*/
	ofstream out;
	int status;

	/* calculate x3 */
	//clock_t t = clock();
	status = hypreSolve(LooRowId, LooColId, Looval, leng_Loo, &rhs[N1 + N2], N3, &x[N1 + N2], 1, 3);   // x = (-omega^2*D_epsoo+Loos)^(-1)*b3
																									   /* pardiso sanity check */
																									   //myint* LooRowId1 = (myint*)calloc((N3 + 1), sizeof(myint));
																									   //status = COO2CSR_malloc(LooRowId, LooColId, Looval, leng_Loo, N3, LooRowId1);
																									   //status = pardisoSolve(LooRowId1, LooColId, Looval, &rhs[N1 + N2], &x[N1 + N2], N3);
																									   //mkl_gmres_A(&rhs[N1 + N2], &x[N1 + N2], LooRowId, LooColId, Looval, leng_Loo, N3);
																									   //cout << "Time for solve A33 is " << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;

																									   /* calculate x2 */
	double* Soosx3 = (double*)calloc(N3, sizeof(double));
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, &x[N1 + N2], beta, Soosx3);   // Soosx3 = (-omega^2*Depsoo+Soos)*x3
	double* v0banSoosx3 = (double*)calloc(N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, Soosx3, beta, v0banSoosx3);   // v0banSoosx3 = V0ban'*(-omega^2*D_epsoo+Soos)*x3
	for (int ind = 0; ind < N2; ++ind) {
		v0banSoosx3[ind] = rhs[N1 + ind] - v0banSoosx3[ind];   // v0banSoosx2 = rhs2-V0ban'*(-omega^2*D_epsoo+Soos)*x3
	}
	//t = clock();
	//status = pardisoSolve(A22RowId1, A22ColId, A22val, v0banSoosx3, &x[N1], N2);   // x2 = A22^(-1)*(rhs2-V0ban'*(-omega^2*D_epsoo+Soos)*x3)
	pardisoSol(A22RowId1, A22ColId, A22val, pt, iparm, v0banSoosx3, &x[N1], N2);
	//status = hypreSolve(A22RowId, A22ColId, A22val, leng_A22, v0banSoosx3, N2, &x[N1], 1, 3);
	//mkl_gmres_A(v0banSoosx3, &x[N1], A22RowId, A22ColId, A22val, leng_A22, N2);
	//cout << "Time to solve A22 is " << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;

	/* calculate x1 */
	double* Doosx3 = (double*)calloc(N3, sizeof(double));
	double* v0bnx2 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &x[N1], beta, v0bnx2);   // v0bnx2 = V0bn*x2
	for (int ind = 0; ind < N3; ++ind) {
		Doosx3[ind] = (x[N1 + N2 + ind] + v0bnx2[ind]) * (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // Doosx3 = (-omega^2*D_epsoo)*(V0bn*x2+x3)
	}
	double* v0datx2x3 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, Doosx3, beta, v0datx2x3);    // v0datx2x3 = V0da'*(-omega^2*D_epsoo)*(V0b*x2+x3)
	for (int ind = 0; ind < N1; ++ind) {
		v0datx2x3[ind] = (rhs[ind] - v0datx2x3[ind]) * v0dan[ind];   // v0datx2x3 = rhs1-V0da'*(-omega^2*D_epsoo)*(V0b*x2+x3)
	}
	//t = clock();
	//mkl_gmres_A(v0datx2x3, x, A11RowId, A11ColId, A11val, leng_A11, N1);
	status = hypreSolve(A11RowId, A11ColId, A11val, leng_A11, v0datx2x3, N1, x, 1, 3);
	/* pardiso sanity check */
	//myint* A11RowId1 = (myint*)calloc((N1 + 1), sizeof(myint));
	//status = COO2CSR_malloc(A11RowId, A11ColId, A11val, leng_A11, N1, A11RowId1);
	//status = pardisoSolve(A11RowId1, A11ColId, A11val, v0datx2x3, x, N1);
	//cout << "Time to solve A11 is " << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;

	for (int ind = 0; ind < N1; ++ind) {
		x[ind] *= v0dn[ind];
		x[ind] /= ((-1) * pow(freq * 2 * M_PI, 2));   // x1 = A11^(-1)*(rhs1-V0da'*(-omega^2*D_epsoo)*(V0b*x2+x3))
	}

	//free(LooRowId1); LooRowId1 = NULL;
	//free(A11RowId1); A11RowId1 = NULL;

	free(Soosx3); Soosx3 = NULL;
	free(v0banSoosx3); v0banSoosx3 = NULL;
	free(Doosx3); Doosx3 = NULL;
	free(v0bnx2); v0bnx2 = NULL;
	free(v0datx2x3); v0datx2x3 = NULL;

}

void pardisoSol(myint* RowId1, myint* ColId, double* val, void* pt[64], myint iparm[64], double* b, double* x, myint N) {
	myint mtype = 11;    /* Real unsymmetric matrix */
	myint nrhs = 1;    /* Number of right hand sides */

					   /* Pardiso control parameters */
	myint maxfct, mnum, phase, error, msglvl, solver;
	double dparm[64];
	int v0csin;
	myint perm;
	myint size = N;

	/* Auxiliary variables */
	char *var;

	msglvl = 0;    /* print statistical information */
	solver = 0;    /* use sparse direct solver */
	error = 0;
	maxfct = 1;
	mnum = 1;
	phase = 33;

	pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, val, RowId1, ColId, &perm, &nrhs, iparm, &msglvl, b, x, &error);
	if (error != 0) {
		printf("\nERROR during numerical factorization: %d\n", error);
		exit(2);
	}

}

int gmres_solveA22(myint* A22RowId, myint* A22ColId, double* A22val, myint leng_A22, myint* A22dbRowId, myint* A22dbColId, double* A22dbval, myint leng_A22db, double* bm, double* x, myint N) {
	/* gmres to solve A22 = Doos + Soos - Ltoos + LtooC
	= Doos + Soos + V0b*V0ba'/mu
	with preconditioner A22db = Doos + Soos - Ltoos + diag(LtooC)
	= Doos + Soos + V0b*V0ba'/mu + V0b*V0ba'/mu - diag(V0b*V0ba'/mu)
	bm : right hand side
	x : results
	N : size of the matrix */

	myint size = 128;
	/*------------------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/*------------------------------------------------------------------------------------*/
	MKL_INT ipar[size];
	ipar[14] = 1000;
	//double b[N];
	//double expected_solution[N];
	//double computed_solution[N];
	//double residual[N];
	double* b = (double*)malloc(N * sizeof(double));
	double* expected_solution = (double*)malloc(N * sizeof(double));
	double* computed_solution = (double*)malloc(N * sizeof(double));
	double* residual = (double*)malloc(N * sizeof(double));
	double dpar[size];
	double* tmp = (double*)calloc(N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1, sizeof(double));
	/*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	ivar = N;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double dvar;
	int status, maxit = 1000;
	ofstream out;
	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector b for future use
	/*---------------------------------------------------------------------------*/
	i = 0;
	for (i = 0; i < N; ++i) {
		b[i] = bm[i];
	}
	/*--------------------------------------------------------------------------
	/* Initialize the initial guess
	/*--------------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		computed_solution[i] = 0;// bm[i];
	}

	/*--------------------------------------------------------------------------
	/* Initialize the solver
	/*--------------------------------------------------------------------------*/
	dfgmres_init(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
	ipar[10] = 1;   // use preconditioner
	ipar[14] = 500;   // restart number
	ipar[7] = 0;
	dpar[0] = 1.0E-2;
	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
ONE: dfgmres(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	//bmn = 0;
	//for (int ind = 0; ind < ivar; ++ind) {
	//	bmn += computed_solution[ind] * computed_solution[ind];
	//}
	//cout << "Right hand side norm is " << bmn << endl;
	if (RCI_request == 0) {
		goto COMPLETE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------*/
	if (RCI_request == 1) {
		/* A22*tmp[ipar[21] - 1] and put the result in tmp[ipar[22] - 1] */
		for (int ind = 0; ind < N; ++ind) {
			tmp[ipar[22] - 1 + ind] = 0;
		}
		sparseMatrixVecMul(A22RowId, A22ColId, A22val, leng_A22, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);

		goto ONE;
	}

	/* do the user-defined stopping test */
	if (RCI_request == 2) {
		ipar[12] = 1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */

		for (myint ind = 0; ind < N; ++ind) {
			residual[ind] = 0;   // before using sparseMatrixVecMul, the resultant vector should be first initialized
		}
		sparseMatrixVecMul(A22RowId, A22ColId, A22val, leng_A22, b, residual);   // A22 multiply x and put the result in residual
		dvar = -1.0E0;
		i = 1;
		daxpy(&ivar, &dvar, bm, &i, residual, &i);
		dvar = cblas_dnrm2(ivar, residual, i) / cblas_dnrm2(ivar, bm, i);    // relative residual
		cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
		if (dvar < dpar[0] || itercount > maxit) goto COMPLETE;
		else goto ONE;
	}

	/* apply the preconditioner on the vector tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1] */
	if (RCI_request == 3) {
		status = hypreSolve(A22dbRowId, A22dbColId, A22dbval, leng_A22db, &tmp[ipar[21] - 1], N, &tmp[ipar[22] - 1], 1, 3);
		goto ONE;
	}

	/* check if the norm of the next generated vector is not zero up to rounding and computational errors. */
	if (RCI_request == 4) {
		if (dpar[6] < 1.0E-12) goto COMPLETE;
		else goto ONE;
		//goto ONE;
	}

	else {
		goto FAILED;
	}

COMPLETE: ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
	for (i = 0; i < ivar; ++i) {
		x[i] = computed_solution[i];
	}
	goto SUCCEDED;
FAILED: cout << "The solver has returned the ERROR code " << RCI_request << endl;

SUCCEDED: free(b); b = NULL;
	free(expected_solution); expected_solution = NULL;
	free(computed_solution); computed_solution = NULL;
	free(residual); residual = NULL;
	free(tmp); tmp = NULL;
	return 0;

}


int combine_x_v0d_uh(double* x, fdtdMesh* sys, complex<double>* xr) {
	/* Compute xr = V0dn * y0 + uh
	x = [y0; uh]
	xr is the destination with both outside and inside edges, vector size is N_edge */
	myint ind;

	for (ind = 0; ind < sys->v0d1num; ++ind) {
		complex<double> v(0, sys->v0d1valo[ind] * x[sys->v0d1ColIdo[ind]] / sys->v0dn[sys->v0d1ColIdo[ind]]);
		xr[sys->mapioR[sys->v0d1RowId[ind]]] += v;
	}
	for (ind = 0; ind < sys->outside; ++ind) {
		complex<double> v(0, x[sys->leng_v0d1 + ind]);
		xr[sys->mapioR[ind]] += v;
	}
	return 0;
}

int combine_x_v0d_v0b_uh(sparse_matrix_t& V0ddt, sparse_matrix_t& V0bt, myint N1, myint N2, myint N3, double* res, double* x) {
	/* calculate solution as x = V0dn*res1 + V0bn*res2 + res3
	and extend x from just outside to nedge */
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;

	double* v0dnRes1 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, res, beta, v0dnRes1);    // v0dnRes1 = V0dn*res1

	double* v0bnRes2 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &res[N1], beta, v0bnRes2);    // v0bnRes2 = V0bn*res2

	for (int ind = 0; ind < N3; ++ind) {
		x[ind] = v0dnRes1[ind] + v0bnRes2[ind] + res[N1 + N2 + ind];    // x = V0dn*res1 + V0bn*res2 + res3
	}

	free(v0dnRes1); v0dnRes1 = NULL;
	free(v0bnRes2); v0bnRes2 = NULL;

	return 0;
}

void extractColumnVectors(myint* RowId, myint* ColId, double* val, double* aval, myint num, myint lengs, myint lenge, myint* eRowId, myint* eColId, double* eval, double* eaval, myint& e_num) {
	/* extract the column vectors from RowId, ColId, val and aval
	num : nnz in the vectors
	lengs : start column # of the extracted vectors
	lenge : end column # of the extracted vectors
	eRowId : extracted vector rowId
	eColId : extracted vector colId
	eval : extracted vector value
	eaval : extracted vector averaged value
	e_num : end nnz of the extracted vectors */

	myint ind = 0, nums;
	while (ind < num && ColId[ind] < lengs) {
		ind++;
	}
	nums = ind;
	while (ind < num && ColId[ind] <= lenge) {
		eRowId[ind - nums] = RowId[ind];
		eColId[ind - nums] = ColId[ind] - lengs;
		eval[ind - nums] = val[ind];
		eaval[ind - nums] = aval[ind];
		ind++;
	}
	e_num = ind - nums;
}

void transposeMatrix(sparse_matrix_t& V, myint row, myint col, sparse_matrix_t& Vt) {
	/* transpose a matrix
	V : sparse matrix handler for V
	row : row # of V
	col : column # of V
	Vt : sparse matrix handler for Vt */
	myint* idenRowId = (myint*)malloc((row + 1) * sizeof(myint));
	myint* idenColId = (myint*)malloc(row * sizeof(myint));
	double* idenval = (double*)malloc(row * sizeof(double));
	myint ind;

	for (ind = 0; ind < row; ++ind) {    // generate the identity matrix
		idenRowId[ind] = ind;
		idenColId[ind] = ind;
		idenval[ind] = 1;
	}
	idenRowId[ind] = ind;
	sparse_status_t s0;
	sparse_matrix_t I;
	s0 = mkl_sparse_d_create_csr(&I, SPARSE_INDEX_BASE_ZERO, row, row, &idenRowId[0], &idenRowId[1], idenColId, idenval);

	s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, V, I, &Vt);    // Vt = V^T*I

	free(idenRowId); idenRowId = NULL;
	free(idenColId); idenColId = NULL;
	free(idenval); idenval = NULL;
	mkl_sparse_destroy(I);
}

int CSR2COO(myint* rowStart, myint* rowEnd, myint Rows, myint* RowId) {
	/* transfer csr rowId to coo rowId */
	myint row = 0, ind = 0;

	while (row < Rows) {
		while (ind < rowEnd[row]) {
			RowId[ind] = row;
			ind++;
		}
		row++;
	}
	return 0;
}

void findNNZ(myint* ArowStart, myint* ArowEnd, myint* AcolId, double* Aval, myint& leng_A, myint ARows, myint ACols, myint** NrowId, myint** NcolId, double** Nval) {
	/* Check whether the generated matrix has zero entry or not,
	if it has exclude the zero entry and regenerate the rowID, colID and val
	in each row, put the column # to be in ascending order */
	myint rows = 0, i = 0, nnz = 0;

	while (rows < ARows) {
		while (i < ArowEnd[rows]) {
			if (Aval[i] != 0) {
				nnz++;
			}
			i++;
		}
		rows++;
	}

	(*NcolId) = (myint*)malloc(nnz * sizeof(myint));
	(*NrowId) = (myint*)malloc(nnz * sizeof(myint));
	(*Nval) = (double*)malloc(nnz * sizeof(double));
	leng_A = nnz;
	rows = 0;
	i = 0;
	nnz = 0;
	ofstream out;
	//out.open("A22.txt", ofstream::out | ofstream::trunc);
	while (rows < ARows) {
		vector<pair<myint, double>> v;
		while (i < ArowEnd[rows]) {
			if (Aval[i] != 0) {
				v.push_back(make_pair(AcolId[i], Aval[i]));
			}
			i++;
		}
		sort(v.begin(), v.end());   // sort the column # in ascending order to use pardiso
		for (auto vi : v) {
			(*NrowId)[nnz] = rows;
			(*NcolId)[nnz] = vi.first;
			(*Nval)[nnz] = vi.second;
			//out << (*NrowId)[nnz] + 1 << " " << (*NcolId)[nnz] + 1 << " ";
			//out << setprecision(15) << (*Nval)[nnz] << endl;
			nnz++;
		}
		rows++;
		v.clear();
	}
	//out.close();

}

void solveAooMatrix(double* Jo, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* res, void* pt[64], myint iparm[64]) {
	/* solve (-omega^2*D_epsoo+Soo)\(Jo) by the 3*3 system */
	double alpha = 1;
	double beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* b1 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, Jo, beta, b1);   // b1 = V0dan'*Jo

	double* b2 = (double*)calloc(N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, Jo, beta, b2);   // b2 = V0ban'*Jo

	double* b = (double*)calloc(N1 + N2 + N3, sizeof(double));
	for (int indi = 0; indi < N1; ++indi) {
		b[indi] = b1[indi];
	}
	for (int indi = 0; indi < N2; ++indi) {
		b[N1 + indi] = b2[indi];
	}
	for (int indi = 0; indi < N3; ++indi) {
		b[N1 + N2 + indi] = Jo[indi];
	}
	/* begin testing */
	//status = hypreSolve(LooRowId, LooColId, Looval, leng_Loo, Jo, sys->outside, xsol, 1, 3);
	//	//status = hypreSolve(LoodbRowId, LoodbColId, Loodbval, leng_Loodb, Jo, sys->outside, Jo, 1, 3);   // HYPRE to solve (Loo-omega^2*D_epsoo)
	//	//	status = hypreSolve(SooRowId, SooColId, Sooval, lengoo, Jo, sys->outside, yo, 1, 3);   // HYPRE to solve (Soo-omega^2*D_epsoo)
	//	//status = mkl_gmres_A(Jo, yo, LooRowId, LooColId, Looval, leng_Loo, sys->outside);   // gmres to solve (Loo-omega^2*D_epsoo)
	//	/* end testing */

	gmres_V0dV0bLoo(freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, b, res, pt, iparm);

	free(b); b = NULL;
	free(b1); b1 = NULL;
	free(b2); b2 = NULL;
}

void solveAooMatrix_schur(double* Jo, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, double* res, void* pt[64], myint iparm[64]) {
	/* solve (-omega^2*D_epsoo+Soo) by the 3*3 system and the schur complement with gmres */
	double alpha = 1;
	double beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;

	double* b1 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, Jo, beta, b1);   // b1 = V0dan'*Jo

	double* b2 = (double*)calloc(N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, Jo, beta, b2);   // b2 = V0ban'*Jo

	double* b = (double*)calloc(N1 + N2 + N3, sizeof(double));
	for (int indi = 0; indi < N1; ++indi) {
		b[indi] = b1[indi];
	}
	for (int indi = 0; indi < N2; ++indi) {
		b[N1 + indi] = b2[indi];
	}
	for (int indi = 0; indi < N3; ++indi) {
		b[N1 + N2 + indi] = Jo[indi];
	}
	double* A12ib12 = (double*)calloc(N1 + N2, sizeof(double));

	solveA11A22(b, A12ib12, freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm);   // A12ib12 = [A11,A12;A21,A22]\[b1;b2]

	/* test to use A22 schur complement approximation to solve [V0a'*Doos*V0, V0a'*Doos*V0b; V0ba'*Doos*V0, V0ba'*(Doos+Soos)*V0b] */
	//double* temp = (double*)calloc(N1 + N2, sizeof(double));
	//solveA11A22_schur_A22(b, temp, freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm);

	//double nn = 0, nnr = 0;
	//for (int ind = 0; ind < N1 + N2; ++ind) {
	//	nn += pow(A12ib12[ind] - temp[ind], 2);
	//	nnr += pow(A12ib12[ind], 2);
	//}
	//cout << "Relative error with schur A22 is " << nn / nnr << endl;
	/* end of testing to use A22 */



	double* A31b = (double*)calloc(N3, sizeof(double));
	double* A32b = (double*)calloc(N3, sizeof(double));
	double* SooA32b = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, A12ib12, beta, A31b);   // A31b = V0dn*b1
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &A12ib12[N1], beta, A32b);   // A32b = V0bn*b2
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, A32b, beta, SooA32b);

	for (int ind = 0; ind < N3; ++ind) {
		A31b[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];
		A31b[ind] += SooA32b[ind];
		A31b[ind] = b[N1 + N2 + ind] - A31b[ind];
	}


	gmres_solveA33schur(A31b, &res[N1 + N2], freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, pt, iparm);   // x3

	double* Doox3 = (double*)calloc(N3, sizeof(double));
	double* Soox3 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, &res[N1 + N2], beta, Soox3);   // Soox3 = (Soo+Doo)*x3
	for (int ind = 0; ind < N3; ++ind) {
		Doox3[ind] = res[N1 + N2 + ind] * (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // Doox3 = Doo*x3
	}
	double* A13A23x3 = (double*)calloc(N1 + N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, Doox3, beta, A13A23x3);   // A13A23x3_1 = V0dan'*Doo*x3
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, Soox3, beta, &A13A23x3[N1]);   // A13A23x3_2 = V0ban'*(Doo+Soo)*x3
	for (int ind = 0; ind < N1 + N2; ++ind) {
		A13A23x3[ind] = b[ind] - A13A23x3[ind];   // A13A23x3 = [b1;b2] - [A13;A23]*x3
	}
	solveA11A22(A13A23x3, res, freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm);

	free(b1); b1 = NULL;
	free(b2); b2 = NULL;
	free(b); b = NULL;
	free(A12ib12); A12ib12 = NULL;
	free(A31b); A31b = NULL;
	free(A32b); A32b = NULL;
	free(SooA32b); SooA32b = NULL;
	free(Doox3); Doox3 = NULL;
	free(Soox3); Soox3 = NULL;
	free(A13A23x3); A13A23x3 = NULL;
}

void solveA11A22_schur_A11(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]) {
	/* solve [V0dan'*Doo*V0dn, V0dan'*Doo*V0bn;
	V0ban'*Doo*V0dn, V0ban'*(Doo+Soo)*V0bn] with the schur complement
	(A11-A12*A22^(-1)*A21)*x1 = b1-A12*A22^(-1)*b2
	A22*x2 = b2-A21*x1 */
	double* A22ib2 = (double*)calloc(N2, sizeof(double));
	double nor = 0;
	for (int ind = 0; ind < N2; ++ind) {   // when the rhs is larger than 0, then pardisoSol can work
		nor += bm[N1 + ind] * bm[N1 + ind];
	}
	if (nor > 0)
		pardisoSol(A22RowId1, A22ColId, A22Val, pt, iparm, &bm[N1], A22ib2, N2);   // A22ib2 = A22\b2
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* v0bA22ib2 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, A22ib2, beta, v0bA22ib2);   // v0bA22ib2 = V0bn*A22\b2
	for (int ind = 0; ind < N3; ++ind) {
		v0bA22ib2[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // v0bA22ib2 = Doo*V0bn*A22\b2
	}
	double* A12A22ib2 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, v0bA22ib2, beta, A12A22ib2);    // A12A22ib2 = A12*A22^(-1)*b2
	double* rhs_schur = (double*)calloc(N1, sizeof(double));
	for (int ind = 0; ind < N1; ++ind) {
		rhs_schur[ind] = bm[ind] - A12A22ib2[ind];   //rhs_schur = b1-A12*A22^(-1)*b2
	}

	gmres_solveA11schur(rhs_schur, x, freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm);   // x1 = (A11-A12*A22^(-1)*A21)^(-1)*(b1-A12*A22^(-1)*b2)

	double* v0dx1 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, x, beta, v0dx1);    // v0dx1 = V0dn*x1
	for (int ind = 0; ind < N3; ++ind) {
		v0dx1[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];    // v0dx1 = Doo*V0dn*x1
	}
	double* b2 = (double*)calloc(N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, v0dx1, beta, b2);   // b2 = V0ban'*Doo*V0dn*x1
	for (int ind = 0; ind < N2; ++ind) {
		b2[ind] = bm[N1 + ind] - b2[ind];   // b2 = b2 - V0ban'*Doo*V0dn*x1
	}
	nor = 0;
	for (int ind = 0; ind < N2; ++ind) {
		nor += b2[ind] * b2[ind];
	}
	if (nor > 0)
		pardisoSol(A22RowId1, A22ColId, A22Val, pt, iparm, b2, &x[N1], N2);

	free(A22ib2); A22ib2 = NULL;
	free(v0bA22ib2); v0bA22ib2 = NULL;
	free(rhs_schur); rhs_schur = NULL;
	free(v0dx1); v0dx1 = NULL;
	free(b2); b2 = NULL;
}

void solveA11A22_schur_A22(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]) {
	/* solve [V0dan'*Doo*V0dn, V0dan'*Doo*V0bn;
	V0ban'*Doo*V0dn, V0ban'*(Doo+Soo)*V0bn] with the schur complement
	A22*x2 = b2-A21*A11^(-1)*b1
	A11*x1 = b1-A12*x2 */
	int status;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double* A11ib1 = (double*)calloc(N1, sizeof(double));
	status = hypreSolve(A11RowId, A11ColId, A11val, leng_A11, bm, N1, A11ib1, 0, 3);    //A11ib1 = A11^(-1)b1
	for (int ind = 0; ind < N1; ++ind) {
		A11ib1[ind] /= (-pow(freq * 2 * M_PI, 2));    //A11ib1 = A11^(-1)b1
	}

	double* V0A11ib1 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, A11ib1, beta, V0A11ib1);   // V0A11ib1 = V0*A11^(-1)*b1
	for (int ind = 0; ind < N3; ++ind) {
		V0A11ib1[ind] *= Doosval[ind] * (-pow(freq * 2 * M_PI, 2));   // V0A11ib1 = Doo*V0*A11^(-1)*b1
	}
	double* A21A11ib1 = (double*)calloc(N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, V0A11ib1, beta, A21A11ib1);   // A21A11ib1 = A21*A11^(-1)*b1
	for (int ind = 0; ind < N2; ++ind) {
		A21A11ib1[ind] = bm[N1 + ind] - A21A11ib1[ind];   // A21A11ib1 = b2-A21*A11^(-1)*b1
	}

	double nor = 0;
	for (int ind = 0; ind < N2; ++ind) {   // when the rhs is larger than 0, then pardisoSol can work
		nor += A21A11ib1[ind] * A21A11ib1[ind];
	}
	if (nor > 0)
		pardisoSol(A22RowId1, A22ColId, A22Val, pt, iparm, A21A11ib1, &x[N1], N2);   // x2 = A22^(-1)*(b2-A21*A11^(-1)*b1)

	double* V0bx2 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &x[N1], beta, V0bx2);   // V0bx2 = V0b*x2
	for (int ind = 0; ind < N3; ++ind) {
		V0bx2[ind] *= Doosval[ind] * (-pow(freq * 2 * M_PI, 2));   // V0bx2 = Doo*V0b*x2
	}
	double* A12x2 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, V0bx2, beta, A12x2);   // A12x2 = A12*x2
	for (int ind = 0; ind < N1; ++ind) {
		A12x2[ind] = bm[ind] - A12x2[ind];   // A12x2 = b1-A12*x2
	}
	status = hypreSolve(A11RowId, A11ColId, A11val, leng_A11, A12x2, N1, x, 0, 3);    //x1 = A11^(-1)(b1 - A12*x2)
	for (int ind = 0; ind < N1; ++ind) {
		x[ind] /= (-pow(freq * 2 * M_PI, 2));
	}

	free(A11ib1); A11ib1 = NULL;
	free(V0A11ib1); V0A11ib1 = NULL;
	free(A21A11ib1); A21A11ib1 = NULL;
	free(V0bx2); V0bx2 = NULL;
	free(A12x2); A12x2 = NULL;
}


int gmres_solveA11schur(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]) {
	/* solve (A11-A12*A22^(-1)*A21)*x1 = b1 - A12*A22^(-1)*b2 */

	myint size = 128;
	/*------------------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/*------------------------------------------------------------------------------------*/
	MKL_INT ipar[size];
	ipar[14] = 1000;
	//double b[N];
	//double expected_solution[N];
	//double computed_solution[N];
	//double residual[N];
	myint N = N1;
	double* b = (double*)malloc(N * sizeof(double));
	double* expected_solution = (double*)malloc(N * sizeof(double));
	double* computed_solution = (double*)malloc(N * sizeof(double));
	double* residual = (double*)malloc(N * sizeof(double));
	double dpar[size];
	double* tmp = (double*)calloc(N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1, sizeof(double));
	/*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	ivar = N;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double dvar;
	int status, maxit = 1000;
	ofstream out;
	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector b for future use
	/*---------------------------------------------------------------------------*/
	i = 0;
	for (i = 0; i < N; ++i) {
		b[i] = bm[i];
	}
	/*--------------------------------------------------------------------------
	/* Initialize the initial guess
	/*--------------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		computed_solution[i] = 0;// bm[i];
	}

	/*--------------------------------------------------------------------------
	/* Initialize the solver
	/*--------------------------------------------------------------------------*/
	dfgmres_init(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
	ipar[10] = 0;   // no preconditioner
	ipar[14] = 1000;   // restart number
	ipar[7] = 0;
	dpar[0] = 1.0E-5;
	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
ONE: dfgmres(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	//bmn = 0;
	//for (int ind = 0; ind < ivar; ++ind) {
	//	bmn += computed_solution[ind] * computed_solution[ind];
	//}
	//cout << "Right hand side norm is " << bmn << endl;
	if (RCI_request == 0) {
		goto COMPLETE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------*/
	if (RCI_request == 1) {
		/* (A11-A12*A22^(-1)*A21)*tmp[ipar[21] - 1] and put the result in tmp[ipar[22] - 1] */
		A11schurtimesV(freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);
		goto ONE;
	}

	/* do the user-defined stopping test */
	if (RCI_request == 2) {
		ipar[12] = 1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */

		for (myint ind = 0; ind < N; ++ind) {
			residual[ind] = 0;   // before using sparseMatrixVecMul, the resultant vector should be first initialized
		}
		A11schurtimesV(freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm, b, residual);
		dvar = -1.0E0;
		i = 1;
		daxpy(&ivar, &dvar, bm, &i, residual, &i);
		dvar = cblas_dnrm2(ivar, residual, i) / cblas_dnrm2(ivar, bm, i);    // relative residual
																			 //cout << "            The relative residual is " << dvar << " with iteration number " << itercount << endl;
		if (dvar < dpar[0] || itercount > maxit) goto COMPLETE;
		else goto ONE;
	}

	/* apply the preconditioner on the vector tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1] */
	if (RCI_request == 3) {
		/* Apply the preconditioner [A11,A12;0,A22] */


		goto ONE;
	}

	/* check if the norm of the next generated vector is not zero up to rounding and computational errors. */
	if (RCI_request == 4) {
		if (dpar[6] < 1.0E-12) goto COMPLETE;
		else goto ONE;
		//goto ONE;
	}

	else {
		goto FAILED;
	}

COMPLETE: ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "            The relative residual is " << dvar << " with iteration number " << itercount << endl << endl;
	for (i = 0; i < ivar; ++i) {
		x[i] = computed_solution[i];
	}
	goto SUCCEDED;
FAILED: cout << "The solver has returned the ERROR code " << RCI_request << endl;

SUCCEDED: free(b); b = NULL;
	free(expected_solution); expected_solution = NULL;
	free(computed_solution); computed_solution = NULL;
	free(residual); residual = NULL;
	free(tmp); tmp = NULL;
	return 0;

}

void A11schurtimesV(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64], double* v1, double* v2) {
	/* do (A11-A12*A22^(-1)*A21)*v1=v2 */
	ofstream out;
	double* A11v1 = (double*)calloc(N1, sizeof(double));
	double* temp = (double*)calloc(N1, sizeof(double));
	for (int ind = 0; ind < N1; ++ind) {
		temp[ind] = v1[ind] / v0dn[ind];
	}
	sparseMatrixVecMul(A11RowId, A11ColId, A11val, leng_A11, temp, A11v1);   // A11v1 = A11*v1 (This A11 is not normalized!)
	for (int ind = 0; ind < N1; ++ind) {
		A11v1[ind] /= v0dan[ind] / ((-1) * pow(2 * M_PI * freq, 2));
	}

	double* v0dv1 = (double*)calloc(N3, sizeof(double));
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, v1, beta, v0dv1);   // v0dv1 = V0dn*v1
	for (int ind = 0; ind < N3; ++ind) {
		v0dv1[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // v0dv1 = Doo*V0dn*v1
	}

	double* A21v1 = (double*)calloc(N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, v0dv1, beta, A21v1);   // A21v1 = A21*v1
	double* A22iA21v1 = (double*)calloc(N2, sizeof(double));
	double nor = 0;
	for (int ind = 0; ind < N2; ++ind) {
		nor += A21v1[ind] * A21v1[ind];
	}
	if (nor > 0)
		pardisoSol(A22RowId1, A22ColId, A22Val, pt, iparm, A21v1, A22iA21v1, N2);    // A22iA21v1 = A22^(-1)*A21*v1


	double* v0bv = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, A22iA21v1, beta, v0bv);   // v0bv = V0b*A22^(-1)*A21*v1
	for (int ind = 0; ind < N3; ++ind) {
		v0bv[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // v0bv = Doo*V0b*A22^(-1)*A21*v1
	}
	double* A12A22A21 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, v0bv, beta, A12A22A21);   // A12A22A21 = A12*A22^(-1)*A21*v1
	for (int ind = 0; ind < N1; ++ind) {
		v2[ind] = A11v1[ind] - A12A22A21[ind];    //(A11-A12*A22^(-1)*A21)*v1
	}
	free(A11v1); A11v1 = NULL;
	free(v0dv1); v0dv1 = NULL;
	free(A21v1); A21v1 = NULL;
	free(A22iA21v1); A22iA21v1 = NULL;
	free(v0bv); v0bv = NULL;
	free(A12A22A21); A12A22A21 = NULL;
}

int solveA11A22(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64]) {
	/* solve [A11,A12;A21,A22] by gmres and use [A11,A12;0,A22] as the preconditioner */
	myint size = 128;
	/*------------------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/*------------------------------------------------------------------------------------*/
	MKL_INT ipar[size];
	ipar[14] = 1000;
	//double b[N];
	//double expected_solution[N];
	//double computed_solution[N];
	//double residual[N];
	myint N = N1 + N2;
	double* b = (double*)malloc(N * sizeof(double));
	double* expected_solution = (double*)malloc(N * sizeof(double));
	double* computed_solution = (double*)malloc(N * sizeof(double));
	double* residual = (double*)malloc(N * sizeof(double));
	double dpar[size];
	double* tmp = (double*)calloc(N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1, sizeof(double));
	/*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	ivar = N;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double dvar;
	int status, maxit = 1000;
	ofstream out;
	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector b for future use
	/*---------------------------------------------------------------------------*/
	i = 0;
	for (i = 0; i < N; ++i) {
		b[i] = bm[i];
	}
	/*--------------------------------------------------------------------------
	/* Initialize the initial guess
	/*--------------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		computed_solution[i] = 0;// bm[i];
	}

	/*--------------------------------------------------------------------------
	/* Initialize the solver
	/*--------------------------------------------------------------------------*/
	dfgmres_init(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
	ipar[10] = 1;   // use preconditioner
	ipar[14] = 1000;   // restart number
	ipar[7] = 0;
	dpar[0] = 1.0E-6;
	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
ONE: dfgmres(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	//bmn = 0;
	//for (int ind = 0; ind < ivar; ++ind) {
	//	bmn += computed_solution[ind] * computed_solution[ind];
	//}
	//cout << "Right hand side norm is " << bmn << endl;
	if (RCI_request == 0) {
		goto COMPLETE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------*/
	if (RCI_request == 1) {
		/* [A11,A12;A21,A22]*tmp[ipar[21] - 1] and put the result in tmp[ipar[22] - 1] */
		A11A22timesV(Doosval, freq, V0ddt, V0ddat, N1, V0bt, V0bat, N2, Soo, N3, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);

		goto ONE;
	}

	/* do the user-defined stopping test */
	if (RCI_request == 2) {
		ipar[12] = 1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */

		for (myint ind = 0; ind < N; ++ind) {
			residual[ind] = 0;   // before using sparseMatrixVecMul, the resultant vector should be first initialized
		}
		A11A22timesV(Doosval, freq, V0ddt, V0ddat, N1, V0bt, V0bat, N2, Soo, N3, b, residual);
		dvar = -1.0E0;
		i = 1;
		daxpy(&ivar, &dvar, bm, &i, residual, &i);
		dvar = cblas_dnrm2(ivar, residual, i) / cblas_dnrm2(ivar, bm, i);    // relative residual
		//cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
		if (dvar < dpar[0] || itercount > maxit) goto COMPLETE;
		else goto ONE;
	}

	/* apply the preconditioner on the vector tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1] */
	if (RCI_request == 3) {
		/* Apply the preconditioner [A11,A12;0,A22] */
		A11A22Precond(freq, Doosval, V0ddat, N1, v0dn, v0dan, V0bt, N2, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);

		goto ONE;
	}

	/* check if the norm of the next generated vector is not zero up to rounding and computational errors. */
	if (RCI_request == 4) {
		if (dpar[6] < 1.0E-12) goto COMPLETE;
		else goto ONE;
		//goto ONE;
	}

	else {
		goto FAILED;
	}

COMPLETE: ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "            The relative residual is " << dvar << " with iteration number " << itercount << endl << endl;
	for (i = 0; i < ivar; ++i) {
		x[i] = computed_solution[i];
	}
	goto SUCCEDED;
FAILED: cout << "The solver has returned the ERROR code " << RCI_request << endl;

SUCCEDED: free(b); b = NULL;
	free(expected_solution); expected_solution = NULL;
	free(computed_solution); computed_solution = NULL;
	free(residual); residual = NULL;
	free(tmp); tmp = NULL;
	return 0;

}

void A11A22Precond(double freq, double* Doosval, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, myint N2, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, void* pt[64], myint iparm[64], double* v1, double* v2) {
	/* apply [V0dan'*Doo*V0dn, V0dan'*Doo*V0bn;
	0              , V0ban'*(Soo+Doo)*V0bn] * v2 = v1 */
	ofstream out;
	int status;

	/* calculate x2 */
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	pardisoSol(A22RowId1, A22ColId, A22Val, pt, iparm, &v1[N1], &v2[N1], N2);
	//status = hypreSolve(A22RowId, A22ColId, A22Val, leng_A22, &v1[N1], N2, &v2[N1], 1, 3);
	//mkl_gmres_A(&v1[N1], &v2[N1], A22RowId, A22ColId, A22Val, leng_A22, N2);
	//cout << "Time to solve A22 is " << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;

	/* calculate x1 */
	double* v0bnx2 = (double*)calloc(N3, sizeof(double));
	double* Doosx3 = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &v2[N1], beta, v0bnx2);   // v0bnx2 = V0bn*x2
	for (int ind = 0; ind < N3; ++ind) {
		Doosx3[ind] = (v0bnx2[ind]) * (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // Doosx3 = (-omega^2*D_epsoo)*(V0bn*x2)
	}
	double* v0datx2 = (double*)calloc(N1, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, Doosx3, beta, v0datx2);    // v0datx2 = V0da'*(-omega^2*D_epsoo)*(V0b*x2)
	for (int ind = 0; ind < N1; ++ind) {
		v0datx2[ind] = (v1[ind] - v0datx2[ind]);// *v0dan[ind];   // v0datx2x3 = rhs1-V0da'*(-omega^2*D_epsoo)*(V0b*x2)
	}
	//t = clock();
	//mkl_gmres_A(v0datx2x3, x, A11RowId, A11ColId, A11val, leng_A11, N1);
	status = hypreSolve(A11RowId, A11ColId, A11val, leng_A11, v0datx2, N1, v2, 0, 3);
	
	/* pardiso sanity check */
	//myint* A11RowId1 = (myint*)calloc((N1 + 1), sizeof(myint));
	//status = COO2CSR_malloc(A11RowId, A11ColId, A11val, leng_A11, N1, A11RowId1);
	//status = pardisoSolve(A11RowId1, A11ColId, A11val, v0datx2x3, x, N1);
	//cout << "Time to solve A11 is " << (clock() - t) * 1.0 / CLOCKS_PER_SEC << endl;

	for (int ind = 0; ind < N1; ++ind) {
		//v2[ind] *= v0dn[ind];
		v2[ind] /= ((-1) * pow(freq * 2 * M_PI, 2));   // x1 = A11^(-1)*(rhs1-V0da'*(-omega^2*D_epsoo)*(V0b*x2+x3))
	}

	//free(LooRowId1); LooRowId1 = NULL;
	//free(A11RowId1); A11RowId1 = NULL;
	free(Doosx3); Doosx3 = NULL;
	free(v0bnx2); v0bnx2 = NULL;
	free(v0datx2); v0datx2 = NULL;
}

void A11A22timesV(double* Doosval, double freq, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, double* v, double* res) {
	/* Do [v0dan'*Doos*v0dn, v0dan'*Doos*v0bn;
	v0ban'*Doos*v0dn, v0ban'*(Doos+Soos)*v0bn] multiply a vector */
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	ofstream out;


	/* first line */
	double* v0dnV1 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, v, beta, v0dnV1);   // v0dnV1 = V0dn*v1
	double* v0bnV2 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &v[N1], beta, v0bnV2);   // v0bnV2 = V0bn*v2

	double* temp = (double*)malloc(N3 * sizeof(double));
	for (int ind = 0; ind < N3; ++ind) {
		temp[ind] = v0dnV1[ind] + v0bnV2[ind];   // temp = V0dn*v1 + V0bn*v2
		temp[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];   // temp = (-omega^2*D_epsoo)*(V0dn*v1 + V0bn*v2)
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, temp, beta, res);   // res1 = V0dan'*(-omega^2*D_epsoo)*(V0dn*v1 + V0bn*v2)

																								  /* second line */
	double* DoosV0dnV1 = (double*)malloc(N3 * sizeof(double));
	for (int ind = 0; ind < N3; ++ind) {
		DoosV0dnV1[ind] = (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind] * v0dnV1[ind];    // DoosV0dnV1 = -omega^2*D_epsoo*V0dn*v1
	}
	double* SoosV0bnV2 = (double*)malloc(N3 * sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, v0bnV2, beta, SoosV0bnV2);   // SoosV0bnV2 = (-omega^2*D_epsoo+Soos)*V0bn*V2
	for (int ind = 0; ind < N3; ++ind) {
		temp[ind] = DoosV0dnV1[ind] + SoosV0bnV2[ind];    // temp = (-omega^2*D_epsoo)*V0dn*V1+(-omega^2*D_epsoo+Soos)*V0bn*V2
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, temp, beta, &res[N1]);    // res2 = V0ban'*((-omega^2*D_epsoo)*V0dn*V1+(-omega^2*D_epsoo+Soos)*(V0bn*V2+V3))


	free(v0dnV1); v0dnV1 = NULL;
	free(v0bnV2); v0bnV2 = NULL;
	free(temp); temp = NULL;
	free(DoosV0dnV1); DoosV0dnV1 = NULL;
	free(SoosV0bnV2); SoosV0bnV2 = NULL;

}

int gmres_solveA33schur(double* bm, double* x, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, void* pt[64], myint iparm[64]) {
	/* gmres(A33-[A31,A32]*[A11,A12;A21,A22]^(-1)*[A13,A23], bm)*/
	myint size = 128;
	/*------------------------------------------------------------------------------------
	/* Allocate storage for the ?par parameters and the solution/rhs/residual vectors
	/*------------------------------------------------------------------------------------*/
	MKL_INT ipar[size];
	ipar[14] = 500;
	//double b[N];
	//double expected_solution[N];
	//double computed_solution[N];
	//double residual[N];
	myint N = N3;
	double* b = (double*)malloc(N * sizeof(double));
	double* expected_solution = (double*)malloc(N * sizeof(double));
	double* computed_solution = (double*)malloc(N * sizeof(double));
	double* residual = (double*)malloc(N * sizeof(double));
	double dpar[size];
	double* tmp = (double*)calloc(N*(2 * ipar[14] + 1) + (ipar[14] * (ipar[14] + 9)) / 2 + 1, sizeof(double));
	/*---------------------------------------------------------------------------
	/* Some additional variables to use with the RCI (P)FGMRES solver
	/*---------------------------------------------------------------------------*/
	MKL_INT itercount;
	MKL_INT RCI_request, i, ivar;
	ivar = N;
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	double dvar;
	int status, maxit = 500;
	ofstream out;
	/*---------------------------------------------------------------------------
	/* Save the right-hand side in vector b for future use
	/*---------------------------------------------------------------------------*/
	i = 0;
	for (i = 0; i < N; ++i) {
		b[i] = bm[i];
	}
	/*--------------------------------------------------------------------------
	/* Initialize the initial guess
	/*--------------------------------------------------------------------------*/
	for (i = 0; i < N; ++i) {
		computed_solution[i] = 0;// bm[i];
	}

	/*--------------------------------------------------------------------------
	/* Initialize the solver
	/*--------------------------------------------------------------------------*/
	dfgmres_init(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
	ipar[10] = 0;   // no preconditioner
	ipar[14] = 500;   // restart number
	ipar[7] = 0;
	dpar[0] = 1.0E-3;
	/*---------------------------------------------------------------------------
	/* Check the correctness and consistency of the newly set parameters
	/*---------------------------------------------------------------------------*/
	dfgmres_check(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	if (RCI_request != 0) goto FAILED;
ONE: dfgmres(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp);
	//bmn = 0;
	//for (int ind = 0; ind < ivar; ++ind) {
	//	bmn += computed_solution[ind] * computed_solution[ind];
	//}
	//cout << "Right hand side norm is " << bmn << endl;
	if (RCI_request == 0) {
		goto COMPLETE;
	}

	/*---------------------------------------------------------------------------
	/* If RCI_request=1, then compute the vector A*tmp[ipar[21]-1]
	/* and put the result in vector tmp[ipar[22]-1]
	/*---------------------------------------------------------------------------*/
	if (RCI_request == 1) {
		/* (A33-[A31,A32]*[A11,A12;A21,A22]^(-1)*[A13;A23])*tmp[ipar[21] - 1] and put the result in tmp[ipar[22] - 1] */
		A33schurTimesV(freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, pt, iparm, &tmp[ipar[21] - 1], &tmp[ipar[22] - 1]);

		goto ONE;
	}

	/* do the user-defined stopping test */
	if (RCI_request == 2) {
		ipar[12] = 1;
		/* Get the current FGMRES solution in the vector b[N] */
		dfgmres_get(&ivar, computed_solution, b, &RCI_request, ipar, dpar, tmp, &itercount);
		/* Compute the current true residual via MKL (Sparse) BLAS routines */

		for (myint ind = 0; ind < N; ++ind) {
			residual[ind] = 0;   // before using sparseMatrixVecMul, the resultant vector should be first initialized
		}
		A33schurTimesV(freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, pt, iparm, b, residual);
		dvar = -1.0E0;
		i = 1;
		daxpy(&ivar, &dvar, bm, &i, residual, &i);
		dvar = cblas_dnrm2(ivar, residual, i) / cblas_dnrm2(ivar, bm, i);    // relative residual
		cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
		if (dvar < dpar[0] || itercount > maxit) goto COMPLETE;
		else goto ONE;
	}

	/* apply the preconditioner on the vector tmp[ipar[21]-1] and put the result in vector tmp[ipar[22]-1] */
	if (RCI_request == 3) {

		goto ONE;
	}

	/* check if the norm of the next generated vector is not zero up to rounding and computational errors. */
	if (RCI_request == 4) {
		if (dpar[6] < 1.0E-12) goto COMPLETE;
		else goto ONE;
		//goto ONE;
	}

	else {
		goto FAILED;
	}

COMPLETE: ipar[12] = 0;
	dfgmres_get(&ivar, computed_solution, bm, &RCI_request, ipar, dpar, tmp, &itercount);
	cout << "The relative residual is " << dvar << " with iteration number " << itercount << endl;
	for (i = 0; i < ivar; ++i) {
		x[i] = computed_solution[i];
	}
	goto SUCCEDED;
FAILED: cout << "The solver has returned the ERROR code " << RCI_request << endl;

SUCCEDED: free(b); b = NULL;
	free(expected_solution); expected_solution = NULL;
	free(computed_solution); computed_solution = NULL;
	free(residual); residual = NULL;
	free(tmp); tmp = NULL;
	return 0;
}

void A33schurTimesV(double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint N1, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint N2, sparse_matrix_t& Soo, myint N3, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, void* pt[64], myint iparm[64], double* v1, double* v2) {
	/* (A33-[A31,A32]*[A11,A12;A21,A22]^(-1)*[A13;A23])*v1=v2 */

	double* A33v1 = (double*)calloc(N3, sizeof(double));
	sparseMatrixVecMul(LooRowId, LooColId, Looval, leng_Loo, v1, A33v1);   // A33v1 = A33*v1

	double* Soov1 = (double*)calloc(N3, sizeof(double));
	double alpha = 1, beta = 0;
	struct matrix_descr descr;
	descr.type = SPARSE_MATRIX_TYPE_GENERAL;
	sparse_status_t s;
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, v1, beta, Soov1);   // Soov1 = (Soo+Doo)*v1
	double* Doov1 = (double*)calloc(N3, sizeof(double));
	for (int ind = 0; ind < N3; ++ind) {
		Doov1[ind] = Doosval[ind] * (-1) * pow(freq * 2 * M_PI, 2) * v1[ind];   // Doov1 = Doo*v1
	}
	double* A13A23v1 = (double*)calloc(N1 + N2, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0ddat, descr, Doov1, beta, A13A23v1);   // A13A23v1[1:N1] = V0dan'*Doo*v1
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0bat, descr, Soov1, beta, &A13A23v1[N1]);   // A13A23v1[N1+1:N1+N2] = V0ban'*(Doo+Soo)*v1

	double* temp = (double*)calloc(N1 + N2, sizeof(double));
	solveA11A22(A13A23v1, temp, freq, Doosval, V0ddt, V0ddat, N1, v0dn, v0dan, V0bt, V0bat, N2, Soo, N3, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, pt, iparm);  // temp=[A11,A12;A21,A22]^(-1)*[A13;A23]*v1

	double* A31b = (double*)calloc(N3, sizeof(double));
	double* A32b = (double*)calloc(N3, sizeof(double));
	double* SooA32b = (double*)calloc(N3, sizeof(double));
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ddt, descr, temp, beta, A31b);   // A31b = V0dn*b1
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0bt, descr, &temp[N1], beta, A32b);   // A32b = V0bn*b2
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, Soo, descr, A32b, beta, SooA32b);   // SooA32b = (Doo+Soo)*V0bn*b2
	for (int ind = 0; ind < N3; ++ind) {
		A31b[ind] *= (-1) * pow(freq * 2 * M_PI, 2) * Doosval[ind];
		A31b[ind] += SooA32b[ind];
		v2[ind] = A33v1[ind] - A31b[ind];
	}

	free(A33v1); A33v1 = NULL;
	free(Soov1); Soov1 = NULL;
	free(Doov1); Doov1 = NULL;
	free(A13A23v1); A13A23v1 = NULL;
	free(temp); temp = NULL;
	free(A31b); A31b = NULL;
	free(A32b); A32b = NULL;
	free(SooA32b); SooA32b = NULL;
}

void solveMatrix(double* Jo, double freq, double* Doosval, sparse_matrix_t& V0ddt, sparse_matrix_t& V0ddat, myint leng_v0d, double* v0dn, double* v0dan, sparse_matrix_t& V0bt, sparse_matrix_t& V0bat, myint leng_v0b, sparse_matrix_t& Soo, myint outside, myint* A11RowId, myint* A11ColId, double* A11val, myint leng_A11, myint* A22RowId, myint* A22RowId1, myint* A22ColId, double* A22Val, myint leng_A22, myint* LooRowId, myint* LooColId, double* Looval, myint leng_Loo, myint* SoiRowId, myint* SoiColId, double* Soival, myint lengoi, myint* SioRowId, myint* SioColId, double* Sioval, myint lengio, myint* mapio, complex<double>* xt, myint nedge, void* pt[64], myint iparm[64]) {
	/* solve (-1i*omega*D_sigii)*xi=-Sio*(Aoo\(1i*Jo))
	and Aoo*xo=1i*Jo-Soi*xi
	Jo : -omega*J
	xt : final total solution */
	ofstream out;

	// AooiJo = Aoo^(-1)*(Jo)
	int status;
	double* res = (double*)calloc(leng_v0b + leng_v0d + outside, sizeof(double));
	double* AooiJo = (double*)calloc(outside, sizeof(double));
	//solveAooMatrix(Jo, freq, Doosval, V0ddt, V0ddat, leng_v0d, v0dn, v0dan, V0bt, V0bat, leng_v0b, Soo, outside, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, res, pt, iparm);   // solve the 3*3 system with the upper triangular matrix as the preconditioner
	solveAooMatrix_schur(Jo, freq, Doosval, V0ddt, V0ddat, leng_v0d, v0dn, v0dan, V0bt, V0bat, leng_v0b, Soo, outside, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, res, pt, iparm);   // solve by the schur complement with gmres
	status = combine_x_v0d_v0b_uh(V0ddt, V0bt, leng_v0d, leng_v0b, outside, res, AooiJo);
	free(res);

	// xi = (1i*omega*D_sigii)\(-Sio*(Aoo\(-1i*omega*Jo)))
	double* xi = (double*)calloc(nedge - outside, sizeof(double));    // initialized the resultant vector to be 0
	sparseMatrixVecMul(SioRowId, SioColId, Sioval, lengio, AooiJo, xi);
	for (int ind = 0; ind < nedge - outside; ++ind) {
		xi[ind] /= (-1) * freq * 2 * M_PI * SIGMA;
	}
	//out.open("xi.txt", ofstream::out | ofstream::trunc);
	//for (int ind = 0; ind < nedge - outside; ++ind) {
	//	out << setprecision(15) << xi[ind] << endl;
	//}
	//out.close();

	// AooiSoixi = Aoo\(-Soi*xi)
	double* Soixi = (double*)calloc(outside, sizeof(double));
	sparseMatrixVecMul(SoiRowId, SoiColId, Soival, lengoi, xi, Soixi);
	for (int ind = 0; ind < outside; ++ind) {
		Soixi[ind] *= -1;
	}
	res = (double*)calloc(leng_v0b + leng_v0d + outside, sizeof(double));
	//solveAooMatrix(Soixi, freq, Doosval, V0ddt, V0ddat, leng_v0d, v0dn, v0dan, V0bt, V0bat, leng_v0b, Soo, outside, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, res, pt, iparm);    // solve the 3*3 system with the upper triangular matrix as the preconditioner
	solveAooMatrix_schur(Soixi, freq, Doosval, V0ddt, V0ddat, leng_v0d, v0dn, v0dan, V0bt, V0bat, leng_v0b, Soo, outside, A11RowId, A11ColId, A11val, leng_A11, A22RowId, A22RowId1, A22ColId, A22Val, leng_A22, LooRowId, LooColId, Looval, leng_Loo, res, pt, iparm);   // solve by the schur complement with gmres

	double* AooiSoixi = (double*)calloc(outside, sizeof(double));
	status = combine_x_v0d_v0b_uh(V0ddt, V0bt, leng_v0d, leng_v0b, outside, res, AooiSoixi);

	// xo = Aoo\(1i*Jo)-Aoo\(Soi*xi)
	complex<double>* xo = (complex<double>*)calloc(outside, sizeof(complex<double>));
	for (int ind = 0; ind < outside; ++ind) {
		xo[ind] = AooiSoixi[ind] + 1i * AooiJo[ind];
	}
	//out.open("xo.txt", ofstream::out | ofstream::trunc);
	//for (int ind = 0; ind < outside; ++ind) {
	//	out << setprecision(15) << xo[ind] << endl;
	//}
	//out.close();

	for (int ind = 0; ind < nedge; ++ind) {
		if (mapio[ind] < outside) {
			xt[ind] = xo[mapio[ind]];
		}
		else {
			xt[ind] = xi[mapio[ind] - outside] + 1i * 0;
		}
	}

	free(res); res = NULL;
	free(AooiJo); AooiJo = NULL;
	free(xi); xi = NULL;
	free(Soixi); Soixi = NULL;
	free(AooiSoixi); AooiSoixi = NULL;
	free(xo); xo = NULL;
}

void DTimesV(fdtdMesh* sys, myint nedge, double freq, complex<double>* u1, complex<double>* u2) {   
	/* u2 = (-omega^2*D_eps+1i*omega*D_sig)*u1 */
	for (int ind = 0; ind < nedge; ++ind) {
		complex<double> epsig;
		if (sys->markEdge[sys->mapEdgeR[ind]])
			epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * freq * 2 * M_PI * SIGMA;
		else
			epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * 0;
		u2[ind] = epsig * u1[ind];   // u2 = (-omega^2*D_eps+1i*omega*D_sig)*
	}
}

void solve_L_schur(complex<double>* u2, complex<double>* uh, myint nedge, double freq, double* Doosval, sparse_matrix_t& V0dt, sparse_matrix_t& V0dat, myint leng_v0d, sparse_matrix_t& V0ct, sparse_matrix_t& V0cat, myint leng_v0c, myint* LrowId, myint* LcolId, double* Lval, myint leng_L) {
	/* solve {(D+L)-D*V0*(V0a'*D*V0)^(-1)*V0a'*D}*uh = u2 = -1i*omega*J-D*V0*(V0a'*D*V0)^(-1)*(V0a'*(-1i*omega*J)) */
	int status, maxIterNum = 500;
	double epsilon = 1e-3;

	status = GMRES(uh, u2, nedge, epsilon, maxIterNum, freq, Doosval, V0dt, V0dat, leng_v0d, V0ct, V0cat, leng_v0c, LrowId, LcolId, Lval, leng_L, )
}


int GMRES(complex<double> *x, complex<double> *b, myint n, double epsilon, int maxIterNum, double freq, double* Doosval, sparse_matrix_t& V0dt, sparse_matrix_t& V0dat, myint leng_v0d, sparse_matrix_t& V0ct, sparse_matrix_t7 V0cat, myint leng_v0c, myint* LrowId, myint* LcolId, double* Lval, myint leng_L) {
	// This function implements GMRES 
	// refer to Iterative Methods Linear System Second Edition Yousef Saad p171 for theory
	// Also refer to https://rtraba.files.wordpress.com/2015/05/cpp_numerical.pdf for some implementing issue
	// However, I believe implementation above is not as good as mine regarding the complexity, since it changes matrix size every iteration
	// And https://en.wikipedia.org/wiki/Generalized_minimal_residual_method, which make easy to refresh the memory.
	/* x is the solution, b is the rhs */

	int i, j, ii, k;

	complex<double> * r0 = (complex<double>*) calloc(n, sizeof(complex<double>));   // r0 = b - A * x0
	complex<double> * v = (complex<double>*) calloc(n, sizeof(complex<double>));
	complex<double> * w = (complex<double> *) calloc(n, sizeof(complex<double>));

	int unisize = 2;
	int Vsize = unisize + 1;
	int Hsize = unisize;
	int Qsize = unisize + 1;
	int gsize = unisize + 1;
	complex<double> * V = (complex<double> *) calloc(n * Vsize, sizeof(complex<double>));
	complex<double> * H = (complex<double> *) calloc((Hsize + 1)*Hsize, sizeof(complex<double>));
	complex<double> * Q = (complex<double> *) calloc(Qsize*Qsize, sizeof(complex<double>));
	complex<double> * tempQ = (complex<double> *) calloc(Qsize*Qsize, sizeof(complex<double>));
	complex<double> * Omega = (complex<double> *) calloc(Qsize*Qsize, sizeof(complex<double>));
	complex<double> * bg = (complex<double> *) calloc(gsize, sizeof(complex<double>));
	complex<double> * bgtemp = (complex<double> *) calloc(gsize, sizeof(complex<double>));
	complex<double> * Htemp = (complex<double> *) calloc(n, sizeof(complex<double>));
	complex<double> * Htemptemp = (complex<double> *) calloc(n, sizeof(complex<double>));
	for (i = 0; i < Qsize; i++) {
		Q[i + i * Qsize] = 1.0;
		Omega[i + i * Qsize] = 1.0;
	}
	// For every iteration, V is n*Vsize, H is (Hsize+1)*Hsize, Q is Qsize*Qsize, where Qsize=Hsize+1
	// Omega is Qsize*Qsize, bg is Qsize
	complex<double> Ht_k_k, Ht_kp1_k, sk, ck;

	// for initial guess, we choose x0 = 0, thus, r0 = b

	for (i = 0; i<n; i++) {
		r0[i] = b[i]; //r0 = b; line 1
	}


	double beta = V2norm(r0, n);    // beta = norm(ro)
	//printf("beta = %.20f \n", beta);
	bg[0] = { beta, 0 };
	for (i = 1; i < gsize; i++) {
		bg[i] = { 0.0, 0.0 };
	}
	double residual = 1.0;
	for (i = 0; i < n; i++) {
		v[i] = r0[i] / beta;
	}
	copyMatrix(v, V, n, 1);

	// for (i = 0; i < n; i++){
	// 	printf("v[%d] = %.20f%+.20fi \n", i, creal(v[i]), cimag(v[i]));
	// }

	k = 0;
	while (residual > epsilon && k < maxIterNum) {
		// 0 resize V, H, J, Jtotal is k > current size
		if (k >= unisize) {
			unisize = 2 * k;
			// 0.1 resize V
			Vsize = unisize + 1;
			V = (double _Complex *) realloc(V, n*Vsize * sizeof(double _Complex));
			// 0.2 resize H
			Hsize = unisize;
			resize_matrix_upper_left_zero_based(&H, k + 1, k, Hsize + 1, Hsize);
			// 0.3 resize Q, Omega
			Qsize = Hsize + 1;
			resize_matrix_upper_left_zero_based(&Q, k + 1, k + 1, Qsize, Qsize);
			for (ii = k + 1; ii < Qsize; ii++) {
				Q[ii + ii * Qsize] = 1.0;
			}
			tempQ = (double _Complex *) realloc(tempQ, Qsize*Qsize * sizeof(double _Complex));
			resize_matrix_and_fill_with_identity(&Omega, Qsize);
			// 0.4 resize bar{g}
			gsize = Qsize;
			bg = (double _Complex *) realloc(bg, Qsize * sizeof(double _Complex));
			for (ii = k + 1; ii < Qsize; ii++) {
				bg[ii] = 0.0;
			}
			bgtemp = (double _Complex *) realloc(bgtemp, Qsize * sizeof(double _Complex));
		}

		// 1, Arnordi process
		H2MV(psm, pcb, V + k*n, w, n, P_inv, A0_inv, FMMorNew);
		L_schur_mv()
		// since we reuse H[0:k-1,k] and V[:,1:k]
		// only compute H[:,k], and V[:,k+1]

		for (i = 0; i <= k; i++) {
			copyMatrix(V + i*n, v, n, 1);
			H[i + k * (Hsize + 1)] = hermitian_inner_product(w, v, n);
			for (ii = 0; ii < n; ii++) {
				w[ii] = w[ii] - H[i + k * (Hsize + 1)] * v[ii];
			}
		}

		H[k + 1 + k * (Hsize + 1)] = V2norm(w, n);

		for (ii = 0; ii < n; ii++) {
			V[ii + (k + 1)*n] = w[ii] / H[k + 1 + k * (Hsize + 1)];
		}

		// 2, Least Square process
		// 2.1 multiply the last column of H by all the previous Q, denote as Htemp		
		copyMatrix(H + k*(Hsize + 1), Htemptemp, k + 2, 1);
		ZGMV(Q, Htemptemp, Htemp, Qsize, Qsize);

		// 2.2 compute the new plane rotation
		Ht_k_k = Htemp[k];
		Ht_kp1_k = Htemp[k + 1];
		sk = Ht_kp1_k / csqrt(cpow(Ht_k_k, 2) + cpow(Ht_kp1_k, 2));
		ck = Ht_k_k / csqrt(cpow(Ht_k_k, 2) + cpow(Ht_kp1_k, 2));
		Omega[k + k * Qsize] = ck;
		Omega[k + (k + 1) * Qsize] = sk;
		Omega[k + 1 + k * Qsize] = -1 * sk;
		Omega[k + 1 + (k + 1) * Qsize] = ck;

		// 2.3 multiply this new plane rotation to previous computed plane rotation, denoted as Q. 
		// Thus (QH=Rm) is upper triangle matrix.
		copyMatrix(bg, bgtemp, gsize, 1);
		ZGMV(Omega, bgtemp, bg, Qsize, Qsize);

		copyMatrix(Q, tempQ, Qsize, Qsize);
		MatrixMatrix(Omega, tempQ, Q, Qsize, Qsize, Qsize);

		// 2.4 compute residual
		residual = cabs(bg[k + 1]) / beta;
		printf("k = %d, residual = %.20f \n", k, residual);

		// 2.5 make Omega (plane rotation) identity again for next iteration.
		Omega[k + k * Qsize] = 1.0;
		Omega[k + (k + 1) * Qsize] = 0.0;
		Omega[k + 1 + k * Qsize] = 0.0;
		Omega[k + 1 + (k + 1) * Qsize] = 1.0;

		k++;
	}

	// 3, After converge, compute result: xm = x0 + V_my_m, where y_m = min(\bar{g} - \bar{Rm}y)
	// 3.0 get Rm:
	double _Complex * Rm = (double _Complex *) calloc(Qsize * Hsize, sizeof(double _Complex));
	MatrixMatrix(Q, H, Rm, Qsize, Qsize, Hsize);
	// 3.1 compute y_m = min(\bar{g} - \bar{Rm}y). \bar{Rm} is upper triangle version of Hm
	double _Complex * y = (double _Complex *) calloc(k, sizeof(double _Complex));
	// This is equivalent to solving an upper triangle linear system
	// Right now \bar{Rm} is of size (k+1) * k.
	// also consider use MKL's pztrtrs function
	for (i = k - 1; i >= 0; i--) {
		y[i] = bg[i];
		for (j = k - 1; j > i; j--) {
			y[i] = y[i] - Rm[i + j * Qsize] * y[j];
		}
		y[i] = y[i] / Rm[i + i * Qsize];
	}
	// 3.2 compute xm = x0 + V_my_m
	ZGMV(V, y, x, n, k);

	// test ||bg - Rmy||
	double _Complex * Ry = (double _Complex *) calloc(Qsize, sizeof(double _Complex));
	ZGMV(Rm, y, Ry, Qsize, k);
	printf("||bg - Rmy|| is %.20f \n", relaErrorFVV(bg, Ry, Hsize + 1));

	// test ||beta e1 - \bar{Hm}y||
	double _Complex * Hy = (double _Complex *) calloc(Hsize + 1, sizeof(double _Complex));
	ZGMV(H, y, Hy, Hsize + 1, k);
	double _Complex * be = (double _Complex *) calloc(Hsize + 1, sizeof(double _Complex));
	be[0] = beta;
	printf("||beta e1 - bar{Hm}y|| is %.20f \n", relaErrorFVV(be, Hy, Hsize + 1));

	// test ||b - Zx||
	double _Complex * Zx = (double _Complex *) calloc(n, sizeof(double _Complex));
	H2MV(psm, pcb, x, Zx, n, P_inv, A0_inv, FMMorNew);
	printf("||b - Zx|| is %.20f \n", relaErrorFVV(b, Zx, n));

	// 4, free memory:
	free(r0); r0 = NULL;
	free(v); v = NULL;
	free(w); w = NULL;
	free(V); V = NULL;
	free(H); H = NULL;
	free(Q); Q = NULL;
	free(tempQ); tempQ = NULL;
	free(Omega); Omega = NULL;
	free(bg); bg = NULL;
	free(bgtemp); bgtemp = NULL;
	free(Htemp); Htemp = NULL;
	free(Htemptemp); Htemptemp = NULL;
	free(Rm); Rm = NULL;
	free(y); y = NULL;

	// return number of iterations
	return k;
}

double V2norm(complex<double>* r0, myint n) {
	/* calculate r0 vector's norm */
	double res = 0;

	for (int i = 0; i < n; ++i) {
		res += r0[i].real() * r0[i].real() + r0[i].imag() * r0[i].imag();
	}
	res = sqrt(res);
	return res;
}