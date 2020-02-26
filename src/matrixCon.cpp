//#include "stdafx.h"
#include <ctime>
#include "fdtd.hpp"
#include "hypreSolver.h"
//#define TIME (1)   // TIME = 1 run time domain
#define FREQ (1)   // define FREQ run frequency domain
//#define SOLVECHIP (1)   // define SOLVECHIP in the frequency domain
//#define GENERATE_V0_SOLUTION
#define SKIP_VH
void addOmegaEpsilon(fdtdMesh* sys, myint* RowId, myint* ColId, double* val, myint leng, myint N, double freq, double* val1);
void solveV0(fdtdMesh* sys, double freq, double* rhs, complex<double>* u0, sparse_matrix_t& v0d, sparse_matrix_t& v0da, sparse_matrix_t& v0c, sparse_matrix_t& v0ca, char ri);   // solve (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0) with real rhs
void V0Multiply(fdtdMesh* sys, sparse_matrix_t& V0dt, sparse_matrix_t& V0ct, myint lengv0d, myint lengv0c, myint row, complex<double>* u1, complex<double>* u2);
void V0aMultiply(fdtdMesh* sys, sparse_matrix_t& V0dat, sparse_matrix_t& V0cat, myint lengv0d, myint lengv0c, myint nedge, complex<double>* u2, complex<double>* u1);
void generateP(fdtdMesh* sys, double freq, myint* LRowId, myint* LColId, double* Lval, myint leng_L, myint N, myint* PRowId, myint* PColId, double* Pval, myint& leng_P);
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


	cout << "Length of V0d1 is " << sys->leng_v0d1 << ", and number of non-zeros in V0d1 is " << sys->v0d1num << endl;
	cout << "Length of V0d1a is " << sys->leng_v0d1a << ", and number of non-zeros in V0d1a is " << sys->v0d1anum << endl;
	cout << "V0d is generated!" << endl;

	/* Generate Ad = V0da1'*D_eps*V0d */
	sys->generateAd(sys->mapd, sys->v0d1num, sys->v0d1anum, sys->leng_v0d1);

	/* Begin to map the inside and outside edges */
	sys->mapEdgeInsideOutside();   // generate mapio and mapioR
	/* End of mapping the inside and outside edges */
	cout << "The number of outside conductor edge is " << sys->outside << endl;
	cout << "The number of inside conductor edge is " << sys->N_edge - sys->bden - sys->outside << endl;

	sys->v0d1valo = (double*)malloc(sys->v0d1num * sizeof(double));
	for (indi = 0; indi < sys->v0d1num; indi++) {
		sys->v0d1valo[indi] = sys->v0d1val[indi];
	}
	for (indi = 0; indi < sys->v0d1num; indi++) {    // compute sqrt(D_eps)*V0d1
		sys->v0d1val[indi] *= sqrt(sys->stackEpsn[(sys->v0d1RowId[indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
	}
	sys->v0d1avalo = (double*)malloc(sys->v0d1anum * sizeof(double));
	for (indi = 0; indi < sys->v0d1anum; indi++) {
		sys->v0d1avalo[indi] = sys->v0d1aval[indi];
	}
	for (indi = 0; indi < sys->v0d1anum; indi++) {
		sys->v0d1aval[indi] *= sqrt(sys->stackEpsn[(sys->v0d1RowId[indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
	}


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
#ifdef FREQ
	/* Begin to set v0d to be the outside v0d */
	for (indi = 0; indi < sys->v0d1num; indi++) {
		//if (sys->markEdge[sys->v0d1RowId[indi]])
		//	cout << "Something wrong with V0d!" << endl;
		sys->v0d1RowId[indi] = sys->mapio[sys->mapEdge[sys->v0d1RowId[indi]]];
	}
	/* End of setting v0d to be the outside v0d */
	V0d_row = sys->outside;
#endif


	//cout << "Number of NNZ in V0d1 is " << sys->v0d1num << endl;

	sparse_status_t s;

	/* V0d^T's csr form handle for MKL */


	sparse_matrix_t V0dt;
	s = mkl_sparse_d_create_csr(&V0dt, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, V0d_row, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1valo);

	/* V0da^T's csr form handle for MKL */
	sparse_matrix_t V0dat;
	s = mkl_sparse_d_create_csr(&V0dat, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, V0d_row, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1avalo);
	sparse_matrix_t V0da1t;   // Approximation of V0da
	s = mkl_sparse_d_create_csr(&V0da1t, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, V0d_row, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1aval1);


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

#ifdef PRINT_VERBOSE_TIMING
	//cout << "Time to generate V0c is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif
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

	free(sys->markNode); sys->markNode = NULL;

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

	/* Print out V0d, V0da, V0c, V0ca */
	ofstream outfile;
	//outfile.open("V0d.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0d1num; ++indi) {
	//	outfile << sys->v0d1RowId[indi] + 1 << " " << sys->v0d1ColIdo[indi] + 1 << " " << sys->v0d1valo[indi] / sys->v0dn[sys->v0d1ColIdo[indi]] << endl;
	//}
	//outfile.close();

	//outfile.open("V0da.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0d1num; ++indi) {
	//	outfile << sys->v0d1RowId[indi] + 1 << " " << sys->v0d1ColIdo[indi] + 1 << " " << sys->v0d1avalo[indi] / sys->v0dan[sys->v0d1ColIdo[indi]] << endl;
	//}
	//outfile.close();

	//outfile.open("V0c.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0cnum; ++indi) {
	//	outfile << sys->v0cRowId[indi] + 1 << " " << sys->v0cColIdo[indi] + 1 << " " << sys->v0cvalo[indi] / sys->v0cn[sys->v0cColIdo[indi]] << endl;
	//}
	//outfile.close();

	//outfile.open("V0ca.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0cnum; ++indi) {
	//	outfile << sys->v0cRowId[indi] + 1 << " " << sys->v0cColIdo[indi] + 1 << " " << sys->v0cavalo[indi] / sys->v0can[sys->v0cColIdo[indi]] << endl;
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
	ofstream out;


	/* Begin to generate Soo, Soi, Sio */
#ifdef FREQ
	myint lengoo = 0, lengoi = 0, lengio = 0;   // nnz for each submatrix
	matrixInsideOutside_count(sys, sys->SRowId, sys->SColId, sys->leng_S, lengoo, lengoi, lengio);   // count nnz in each submatrix
	myint* SooRowId = (myint*)malloc(lengoo * sizeof(myint));
	myint* SooColId = (myint*)malloc(lengoo * sizeof(myint));
	double* Sooval = (double*)malloc(lengoo * sizeof(double));
	myint* SoiRowId = (myint*)malloc(lengoi * sizeof(myint));
	myint* SoiColId = (myint*)malloc(lengoi * sizeof(myint));
	double* Soival = (double*)malloc(lengoi * sizeof(double));
	myint* SioRowId = (myint*)malloc(lengio * sizeof(myint));
	myint* SioColId = (myint*)malloc(lengio * sizeof(myint));
	double* Sioval = (double*)malloc(lengio * sizeof(double));
	matrixInsideOutside(sys, sys->SRowId, sys->SColId, sys->Sval, sys->leng_S, SooRowId, SooColId, Sooval, SoiRowId, SoiColId, Soival, SioRowId, SioColId, Sioval);   // assign each matrix's rowId, colId, val
	/* Transfer Soo to CSR format with SooRowId1 */
	myint* SooRowId1 = (myint*)malloc((sys->outside + 1) * sizeof(myint));
	status = COO2CSR_malloc(SooRowId, SooColId, Sooval, lengoo, sys->outside, SooRowId1);
	/* End of transferring Soo to CSR format with SooRowId1 */
	double *Sooval_old; //dj
	Sooval_old = (double*)malloc(lengoo * sizeof(double));
	for (int ind = 0; ind < lengoo; ++ind) {
		Sooval_old[ind] = Sooval[ind];
	}


	/* calculate D_epsoo+dt^2*Soo */
	//outfile.open("Soo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < lengoo; ++ind) {
	//	//Sooval[ind] = Sooval[ind] * pow(DT, 2);
	//	//if (SooRowId[ind] == SooColId[ind]) {
	//		//Sooval[ind] += sys->getEps(sys->mapEdgeR[sys->mapioR[SooRowId[ind]]]);
	//	//}
	//	outfile << SooRowId[ind] << " " << SooColId[ind] << " ";
	//	outfile << setprecision(15) << Sooval[ind] << endl;
	//}
	//outfile.close();
	/* End of generating Soo, Soi, Sio */

	/* Begin to generate Loo */
	//sparse_matrix_t A;
	//sparse_status_t s0;
	//s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, V0dt, V0da1t, &A);   // A = v0d*v0da'/MU
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

	//outfile.open("Soo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < lengoo; ++ind) {
	//	outfile << SooRowId[ind] + 1 << " " << SooColId[ind] + 1 << " ";
	//	outfile << setprecision(15) << Sooval[ind] << endl;
	//}
	//outfile.close();

	//if (status != 0) {
	//	return status;
	//}
	//myint* LooRowId = NULL;
	//myint* LooRowId1 = NULL;
	//myint* LooColId = NULL;
	//double* Looval = NULL;
	//myint leng_Loo;
	//sparseMatrixSum(sys, A, SooRowId1, SooColId, Sooval, sys->outside, &LooRowId, &LooColId, &Looval, leng_Loo);   //do the sparse matrix summation and put the result in Loo
	//LooRowId1 = (myint*)calloc(sys->outside + 1, sizeof(myint));
	//COO2CSR_malloc(LooRowId, LooColId, Looval, leng_Loo, sys->outside, LooRowId1);
	//double* Looval_old = (double*)calloc(leng_Loo, sizeof(double));
	//for (int ind = 0; ind < leng_Loo; ind++) {
	//	Looval_old[ind] = Looval[ind];
	//}
	//free(Looval); Looval = NULL;
	//sparseMatrixSum1(sys, AoorowId1, AcolId, Aval, SooRowId1, SooColId, Sooval, sys->outside);
	
#endif
	///* output the Loo */
	////outfile.open("Loo.txt", std::ofstream::out | std::ofstream::trunc);
	////for (int ind = 0; ind < leng_Loo; ++ind) {
	////	outfile << LooRowId[ind] + 1 << " " << LooColId[ind] + 1 << " ";
	////	outfile << setprecision(15) << Looval[ind] << endl;
	////	if (Looval[ind] == 0)
	////		cout << "Wrong!\n";
	////}
	////outfile.close();
	////outfile.open("epsoo.txt", std::ofstream::out | std::ofstream::trunc);
	////for (int ind = 0; ind < sys->outside; ++ind) {
	////	outfile << sys->getEps(sys->mapEdgeR[sys->mapioR[ind]]) << endl;
	////}
	////outfile.close();
	///* End of generating Loo */

	///* free the space of V0d*V0da' */
	////free(AoorowId1); AoorowId1 = NULL;
	//mkl_sparse_destroy(A);
	/* end of free of the space V0d*V0da' */



	/* generate Laplacian matrix for this mesh */
	myint nedge = sys->N_edge - sys->bden;
	/* Generate the Laplacian matrix */
	myint leng = 0;
	myint *LrowId, *LcolId;
	double *Lval;
	leng = generateLaplacian_count(sys);   // count how many nnz in the matrix
	LrowId = (myint*)malloc(leng * sizeof(myint));
	LcolId = (myint*)malloc(leng * sizeof(myint));
	Lval = (double*)malloc(leng * sizeof(double));
	status = generateLaplacian(sys, LrowId, LcolId, Lval);   

	outfile.open("L.txt", std::ofstream::out | std::ofstream::trunc);
	for (int ind = 0; ind < leng; ++ind) {
		outfile << LrowId[ind] + 1 << " " << LcolId[ind] + 1 << " ";
		outfile << setprecision(15) << Lval[ind] << endl;
	} 
	outfile.close();
	cout << "L is generated!\n";
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

	/* Print out current */
	out.open("J.txt", std::ofstream::out | std::ofstream::trunc);
	for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
		double* cur = (double*)calloc((sys->N_edge - sys->bden), sizeof(double));   // current
		for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
			for (int indEdge = 0; indEdge < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
				/* Set current density for all edges within sides in port to prepare solver */
				cur[sys->mapEdge[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]]] = sys->portCoor[sourcePort].portDirection[sourcePortSide];
			}
		}
		for (int ind = 0; ind < sys->N_edge - sys->bden; ++ind) {
			out << cur[ind] << endl;
		}
		free(cur); cur = NULL;
	}
	out.close();
	cout << "current is generated!\n";

	
	/* Frequency domain to solve inside outside part */
	for (int ind = 0; ind < sys->nfreq; ++ind) {
#ifdef FREQ
		double freq = sys->freqNo2freq(ind);

		/* generate Loo-omega^2*D_epsoo */
		//double* Looval = (double*)malloc(leng_Loo * sizeof(double));
		//addOmegaEpsilon(sys, LooRowId, LooColId, Looval_old, leng_Loo, sys->outside, freq, Looval);
		/* end of generating Loo-omega^2*D_epsoo */

		/* Frequency domain to solve the inside and outside part
		[-omega^2*D_epsoo+Soo, Soi;
		Sio,                  i*omega*D_sigii+Sii] */
		/* Begin to output (-omega^2*D_epsoo+Soo) */
		//cout << endl << "Frequency oo is " << freq << endl << endl;
		//addOmegaEpsilon(sys, SooRowId, SooColId, Sooval_old, lengoo, sys->outside, freq, Sooval);
		/* End of outputting (-omega^2*D_epsoo+Soo) */

		/* debug testing to see the performance of different solvers */
		//double* bo = (double*)calloc(sys->outside, sizeof(double));
		//double* yo = (double*)calloc(sys->outside, sizeof(double));
		//for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
		//	for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
		//		for (int indEdge = 0; indEdge < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
		//			/* Set current density for all edges within sides in port to prepare solver */
		//			bo[sys->mapio[sys->mapEdge[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]]]] = sys->portCoor[sourcePort].portDirection[sourcePortSide];
		//		}
		//	}

		//	//status = hypreSolve(LooRowId, LooColId, Looval, leng_Loo, bo, sys->outside, yo, 1, 3);   // HYPRE to solve (Loo-omega^2*D_epsoo)
		//	status = hypreSolve(SooRowId, SooColId, Sooval, lengoo, bo, sys->outside, yo, 1, 3);   // HYPRE to solve (Soo-omega^2*D_epsoo)
		//	//status = mkl_gmres_A(bo, yo, LooRowId, LooColId, Looval, leng_Loo, sys->outside);   // gmres to solve (Loo-omega^2*D_epsoo)
		//}
		/* End of debugging to see the performance of different solvers */


		/* Only solve oo part with pardiso */
		//reference_oo(sys, ind, SooRowId1, SooColId, Sooval);
		/* End of only solving oo part with pardiso */

		/* Solve (-omega^2*D_eps+1i*omega*D_sig+S)*xr = -1i*omega*J */
		cout << endl << "Frequency is " << freq << endl << endl;
		reference(sys, ind, sys->SRowId, sys->SColId, sys->Sval);
		/* End of solving (-omega^2*D_eps+1i*omega*D_sig+S)*xr = -1i*omega*J */

		/* start to solve with inside and outside */
		//complex<double>* y;
		//y = (complex<double>*)calloc(nedge * sys->numPorts, sizeof(complex<double>));
		//solveFreqIO(sys, ind, y, V0dt, V0dat, SooRowId, SooRowId1, SooColId, Sooval, lengoo, SioRowId, SioColId, Sioval, lengio, LooRowId1, LooColId, Looval, leng_Loo);
		//for (int sourcePort = 0; sourcePort < sys->numPorts; ++sourcePort) {
		//	sys->Construct_Z_V0_Vh(&y[sourcePort * nedge], ind, sourcePort);
		//}
		//free(y); y = NULL;
		/* end of solving with inside and outside */

		/* Solve [V0a'*D*V0, V0a'*D;
				  D*V0,      D+L] with preconditioner which is the L part */
				  //solveV0L(sys, ind, LrowId, LcolId, Lval, leng, (sys->leng_v0d1 + sys->leng_v0c + sys->N_edge - sys->bden) * 2, V0dt, V0dat, V0ct, V0cat);
				  /* End of solving */

				  //free(Looval); Looval = NULL;

#endif

#ifdef SOLVECHIP
		double freq = sys->freqNo2freq(ind);
		cout << "Frequenc SOLVECHIP is " << freq << endl;
		/* Start to generate matrix [-omega^2*D_eps+L, -omega*D_sig;
		omega*D_sig,      -omega^2*D_eps+L] */
		myint leng_L1 = (leng + sys->inside - sys->outside) * 2;
		myint* L1RowId = (myint*)malloc(leng_L1 * sizeof(myint));
		myint* L1ColId = (myint*)malloc(leng_L1 * sizeof(myint));
		double* L1val = (double*)malloc(leng_L1 * sizeof(double));
		generateL1(sys, freq, LrowId, LcolId, Lval, leng, nedge, L1RowId, L1ColId, L1val, leng_L1);
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
		//myint leng_S1 = (sys->leng_S + sys->inside - sys->outside) * 2;
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

			///* Solve -omega^2*D_eps+1i*omega*D_sig+S by gmres */
			//double* xrri = (double*)calloc(2 * nedge, sizeof(double));
			//double* bri = (double*)calloc(2 * nedge, sizeof(double));
			//for (int ind = 0; ind < nedge; ++ind) {
			//	bri[nedge + ind] = J[ind];
			//}
			//mkl_gmres_A(bri, xrri, S1RowId, S1ColId, S1val, leng_S1, 2 * nedge);
			//complex<double>* xr = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//for (int ind = 0; ind < nedge; ++ind) {
			//	xr[ind] = xrri[ind] + 1i * xrri[ind + nedge];
			//}
			//sys->Construct_Z_V0_Vh(xr, ind, sourcePort);
			///* End of solving -omega^2*D_eps+1i*omega*D_sig+S by gmres */

			/* Solve [V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0, V0a'*(-omega^2*D_eps+1i*omega*D_sig);
					  (-omega^2*D_eps+1i*omega*D_sig)*V0,      -omega^2*D_eps+1i*omega*D_sig+L] */
			double* V0aJ = (double*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(double));   // v0aJ = V0a'*(-omega*J)   question? can whole array put in multiply?
			alpha = 1;
			beta = 0;
			descr.type = SPARSE_MATRIX_TYPE_GENERAL;
			s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, J, beta, V0aJ);
			s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, J, beta, &V0aJ[sys->leng_v0d1]);
			for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
				V0aJ[ind] /= sys->v0dan[ind];
			}
			for (int ind = 0; ind < sys->leng_v0c; ++ind) {
				V0aJ[sys->leng_v0d1 + ind] /= sys->v0can[ind];   // v0aJ= = V0a'*(-omega*J)
			}
			complex<double>* u1 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			solveV0(sys, freq, V0aJ, u1, V0dt, V0dat, V0ct, V0cat, 'i');   // u1 = (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\(V0a'*1i*v0aJ)
			complex<double>* u2 = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			V0Multiply(sys, V0dt, V0ct, sys->leng_v0d1, sys->leng_v0c, nedge, u1, u2);   // u2 = V0*(V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\(V0a'*1i*v0aJ)
			double* u2ri = (double*)calloc(nedge * 2, sizeof(double));
			for (int ind = 0; ind < nedge; ++ind) {
				complex<double> epsig;
				if (sys->markEdge[sys->mapEdgeR[ind]])
					epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * freq * 2 * M_PI * SIGMA;
				else
					epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * 0;
				u2[ind] *= epsig;   // u2 = (-omega^2*D_eps+1i*omega*D_sig)*V0
				u2[ind] = 1i * (J[ind] - u2[ind].imag()) - u2[ind].real();   // u2 = -1i*omega*J-(-omega^2*D_eps+1i*omega*D_sig)*V0*(V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\(V0a'*1i*v0aJ)
				u2ri[ind] = u2[ind].real();   // u2.real()
				u2ri[nedge + ind] = u2[ind].imag();   // u2.imag()
			}
			/* solve [-omega^2*D_eps+L, -omega*D_sig;
					omega*D_sig,      -omega^2*D_eps+L]*[xr; xi] = [u2r; u2i] */
			double* xhri = (double*)calloc(nedge * 2, sizeof(double));
			//hypreSolve(L1RowId, L1ColId, L1val, leng_L1, u2ri, nedge * 2, xhri, 1, 4);
			mkl_gmres_A(u2ri, xhri, L1RowId, L1ColId, L1val, leng_L1, 2 * nedge);
			complex<double>* xh = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			for (int ind = 0; ind < nedge; ++ind) {
				xh[ind] = xhri[ind] + 1i * xhri[nedge + ind];   // xh is got
			}

			/* (-omega^2*D_eps+1i*omega*D_sig)*xh */
			for (int ind = 0; ind < nedge; ++ind) {
				complex<double> epsig;
				if (sys->markEdge[sys->mapEdgeR[ind]])
					epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * freq * 2 * M_PI * SIGMA;
				else
					epsig = -pow(freq * 2 * M_PI, 2) * sys->getEps(sys->mapEdgeR[ind]) + 1i * 0;
				u2[ind] = xh[ind] * epsig;   // u2 = (-omega^2*D_eps+1i*omega*D_sig)*xh
				u2[ind] = 1i * (J[ind] - u2[ind].imag()) - u2[ind].real();   // u2 = -1i*omega*J-(-omega^2*D_eps+1i*omega*D_sig)*xh
			}
			V0aMultiply(sys, V0dat, V0cat, sys->leng_v0d1, sys->leng_v0c, nedge, u2, u1);   // u1 = V0a'*(-1i*omega*J-(-omega^2*D_eps+1i*omega*D_sig)*xh)
			double* tempr = (double*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(double));
			double* tempi = (double*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(double));
			complex<double>* u01 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			complex<double>* u02 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			for (int ind = 0; ind < sys->leng_v0d1 + sys->leng_v0c; ++ind) {
				tempr[ind] = u1[ind].real();
				tempi[ind] = u1[ind].imag();
			}
			solveV0(sys, freq, tempr, u01, V0dt, V0dat, V0ct, V0cat, 'r');   // u01 = (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\real(V0a'*(1i*v0aJ-(-omega^2*D_eps+1i*omega*D_sig)*xh))
			solveV0(sys, freq, tempi, u02, V0dt, V0dat, V0ct, V0cat, 'i');   // u01 = (V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0)\imag(V0a'*(1i*v0aJ-(-omega^2*D_eps+1i*omega*D_sig)*xh))
			complex<double>* u0 = (complex<double>*)calloc(sys->leng_v0d1 + sys->leng_v0c, sizeof(complex<double>));
			for (int ind = 0; ind < sys->leng_v0d1 + sys->leng_v0c; ++ind) {
				u0[ind] = u01[ind].real() + u02[ind].real() + 1i * (u01[ind].imag() + u02[ind].imag());
			}
			complex<double>* x0 = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			complex<double>* x = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			V0Multiply(sys, V0dt, V0ct, sys->leng_v0d1, sys->leng_v0c, nedge, u0, x0);   // x0 = V0*u0
			for (int ind = 0; ind < nedge; ++ind) {
				x[ind] = x0[ind].real() + xh[ind].real() + 1i * (x0[ind].imag() + xh[ind].imag());   // x = V0*u0+xh
			}


			sys->Construct_Z_V0_Vh(x, ind, sourcePort);

			///* Calculate the error norm */
			//xr = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			//sys->reference1(ind, sourcePort, xr);
			//double error = 0, b = 0;
			//for (int ind = 0; ind < nedge; ++ind) {
			//	error += pow(xr[ind].real() - x[ind].real(), 2) + pow(xr[ind].imag() - x[ind].imag(), 2);
			//	b += pow(xr[ind].real(), 2) + pow(xr[ind].imag(), 2);
			//}
			//cout << "Frequency " << freq << " error is " << error / b << endl;

			//free(xr); xr = NULL;
			///* End of calculating the error norm */

			free(J); J = NULL;
			free(V0aJ); V0aJ = NULL;
			free(u1); u1 = NULL;
			free(u2); u2 = NULL;
			free(u2ri); u2ri = NULL;
			free(xhri); xhri = NULL;
			free(xh); xh = NULL;
			free(tempr); tempr = NULL;
			free(tempi); tempi = NULL;
			free(u01); u01 = NULL;
			free(u02); u02 = NULL;
			free(u0); u0 = NULL;
			free(x0); x0 = NULL;
			free(x); x = NULL;
			/* End of solving [V0a'*(-omega^2*D_eps+1i*omega*D_sig)*V0, V0a'*(-omega^2*D_eps+1i*omega*D_sig);
			(-omega^2*D_eps+1i*omega*D_sig)*V0,      -omega^2*D_eps+1i*omega*D_sig+L] */
		}
		free(L1RowId); L1RowId = NULL;
		free(L1ColId); L1ColId = NULL;
		free(L1val); L1val = NULL;

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
    mkl_sparse_destroy(V0ct);
    mkl_sparse_destroy(V0cat);
	
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
#endif

    //HYPRE_IJMatrixDestroy(ad);
    //HYPRE_IJMatrixDestroy(ac);

    return 0;
}




int COO2CSR(vector<int> &rowId, vector<int> &ColId, vector<double> &val) {
    int i;
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

void matrixInsideOutside_count(fdtdMesh* sys, myint* rowId, myint* colId, myint leng, myint& lengoo, myint& lengoi, myint& lengio) {
	/* count each submatrix nnz number
	   sys : used to provide inside outside edge numbers
	   rowId : matrix rowId
	   colId : matrix colId
	   leng: : nnz in this matrix
	   lengoo : oo nnz
	   lengoi : oi nnz
	   lengio : io nnz */

	for (myint ind = 0; ind < leng; ++ind) {
		if (sys->mapio[rowId[ind]] >= 0 && sys->mapio[rowId[ind]] < sys->outside && sys->mapio[colId[ind]] >= 0 && sys->mapio[colId[ind]] < sys->outside) {
			lengoo++;
		}
		else if (sys->mapio[rowId[ind]] >= 0 && sys->mapio[rowId[ind]] < sys->outside && sys->mapio[colId[ind]] >= sys->outside && sys->mapio[colId[ind]] < sys->inside) {
			lengoi++;
		}
		else if (sys->mapio[rowId[ind]] >= sys->outside && sys->mapio[rowId[ind]] < sys->inside && sys->mapio[colId[ind]] >= 0 && sys->mapio[colId[ind]] < sys->outside) {
			lengio++;
		}
	}
	return;
}

void matrixInsideOutside(fdtdMesh* sys, myint* rowId, myint* colId, double* val, myint leng, myint* ooRowId, myint* ooColId, double* ooval, myint* oiRowId, myint* oiColId, double* oival, myint* ioRowId, myint* ioColId, double* ioval) {
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
		if (sys->mapio[rowId[ind]] >= 0 && sys->mapio[rowId[ind]] < sys->outside && sys->mapio[colId[ind]] >= 0 && sys->mapio[colId[ind]] < sys->outside) {
			ooRowId[lengoo] = sys->mapio[rowId[ind]];
			ooColId[lengoo] = sys->mapio[colId[ind]];
			ooval[lengoo] = val[ind];
			lengoo++;
		}
		else if (sys->mapio[rowId[ind]] >= 0 && sys->mapio[rowId[ind]] < sys->outside && sys->mapio[colId[ind]] >= sys->outside && sys->mapio[colId[ind]] < sys->inside) {
			oiRowId[lengoi] = sys->mapio[rowId[ind]];
			oiColId[lengoi] = sys->mapio[colId[ind]] - sys->outside;
			oival[lengoi] = val[ind];
			lengoi++;
		}
		else if (sys->mapio[rowId[ind]] >= sys->outside && sys->mapio[rowId[ind]] < sys->inside && sys->mapio[colId[ind]] >= 0 && sys->mapio[colId[ind]] < sys->outside) {
			ioRowId[lengio] = sys->mapio[rowId[ind]] - sys->outside;
			ioColId[lengio] = sys->mapio[colId[ind]];
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

void solveV0(fdtdMesh* sys, double freq, double* rhs, complex<double>* u0, sparse_matrix_t& v0dt, sparse_matrix_t& v0dat, sparse_matrix_t& v0ct, sparse_matrix_t& v0cat, char ri){
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
	myint nedge = sys->N_edge - sys->bden;

	double* temp = (double*)calloc(sys->leng_v0d1, sizeof(double));
	double* temp1 = (double*)calloc(nedge, sizeof(double));
	double* temp2 = (double*)calloc(sys->leng_v0c, sizeof(double));
	double* u0c = (double*)calloc(sys->leng_v0c, sizeof(double));
	double* u0di = (double*)calloc(sys->leng_v0d1, sizeof(double));
	double* u0dr = (double*)calloc(sys->leng_v0d1, sizeof(double));
	status = hypreSolve(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, rhs, sys->leng_v0d1, temp, 0, 3);   // temp = (V0da'*D_eps*V0d)\bd
	for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
		temp[ind] /= sys->v0dn[ind];
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, v0dt, descr, temp, beta, temp1);   // temp1 = V0d*(V0da'*D_eps*V0d)\bd
	for (int ind = 0; ind < nedge; ++ind) {
		temp1[ind] *= sys->getEps(sys->mapEdgeR[ind]);   // temp1 = D_eps*V0d*(V0da'*D_eps*V0d)\bd
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, v0cat, descr, temp1, beta, temp2);   // temp2 = V0ca'*D_eps*V0d*(V0da'*D_eps*V0d)\bd
	for (int ind = 0; ind < sys->leng_v0c; ++ind) {
		temp2[ind] /= sys->v0can[ind];
		temp2[ind] = rhs[sys->leng_v0d1 + ind] - temp2[ind];   // temp2 = bc-V0ca'*D_eps*V0d*(V0da'*D_eps*V0d)\bd
	}

	for (int indi = 0; indi < sys->leng_v0c; ++indi) {   // Ac is not normalized with V0ca and V0c
		temp2[indi] *= sys->v0can[indi];
	}
	status = hypreSolve(sys->AcRowId, sys->AcColId, sys->Acval, sys->leng_Ac, temp2, sys->leng_v0c, u0c, 0, 3);
	for (int indi = 0; indi < sys->leng_v0c; ++indi) {
		u0c[indi] /= -((freq * 2 * M_PI) / sys->v0cn[indi]);    // u0c = (V0ca'*(omega*D_sig)*V0c)(bc-V0ca'*D_eps*V0d*(V0da'*D_eps*V0d)\bd), imaginary
		temp2[indi] = u0c[indi];
	}

	for (int ind = 0; ind < sys->leng_v0c; ++ind) {
		temp2[ind] /= sys->v0cn[ind];
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, v0ct, descr, temp2, beta, temp1);   // temp1 = V0c*u0c, imaginary
	for (int ind = 0; ind < nedge; ++ind) {
		temp1[ind] *= sys->getEps(sys->mapEdgeR[ind]);   // temp1 = omega^2*D_eps*V0c*u0c
	}
	s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, v0dat, descr, temp1, beta, temp);
	for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
		temp[ind] /= sys->v0dan[ind];
	}
	status = hypreSolve(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, temp, sys->leng_v0d1, u0di, 0, 3);   // u0di = (V0da'*D_eps*V0d)\(V0da'*D_eps*V0c*u0c), imaginary
	for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
		u0di[ind] *= -1;   // u0di = -(V0da'*D_eps*V0d)\(V0da'*D_eps*V0c*u0c), imaginary
	}
	status = hypreSolve(sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, rhs, sys->leng_v0d1, u0dr, 0, 3);   // u0dr = (V0da'*D_eps*V0d)\(bd), real
	for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
		u0dr[ind] /= (-pow(freq * 2 * M_PI, 2));   // u0dr = (V0da'*(-omega^2*D_eps)*V0d)\(bd), real
	}

	/* Assign to the final solution */
	if (ri == 'r') {
		for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
			u0[ind] = u0dr[ind] + 1i * u0di[ind];
		}
		for (int ind = 0; ind < sys->leng_v0c; ++ind) {
			u0[sys->leng_v0d1 + ind] = 0 + 1i * u0c[ind];
		}
	}
	else if (ri == 'i') {
		for (int ind = 0; ind < sys->leng_v0d1; ++ind) {
			u0[ind] = 1i * u0dr[ind] - u0di[ind];
		}
		for (int ind = 0; ind < sys->leng_v0c; ++ind) {
			u0[sys->leng_v0d1 + ind] = - u0c[ind];
		}
	}

	free(temp); temp = NULL;
	free(temp1); temp1 = NULL;
	free(temp2); temp2 = NULL;
	free(u0c); u0c = NULL;
	free(u0di); u0di = NULL;
	free(u0dr); u0dr = NULL;
}

void V0Multiply(fdtdMesh* sys, sparse_matrix_t& V0dt, sparse_matrix_t& V0ct, myint lengv0d, myint lengv0c, myint row, complex<double>* u1, complex<double>* u2) {
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
		tempr[ind] = u1[ind].real() / sys->v0dn[ind];
		tempi[ind] = u1[ind].imag() / sys->v0dn[ind];
	}
	for (int ind = lengv0d; ind < lengv0d + lengv0c; ++ind) {
		tempr[ind] = u1[ind].real() / sys->v0cn[ind - lengv0d];
		tempi[ind] = u1[ind].imag() / sys->v0cn[ind - lengv0d];
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

