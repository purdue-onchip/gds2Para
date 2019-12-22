//#include "stdafx.h"
#include <ctime>
#include "fdtd.hpp"
#include "hypreSolver.h"


static bool comp(pair<double, int> a, pair<double, int> b) {
    return a.first <= b.first;
};

int paraGenerator(fdtdMesh *sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi) {


    int indi, indj, mark, k, l, n;
    int status = 0;
    int count = 0;
    int xcol = 0;
    vector<int> rowId;
    vector<int> colId;
    vector<double> val;
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
	/* Begin to set v0d to be the outside v0d */
	for (indi = 0; indi < sys->v0d1num; indi++) {
		sys->v0d1RowId[indi] = sys->mapio[sys->mapEdge[sys->v0d1RowId[indi]]];
	}
	/* End of setting v0d to be the outside v0d */

    /*sys->v0d1aColIdo = (myint*)malloc(sys->v0d1anum * sizeof(myint));
    for (indi = 0; indi < sys->v0d1anum; indi++)
    sys->v0d1aColIdo[indi] = sys->v0d1aColId[indi];
    free(sys->v0d1aColId); sys->v0d1aColId = (myint*)malloc((sys->sys->leng_v0d1a + 1) * sizeof(myint));
    status = COO2CSR_malloc(sys->v0d1aColIdo, sys->v0d1aRowId, sys->v0d1aval, sys->v0d1anum, sys->sys->leng_v0d1a, sys->v0d1aColId);
    if (status != 0)
    return status;*/

    //cout << "Number of NNZ in V0d1 is " << sys->v0d1num << endl;

    sparse_status_t s;

    /* V0d^T's csr form handle for MKL */
    sparse_matrix_t V0dt;
    s = mkl_sparse_d_create_csr(&V0dt, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, sys->outside/*sys->N_edge*/, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1valo);

    /* V0da^T's csr form handle for MKL */
    sparse_matrix_t V0dat;
    s = mkl_sparse_d_create_csr(&V0dat, SPARSE_INDEX_BASE_ZERO, sys->leng_v0d1, sys->outside/*sys->N_edge*/, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1avalo);

#ifdef V0_new_schema
    double * drhs = new double[sys->leng_v0d1 * leng_v0d2];
    mkl_sparse_d_mm(SPARSE_OPERARION_NON_TRANSPOSE, 1, V0dat, SPARSE_MATRIX_TYPE_GENERAL, SPARSE_LAYOUT_COLUMN_MAJOR, temp, leng_v0d2, (sys->N_edge - sys->bden), 0, drhs, sys->leng_v0d1);
    delete[] temp;
    delete[] drhs;
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

#ifdef PRINT_VERBOSE_TIMING
    cout << "Time to generate V0c is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
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
    cout << "Time to generate Ac is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
    cout << "Number of non-zeros in Ac is " << sys->leng_Ac << endl;
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
	//	outfile << sys->v0d1RowId[indi] + 1 << " " << sys->v0d1ColIdo[indi] + 1 << " " << sys->v0d1valo[indi] << endl;
	//}
	//outfile.close();

	//outfile.open("V0da.txt", std::ofstream::trunc | std::ofstream::out);
	//for (indi = 0; indi < sys->v0d1num; ++indi) {
	//	outfile << sys->v0d1RowId[indi] + 1 << " " << sys->v0d1ColIdo[indi] + 1 << " " << sys->v0d1avalo[indi] << endl;
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
    ofstream out;

	/* Begin to map the inside and outside edges */
	sys->mapEdgeInsideOutside();   // generate mapio and mapioR
	/* End of mapping the inside and outside edges */

	/* Begin to generate Soo, Soi, Sio */
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
	/* calculate D_epsoo+dt^2*Soo */
	//outfile.open("Moo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < lengoo; ++ind) {
	//	Sooval[ind] = Sooval[ind] * pow(DT, 2);
	//	if (SooRowId[ind] == SooColId[ind]) {
	//		Sooval[ind] += sys->getEps(sys->mapEdgeR[sys->mapioR[SooRowId[ind]]]);
	//	}
	////	outfile << SooRowId[ind] << " " << SooColId[ind] << " " << Sooval[ind] << endl;
	//}
	//outfile.close();
	/* End of generating Soo, Soi, Sio */

	/* Begin to generate Loo */
	myint* AoorowId, *AoocolId;
	double* Aooval;   // the operator gradient divergence
	sparse_matrix_t A;
	sparse_status_t s0;
	s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, V0dt, V0dat, &A);   // A = v0d*v0da'
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	myint ARows, ACols, leng_Aoo = 0;
	MKL_INT *ArowStart, *ArowEnd, *AcolId;
	double *Aval;
	s0 = mkl_sparse_d_export_csr(A, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
	myint leng_A = ArowEnd[ARows - 1];    // how many non-zeros are in A
	matrixOutside_count(sys, ArowStart, ArowEnd, AcolId, Aval, ARows, leng_Aoo);    // count how many nnz in Aoo
	AoorowId = (myint*)malloc(leng_Aoo * sizeof(myint));
	AoocolId = (myint*)malloc(leng_Aoo * sizeof(myint));
	Aooval = (double*)malloc(leng_Aoo * sizeof(double));
	matrixOutside(sys, ArowStart, ArowEnd, AcolId, Aval, ARows, AoorowId, AoocolId, Aooval, 1 / MU);   // assign Aoo rowId, colId, val
	myint* AoorowId1 = (myint*)malloc((sys->outside + 1) * sizeof(myint));
	status = COO2CSR_malloc(AoorowId, AoocolId, Aooval, leng_Aoo, sys->outside, AoorowId1);
	//outfile.open("Aoo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < leng_Aoo; ++ind) {
	//	outfile << AoorowId[ind] << " " << AoocolId[ind] << " " << Aooval[ind] << endl;
	//}
	//outfile.close();

	myint* SooRowId1 = (myint*)malloc((sys->outside + 1) * sizeof(myint));
	status = COO2CSR_malloc(SooRowId, SooColId, Sooval, lengoo, sys->outside, SooRowId1);
	//outfile.open("Soo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < lengoo; ++ind) {
	//	outfile << SooRowId[ind] << " " << SooColId[ind] << " " << Sooval[ind] << endl;
	//}
	//outfile.close();

	if (status != 0) {
		return status;
	}
	sparseMatrixSum(sys, AoorowId1, AoocolId, Aooval, SooRowId1, SooColId, Sooval, sys->outside);   //do the sparse matrix summation and put the result in Loo
	/* output the Loo and Soo */
	//outfile.open("Loo.txt", std::ofstream::out | std::ofstream::trunc);
	//for (int ind = 0; ind < sys->leng_Loo; ++ind) {
	//	outfile << sys->LooRowId[ind] << " " << sys->LooColId[ind] << " " << sys->Looval[ind] << endl;
	//}
	//outfile.close();
	/* End of generating Loo */
	myint nedge = sys->N_edge - sys->bden;
	/* HYPRE solves for each port are messy */
    for (sourcePort = 0; sourcePort < sys->numPorts; sourcePort++) {
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

		/* debug testing to see the performance of different solvers */
		double* bo = (double*)calloc(sys->outside, sizeof(double));
		double* yo = (double*)calloc(sys->outside, sizeof(double));
		for (int ind = 0; ind < sys->N_edge; ++ind) {
			if (sys->markEdge[ind] == 0) {
				bo[sys->mapio[ind]] = 1;   // test right hand side
			}
		}
		//status = hypreSolve(sys, sys->LooRowId, sys->LooColId, sys->Looval, sys->leng_Loo, bo, sys->outside, yo, 1, 3);   // HYPRE to solve (dt^2*Loo+D_epsoo)
		//status = hypreSolve(sys, SooRowId, SooColId, Sooval, lengoo, bo, sys->outside, yo, 1, 3);   // HYPRE to solve (dt^2*Soo+D_epsoo)
		status = mkl_gmres_A(sys, bo, yo, sys->LooRowId, sys->LooColId, sys->Looval, sys->leng_Loo, sys->outside);   // gmres to solve (dt^2*Loo+D_epsoo)
		/* End of debugging to see the performance of different solvers */
		
		/* backward difference */
		//find_Vh_back(sys, sourcePort, V0ct, V0cat, V0dt, V0dat);
		

		/* central difference */
		//find_Vh_central(sys, sourcePort);

		/* Frequency domain to solve the inside and outside part
		[-omega^2*D_epsoo+Soo, Soi;
		 Sio,                  i*omega*D_sigii+Sii] */
		complex<double>* y;
		for (int ind = 0; ind < sys->nfreq; ++ind) {
			y = (complex<double>*)calloc(nedge, sizeof(complex<double>));
			solveFreqIO(sys, sourcePort, ind, y, V0dt, V0dat, SioRowId, SioColId, Sioval, lengio);
			free(y);
		}


#ifdef GENERATE_V0_SOLUTION


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

        //status = hypreSolve(sys, ad, parcsr_ad, sys->leng_Ad, v0daJ, sys->leng_v0d1, y0d);
        status = hypreSolve(sys, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, v0daJ, sys->leng_v0d1, y0d, 1, 3);
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

            //status = hypreSolve(sys, ia, ja, a, (sys->cindex[indi] - sys->cindex[indi - 1]), crhss, (sys->acu_cnno[indi] - sys->acu_cnno[indi - 1]), y0cs);

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

            //status = hypreSolve(sys, ac, parcsr_ac, sys->leng_Ac, v0caJ, sys->leng_v0c, y0c);
            status = hypreSolve(sys, sys->AcRowId, sys->AcColId, sys->Acval, sys->leng_Ac, v0caJ, sys->leng_v0c, y0c, 1, 3);

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

        
        status = hypreSolve(sys, sys->AdRowId, sys->AdColId, sys->Adval, sys->leng_Ad, dRhs2, sys->leng_v0d1, y0d2, 1, 3);
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
        //sys->Construct_Z_V0(yd, sourcePort);
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

        for (indi = 0; indi < sys->nfreq; indi++){
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
            for (myint j = 0; j < sys->leng_Vh; j++){    // calculate (A+C)*V_re1
                inde = 0;
                while (inde < sys->leng_S){
                    start = sys->SRowId[inde];
                    while (inde < sys->leng_S && sys->SRowId[inde] == start){

                        tmp[j * (sys->N_edge - sys->bden) + sys->SRowId[inde]].real += sys->Sval[inde] * Vh[j * (sys->N_edge - sys->bden) + sys->SColId[inde]].real;
                        tmp[j * (sys->N_edge - sys->bden) + sys->SRowId[inde]].imag += sys->Sval[inde] * Vh[j * (sys->N_edge - sys->bden) + sys->SColId[inde]].imag;


                        inde++;
                    }
                }

                for (inde = 0; inde < sys->N_edge - sys->bden; inde++){
                    if (sys->markEdge[sys->mapEdgeR[inde]] != 0){
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
            for (inde = 0; inde < sys->N_edge; inde++){
                if (sys->J[inde] != 0) {
                    J[sys->mapEdge[inde]].imag = -1 * freq * 2 * M_PI;
                }
            }
            
            status = matrix_multi('T', Vh, (sys->N_edge - sys->bden), sys->leng_Vh, J, (sys->N_edge - sys->bden), 1, rhs_h);    // -1i*omega*V_re1'*J
            free(J); J = NULL;
            

            /* V_re1'*A*u */
            free(tmp);
            tmp = (lapack_complex_double*)calloc((sys->N_edge - sys->bden), sizeof(lapack_complex_double));
            for (inde = 0; inde < sys->N_edge - sys->bden; inde++){
                if (sys->markEdge[sys->mapEdgeR[inde]] != 0){
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
            for (inde = 0; inde < sys->leng_Vh; inde++){
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
            for (inde = 0; inde < sys->N_edge - sys->bden; inde++){
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
            for (inde = 0; inde < sys->N_edge - sys->bden; inde++){
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
#endif

		free(sys->J); sys->J = NULL;
        /*free(u0); u0 = NULL;*/
        //free(sys->y); sys->y = NULL;
        xcol++;
    }
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
    sys->print_z_V0();
#endif

#ifdef PRINT_V0_Vh_Z_PARAM
    sys->print_z_V0_Vh();
#endif
    sys->cindex.clear();
    sys->acu_cnno.clear();
    sys->stackEpsn.clear();
    
    mkl_sparse_destroy(V0dt);
    mkl_sparse_destroy(V0dat);
    mkl_sparse_destroy(V0ct);
    mkl_sparse_destroy(V0cat);
	//mkl_sparse_destroy(A);

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

	for (int ind = 0; ind < leng; ++ind) {
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
	for (int ind = 0; ind < leng; ++ind) {
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
	for (int ind = 0; ind < ARows; ++ind) {
		if (sys->mapio[ind] < sys->outside) {
			for (int indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
				if (sys->mapio[AcolId[indi]] >= 0 && sys->mapio[AcolId[indi]] < sys->outside) {
					leng_Aoo++;
				}
			}
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
	for (int ind = 0; ind < ARows; ++ind) {
		if (sys->mapio[ind] < sys->outside) {
			for (int indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
				if (sys->mapio[AcolId[indi]] >= 0 && sys->mapio[AcolId[indi]] < sys->outside) {
					AoorowId[leng_Aoo] = sys->mapio[ind];
					AoocolId[leng_Aoo] = sys->mapio[AcolId[indi]];
					Aooval[leng_Aoo] = Aval[indi] * scale;
					leng_Aoo++;
				}
			}
		}
	}
}

void sparseMatrixSum(fdtdMesh* sys, myint* arowId1, myint* acolId, double* aval, myint* browId1, myint* bcolId, double* bval, myint Rows) {
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

	s0 = mkl_sparse_d_add(operation, b, alpha, a, &r);   // add the two sparse matrices

	myint ARows, ACols;
	MKL_INT *ArowStart, *ArowEnd;
	MKL_INT *AcolId;
	double *Aval;
	sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
	s0 = mkl_sparse_d_export_csr(r, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
	sys->leng_Loo = ArowEnd[ARows - 1];

	sys->LooRowId = (myint*)malloc(sys->leng_Loo * sizeof(myint));
	sys->LooColId = (myint*)malloc(sys->leng_Loo * sizeof(myint));
	sys->Looval = (double*)malloc(sys->leng_Loo * sizeof(double));

	myint count = 0;
	for (int ind = 0; ind < ARows; ++ind) {
		for (int indi = ArowStart[ind]; indi < ArowEnd[ind]; ++indi) {
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