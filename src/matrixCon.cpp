//#include "stdafx.h"
#include <ctime>
#include "fdtd.hpp"
#include "hypreSolver.h"


static bool comp(pair<double, int> a, pair<double, int> b) {
    return a.first <= b.first;
};

int paraGenerator(fdtdMesh *sys, tf::Subflow subflow) {
    tf::Task taskPG0 = subflow.placeholder();
    clock_t t0 = clock();
    double *temp;
    int bdl = 0, bdu = 0;   // # of boundary
#ifdef LOWER_BOUNDARY_PEC
    bdl = 1;
#endif
#ifdef UPPER_BOUNDARY_PEC
    bdu = 1;
#endif

    /* Construct V0d with row id, col id and its val */
    //shared_ptr<myint> leng_v0d1(new myint);
    myint *leng_v0d1 = (myint*)calloc(1, sizeof(myint));
    myint *v0d1num = (myint*)calloc(1, sizeof(myint)); // Number of v0d1 vectors, which are nodes outside the conductors
    myint *leng_v0d1a = (myint*)calloc(1, sizeof(myint));
    myint *v0d1anum = (myint*)calloc(1, sizeof(myint));
    myint *leng_v0d2 = (myint*)calloc(1, sizeof(myint));
    myint *leng_v0d2a = (myint*)calloc(1, sizeof(myint));
    myint *v0d2num = (myint*)calloc(1, sizeof(myint));
    myint *v0d2anum = (myint*)calloc(1, sizeof(myint));
    myint *leng_Ad = (myint*)calloc(1, sizeof(myint));
    myint *mapd = (myint*)calloc(sys->N_node, sizeof(myint));
    myint *mapc = (myint*)calloc(sys->N_node, sizeof(myint));
    double block1_x, block1_y, block2_x, block2_y, block3_x, block3_y;
    double sideLen = 0.; // Around the conductor, 10 um is considered with rapid potential change

    block1_x = 0.;// (sys->xlim2 - sys->xlim1) / 5 * sys->lengthUnit;
    block1_y = 0.;// (sys->ylim2 - sys->ylim1) / 5 * sys->lengthUnit;
    block2_x = 0.;// (sys->xlim2 - sys->xlim1) / 5 * sys->lengthUnit;:q
    block2_y = 0.;// (sys->ylim2 - sys->ylim1) / 5 * sys->lengthUnit;
    block3_x = 0.;
    block3_y = 0.;
    taskPG0.name("Initialize Vars for paraGenerator()");
    /*sparse_matrix_t mat1;
    myint *rowin = (myint*)calloc(5, sizeof(myint));
    myint *colin = (myint*)calloc(5, sizeof(myint));
    double *valin = (double*)calloc(5, sizeof(double));
    rowin[0] = 0;
    colin[0] = 0;
    valin[0] = 2.2;
    rowin[1] = 0;
    colin[1] = 1;
    valin[1] = 7.3;
    rowin[2] = 1;
    colin[2] = 1;
    valin[2] = 4.1;
    rowin[3] = 2;
    colin[3] = 2;
    valin[3] = 5.9;
    sparse_status_t stat = mkl_sparse_d_create_csr(&mat1, SPARSE_INDEX_BASE_ZERO, 3, 3, &rowin[0], &rowin[1], colin, valin);
    free(rowin);
    free(colin);
    free(valin);*/
    /*sparse_index_base_t indexing;
    myint rows = 0;
    myint cols = 0;
    myint *rows_start, *rows_end, *col_indx;
    double *values;
    sparse_status_t check = mkl_sparse_d_export_csr(mat1, &indexing, &rows, &cols, &rows_start, &rows_end, &col_indx, &values);
    cout << "check = " << check << endl;
    cout << "MKL Status of mat1: rows = " << rows << " and cols = " << cols << endl;
    cout << "MKL Status of mat1: rows_start[rows] = " << rows_start[rows] << " and col_indx[rows - 1] = " << col_indx[rows - 1] << endl;
    cout << "MKL Status of mat1: values[0] = " << values[0] << endl;*/


    /* Generate V0d */
    tf::Task taskPG1 = subflow.emplace([=]() {
        cout << "Start generating V0d" << endl;
#ifdef PRINT_V0D_BLOCKS
        cout << " V0d's block1_x and block1_y are " << block1_x << " " << block1_y << endl;
        cout << " V0d's block2_x and block2_y are " << block2_x << " " << block2_y << endl;
        cout << " V0d's block3_x and block3_y are " << block3_x << " " << block3_y << endl;
#endif

        sys->merge_v0d1(block1_x, block1_y, block2_x, block2_y, block3_x, block3_y, *v0d1num, *leng_v0d1, *v0d1anum, *leng_v0d1a, mapd, sideLen, bdl, bdu);

#ifdef PRINT_VERBOSE_TIMING
        cout << " Time to generate V0d is " << (clock() - t0) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif
        cout << " Length of V0d1 is " << *leng_v0d1 << ", and number of non-zeros in V0d1 is " << *v0d1num << endl;
        cout << " Length of V0d1a is " << *leng_v0d1a << ", and number of non-zeros in V0d1a is " << *v0d1anum << endl;
        cout << "V0d is generated!" << endl;
    });
    taskPG0.precede(taskPG1);
    taskPG1.name("Generate V_{0d}");

    /* Generate Ad = V0da^T*D_eps*V0d = D_eps,0 */
    tf::Task taskPG2 = subflow.emplace([=]() {
        clock_t t1 = clock();
        sys->generateAd(mapd, *v0d1num, *v0d1anum, *leng_v0d1, *leng_Ad);

#ifdef PRINT_VERBOSE_TIMING
        cout << "Time to generate Ad is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
        cout << "Number of non-zeros in Ad is " << *leng_Ad << endl;
#endif

        sys->v0d1valo = (double*)malloc(*v0d1num * sizeof(double));
        for (myint indi = 0; indi < *v0d1num; indi++) {
            sys->v0d1valo[indi] = sys->v0d1val[indi];
        }
        for (myint indi = 0; indi < *v0d1num; indi++) {    // compute sqrt(D_eps)*V0d1
            sys->v0d1val[indi] *= sqrt(sys->stackEpsn[(sys->v0d1RowId[indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
        }
        sys->v0d1avalo = (double*)malloc(*v0d1anum * sizeof(double));
        for (myint indi = 0; indi < *v0d1anum; indi++) {
            sys->v0d1avalo[indi] = sys->v0d1aval[indi];
        }
        for (myint indi = 0; indi < *v0d1anum; indi++) {
            sys->v0d1aval[indi] *= sqrt(sys->stackEpsn[(sys->v0d1RowId[indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
        }


        myint leng_v0d = *leng_v0d1;
        sys->v0d1ColIdo = (myint*)malloc(*v0d1num * sizeof(myint));
        for (myint indi = 0; indi < *v0d1num; indi++) {
            sys->v0d1ColIdo[indi] = sys->v0d1ColId[indi];
        }
        free(sys->v0d1ColId);
        sys->v0d1ColId = (myint*)malloc((*leng_v0d1 + 1) * sizeof(myint));
        int status = COO2CSR_malloc(sys->v0d1ColIdo, sys->v0d1RowId, sys->v0d1val, *v0d1num, *leng_v0d1, sys->v0d1ColId);
        if (status != 0) {
            cerr << "Could not make V0d1ColIdo into a sparse matrix (status = " << status << ")" << endl;
            exit(1);
        }
        free(sys->v0d1val); sys->v0d1val = NULL;
        free(sys->v0d1ColIdo); sys->v0d1ColIdo = NULL;

        /*sys->v0d1aColIdo = (myint*)malloc(v0d1anum * sizeof(myint));
        for (indi = 0; indi < v0d1anum; indi++) {
            sys->v0d1aColIdo[indi] = sys->v0d1aColId[indi];
        }
        free(sys->v0d1aColId); sys->v0d1aColId = (myint*)malloc((leng_v0d1a + 1) * sizeof(myint));
        status = COO2CSR_malloc(sys->v0d1aColIdo, sys->v0d1aRowId, sys->v0d1aval, v0d1anum, leng_v0d1a, sys->v0d1aColId);
        if (status != 0) {
            cerr << "Could not make V0d1aColIdo into a sparse matrix (status = " << status << ")" << endl;
            exit(1);
        }*/
        free(sys->v0d1aval); sys->v0d1aval = NULL;
    });
    taskPG2.name("Generate A_{d}");
    taskPG1.precede(taskPG2);


    /* V0d^T's csr form handle for MKL */
    sparse_matrix_t *V0dt = (sparse_matrix_t*)calloc(1, sizeof(sparse_matrix_t));
    tf::Task taskPG3 = subflow.emplace([=]() {
        sparse_status_t s = mkl_sparse_d_create_csr(V0dt, SPARSE_INDEX_BASE_ZERO, *leng_v0d1, sys->N_edge, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1valo);
        if (s != 0) {
            cerr << "Could not create sparse matrix V0dt (status = " << s << ")" << endl;
        }
    });
    taskPG3.name("V_{0d}^T to CSR Form");
    taskPG2.precede(taskPG3);

    /* V0da^T's csr form handle for MKL */
    sparse_matrix_t *V0dat = (sparse_matrix_t*)calloc(1, sizeof(sparse_matrix_t));
    tf::Task taskPG4 = subflow.emplace([=]() {
        sparse_status_t s = mkl_sparse_d_create_csr(V0dat, SPARSE_INDEX_BASE_ZERO, *leng_v0d1, sys->N_edge, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1avalo);
        if (s != 0) {
            cerr << "Could not create sparse matrix V0dat (status = " << s << ")" << endl;
        }
    });
    taskPG4.name("V_{0da}^T to CSR Form");
    taskPG2.precede(taskPG4);


    /* Construct V0c with row id, col id and its val */
    myint *leng_v0c = (myint*)calloc(1, sizeof(myint));
    myint *v0cnum = (myint*)calloc(1, sizeof(myint));
    myint *leng_v0ca = (myint*)calloc(1, sizeof(myint));
    myint *v0canum = (myint*)calloc(1, sizeof(myint));
    myint *leng_Ac = (myint*)calloc(1, sizeof(myint));

    block1_x = 0;// (sys->xlim2 - sys->xlim1) / 3 * sys->lengthUnit;// (sys->xn[sys->nx - 1] - sys->xn[0]) / 10;
    block1_y = 0;// (sys->ylim2 - sys->ylim1) / 3 * sys->lengthUnit;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 10;
    block2_x = 0;// (sys->xlim2 - sys->xlim1) / 5 * sys->lengthUnit;
    block2_y = 0;// (sys->ylim2 - sys->ylim1) / 5 * sys->lengthUnit;

    /* Generate V0c */
    tf::Task taskPG6 = subflow.emplace([=]() {
        cout << "Start generating V0c" << endl;
        clock_t t1 = clock();
        sys->cindex.push_back(-1);    // Last index in the sparse form for each conductor in V0c, indi denotes ith conductor (starting from 1)
        sys->acu_cnno.push_back(0);    // How many v0c are in each conductor (cumulative), indi denotes the ith conductor (starting from 1)
        sys->merge_v0c(block1_x, block1_y, block2_x, block2_y, *v0cnum, *leng_v0c, *v0canum, *leng_v0ca, mapc, bdl, bdu);

#ifdef PRINT_VERBOSE_TIMING
        cout << " Time to generate V0c is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif
        cout << " Length of V0c is " << *leng_v0c << " number of non-zeros in V0c is " << *v0cnum << endl;
        cout << " Length of V0ca is " << *leng_v0ca << " number of non-zeros in V0ca is " << *v0canum << endl;
        cout << "V0c is generated!" << endl;
    });
    taskPG6.name("Generate V_{0c}");
    taskPG0.precede(taskPG6);

    /* Generate Ac = V0ca^T*D_sig*V0c = D_sig,0 */
    tf::Task taskPG7 = subflow.emplace([=]() {
        clock_t t1 = clock();
        sys->generateAc(mapc, *v0cnum, *v0canum, *leng_v0c, *leng_Ac);

#ifdef PRINT_VERBOSE_TIMING
        cout << "Time to generate Ac is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
        cout << "Number of non-zeros in Ac is " << *leng_Ac << endl;
#endif

        sys->v0cvalo = (double*)malloc(*v0cnum * sizeof(double));
        for (myint indi = 0; indi < *v0cnum; indi++) {
            sys->v0cvalo[indi] = sys->v0cval[indi];    // v0cvalo is the v0c values without D_sig
        }
        for (myint indi = 0; indi < *v0cnum; indi++) {
            if (sys->markEdge[sys->v0cRowId[indi]] != 0) {
                sys->v0cval[indi] *= sqrt(SIGMA);       // Compute the sparse form of D_sig*V0c
            }
        }
        sys->v0cavalo = (double*)malloc(*v0canum * sizeof(double));
        for (myint indi = 0; indi < *v0canum; indi++) {
            sys->v0cavalo[indi] = sys->v0caval[indi];
        }
        for (myint indi = 0; indi < *v0canum; indi++) {
            if (sys->markEdge[sys->v0cRowId[indi]] != 0) {
                sys->v0caval[indi] *= sqrt(SIGMA);
            }
        }

        sys->v0cColIdo = (myint*)malloc(*v0cnum * sizeof(myint));
        for (myint indi = 0; indi < *v0cnum; indi++) {
            sys->v0cColIdo[indi] = sys->v0cColId[indi];
        }
        free(sys->v0cColId); sys->v0cColId = (myint*)malloc((*leng_v0c + 1) * sizeof(myint));
        int status = COO2CSR_malloc(sys->v0cColIdo, sys->v0cRowId, sys->v0cval, *v0cnum, *leng_v0c, sys->v0cColId);
        if (status != 0) {
            cerr << "Could not make V0cColIdo into a sparse matrix (status = " << status << ")" << endl;
            exit(1);
        }
        free(sys->v0cColIdo); sys->v0cColIdo = NULL;
        free(sys->v0cval); sys->v0cval = NULL;

        /*sys->v0caColIdo = (myint*)malloc(v0canum * sizeof(myint));
        for (indi = 0; indi < v0canum; indi++) {
            sys->v0caColIdo[indi] = sys->v0caColId[indi];
        }
        free(sys->v0caColId); sys->v0caColId = (myint*)malloc((leng_v0ca + 1)*sizeof(myint));
        status = COO2CSR_malloc(sys->v0caColIdo, sys->v0caRowId, sys->v0caval, v0canum, leng_v0ca, sys->v0caColId);
        if (status != 0) {
            cerr << "Could not make V0caColIdo into a sparse matrix (status = " << status << ")" << endl;
            exit(1);
        }
        free(sys->v0caColIdo); sys->v0caColIdo = NULL;*/
        free(sys->v0caval); sys->v0caval = NULL;


        //status = mklMatrixMulti(sys, leng_Ac, sys->v0caRowId, sys->v0caColId, sys->v0caval, sys->N_edge, leng_v0c, sys->v0cRowId, sys->v0cColId, sys->v0cval, 2);

        /*for (indi = 0; indi < sys->numCdt + 1; indi++) {
            cout << sys->acu_cnno[indi] << " ";
        }
        cout << endl;
        for (indi = 0; indi < sys->numCdt + 1; indi++) {
            cout << sys->cindex[indi] << " ";
        }
        cout << endl;*/

        /* Compute the matrix V0c'*D_sig*V0c */
        /*status = matrixMulti(sys->v0caColId, sys->v0caRowId, sys->v0caval, sys->v0cRowId, sys->v0cColId, sys->v0cval, sys->AcRowId, sys->AcColId, sys->Acval);
        if (status != 0) {
            exit(1);
        }*/
    });
    taskPG7.name("Generate A_{c}");
    taskPG6.precede(taskPG7);

    /* V0ca^T's csr form handle for MKL */
    sparse_matrix_t *V0cat = (sparse_matrix_t*)calloc(1, sizeof(sparse_matrix_t));
    tf::Task taskPG8 = subflow.emplace([=]() {
        sparse_status_t s = mkl_sparse_d_create_csr(V0cat, SPARSE_INDEX_BASE_ZERO, *leng_v0c, sys->N_edge, &sys->v0cColId[0], &sys->v0cColId[1], sys->v0cRowId, sys->v0cavalo);
        if (s != 0) {
            cerr << "Could not create sparse matrix V0cat (status = " << s << ")" << endl;
        }
    });
    taskPG8.name("V_{0ca}^T to CSR Form");
    taskPG7.precede(taskPG8);

    /* V0c^T's csr form handle for MKL */
    sparse_matrix_t *V0ct = (sparse_matrix_t*)calloc(1, sizeof(sparse_matrix_t));
    tf::Task taskPG9 = subflow.emplace([=]() {
        sparse_status_t s = mkl_sparse_d_create_csr(V0ct, SPARSE_INDEX_BASE_ZERO, *leng_v0c, sys->N_edge, &sys->v0cColId[0], &sys->v0cColId[1], sys->v0cRowId, sys->v0cvalo);
        if (s != 0) {
            cerr << "Could not create sparse matrix V0ct (status = " << s << ")" << endl;
        }
    });
    taskPG9.name("V_{0c}^T to CSR Form");
    taskPG7.precede(taskPG9);


    /* Prepare Z-parameters storage  */
    tf::Task taskPG10 = subflow.emplace([=]() {
        //sys->Y = (complex<double>*)calloc(sys->numPorts * sys->numPorts * sys->nfreq, sizeof(complex<double>));
        //sys->x.assign(sys->numPorts * sys->numPorts, complex<double>(0., 0.)); // Use complex double constructor to assign initial output matrix for single-frequency solve for V0 solution
        sys->x.assign(sys->numPorts * sys->numPorts * sys->nfreq, complex<double>(0., 0.)); // Use complex double constructor to assign initial output matrix for single-frequency solve
    });
    taskPG10.name("Prepare Z-param storage");
    taskPG0.precede(taskPG10);


    /* Prepare for port solving */
    int *argc;
    char ***argv;
    lapack_int *iter;
    int *xcol = (int*)calloc(1, sizeof(int));
    //char matdescra[6];
    //matdescra[0] = 'G'; matdescra[3] = 'C';    // general matrix multi, 0-based indexing

    /* HYPRE solves for each port are messy */
    auto[PG11Start, PG11End] = subflow.parallel_for(0, sys->numPorts, 1, [=](int sourcePort) {
        cout << " Port direction for port" << sourcePort + 1 << " is " << sys->portCoor[sourcePort].portDirection[0] << endl;
#ifdef GENERATE_V0_SOLUTION
        clock_t t1 = clock();
        struct matrix_descr descr;
        double *Jport = (double*)calloc(sys->N_edge, sizeof(double));
        for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
            for (int indEdge = 0; indEdge < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
                /* Set current density for all edges within sides in port to prepare solver */
                //cout << " port #" << sourcePort + 1 << ", side #" << sourcePortSide + 1 << ", edge #" << sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge] << ": J = " << sys->portCoor[sourcePort].portDirection[sourcePortSide] << " A/m^2" << endl;
                Jport[sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]] = sys->portCoor[sourcePort].portDirection[sourcePortSide];
            }
        }

        /* MKL Level 2 (matrix-vector or "mv") operations: y <- alpha * mat_A * vec_x + beta * vec_y
        y = V0daJ: Using sparse storage, without matrix transposition, on double-precision numbers */
        double *v0daJ = (double*)calloc(*leng_v0d1, sizeof(double));
        double alpha = 1.0;
        double beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        sparse_status_t s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, *V0dat, descr, Jport, beta, v0daJ); // +V0da^T*J
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply V0da^T*J (status = " << s << ")" << endl;
        }
        for (myint indi = 0; indi < *leng_v0d1; indi++) {
            v0daJ[indi] *= -1.0; // -V0da^T*J
        }

        /* Solve V0d system */
        double *y0d = (double*)calloc(*leng_v0d1, sizeof(double));
        t1 = clock();
        //cout << " About to call hypreSolve() for y0d on port" << sourcePort + 1 << endl;
        int status = hypreSolve(sys, sys->AdRowId, sys->AdColId, sys->Adval, *leng_Ad, v0daJ, *leng_v0d1, y0d); // Numerically find D_eps,0^-1*(-V0da^T*J)
        //cout << " Finished hypreSolve() for y0d on port" << sourcePort + 1 << endl;

#ifndef SKIP_PARDISO
        t1 = clock();
        status = solveV0dSystem(sys, v0daJ, y0d, *leng_v0d1);
        cout << " Pardiso solve time " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
#endif

        for (myint indi = 0; indi < *leng_v0d1; indi++) {
            y0d[indi] /= (2 * M_PI * sys->freqStart * sys->freqUnit);    // y0d is imaginary; 1/omega * D_eps,0\(-V0da^T*J)
        }
        free(v0daJ); v0daJ = NULL;
        /* End of V0d solving */

        /* Get first vector term with V0d*y0d */
        /* MKL Level 2 (matrix-vector or "mv") operations: y <- alpha * mat_A * vec_x + beta * vec_y
        y = ydt: Using sparse storage, with matrix transposition, on double-precision numbers */
        double *ydt = (double*)calloc(sys->N_edge, sizeof(double));
        double *yd1 = (double*)malloc(sys->N_edge * sizeof(double));
        alpha = 1.0;
        beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, *V0dt, descr, y0d, beta, ydt);    // Scaled first vector term: 1/omega * V0d*(D_eps,0\(-V0da^T*J))
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply ydt (status = " << s << ")" << endl;
        }

        //lapack_complex_double *u0 = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * 2, sizeof(lapack_complex_double));
        //lapack_complex_double *u0a = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * 2, sizeof(lapack_complex_double));
        //double nn = 0;
        //double nna = 0;
        //for (indi = 0; indi < sys->N_edge; indi++) {
        //    nn += ydt[indi] * ydt[indi];
        //    nna += ydat[indi] * ydat[indi];
        //}
        //nn = sqrt(nn);
        //nna = sqrt(nna);
        //for (indi = sys->N_edge_s; indi < sys->N_edge - sys->N_edge_s; indi++) {
        //    u0[indi - sys->N_edge_s].real = ydt[indi] / nn;    // u0d is one vector in V0
        //    u0a[indi - sys->N_edge_s].real = ydat[indi] / nna;    // u0da
        //}
        for (myint indi = 0; indi < sys->N_edge; indi++) {
            yd1[indi] = ydt[indi]; // Save for later; yd1 = 1/omega * V0d*(D_eps,0\(-V0da^T*J))
            ydt[indi] *= -1.0 * (2 * M_PI * sys->freqStart * sys->freqUnit) * sys->stackEpsn[(indi + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0; // yd^T = D_eps*V0d*(D_eps,0\(+V0da^T*J))
        }
        free(y0d); y0d = NULL;


        /* MKL Level 2 (matrix-vector or "mv") operations: y <- alpha * mat_A * vec_x + beta * vec_y
        y = v0caJ: Using sparse storage, without matrix transposition, on double-precision numbers */
        double *v0caJ = (double*)calloc(*leng_v0c, sizeof(double));
        alpha = 1.0;
        beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, *V0cat, descr, Jport, beta, v0caJ); // +V0ca^T*J
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply V0ca^T*J (status = " << s << ")" << endl;
        }
        for (myint indi = 0; indi < *leng_v0c; indi++) {
            v0caJ[indi] *= -1.0; // -V0ca^T*J
        }

        /* Compute conductor region right hand side */
        double *crhs = (double*)calloc(*leng_v0c, sizeof(double));
        alpha = 1.0;
        beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, *V0cat, descr, ydt, beta, crhs); // crhs is V0ca^T*D_eps*V0d*(D_eps,0\(V0da^T*J))
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply crhs (status = " << s << ")" << endl;
        }
        free(ydt); ydt = NULL;

        double v0caJn = 0.0, crhsn = 0.0;
        for (myint indi = 0; indi < *leng_v0c; indi++) {
            v0caJ[indi] += crhs[indi]; // V0ca^T*D_eps*V0d*(D_eps,0\(-V0da^T*J)) - V0ca^T*J = V0ca^T*[D_eps*V0d*(D_eps,0\(V0da^T*J)) - I]*J
            v0caJn += v0caJ[indi] * v0caJ[indi]; // Norm for part of third term (resistive component due to conductors)
            crhsn += crhs[indi] * crhs[indi]; // Norm for conductor right hand side
        }
        v0caJn = sqrt(v0caJn);
        crhsn = sqrt(crhsn);
        free(crhs); crhs = NULL;

        /* Solve V0c system block-by-block */
        double *y0c = (double*)calloc(*leng_v0c, sizeof(double));
        //cout << " Time between the first and the second HYPRE solves is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;
        for (myint indi = 1; indi <= 1; indi++) {
#ifndef SKIP_PARDISO
            //pardisoSolve_c(sys, &v0caJ[sys->acu_cnno[indi - 1]], &y0c[sys->acu_cnno[indi - 1]], sys->acu_cnno[indi - 1], sys->acu_cnno[indi] - 1, sys->cindex[indi - 1] + 1, sys->cindex[indi]);
#endif

            //double *a = (double*)malloc((sys->cindex[indi] - sys->cindex[indi - 1]) * sizeof(double));
            //int *ia = (int*)malloc((sys->cindex[indi] - sys->cindex[indi - 1]) * sizeof(int));
            //int *ja = (int*)malloc((sys->cindex[indi] - sys->cindex[indi - 1]) * sizeof(int));
            //double *crhss = (double*)malloc((sys->acu_cnno[indi] - sys->acu_cnno[indi - 1]) * sizeof(double));
            //double *y0cs = (double*)calloc((sys->acu_cnno[indi] - sys->acu_cnno[indi - 1]), sizeof(double));

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
            //free(a); a = NULL;
            //free(ia); ia = NULL;
            //free(ja); ja = NULL;
            //free(crhss); crhss = NULL;
            //free(y0cs); y0cs = NULL;

            //cout << " About to call hypreSolve() for y0c on port" << sourcePort + 1 << endl;
            status = hypreSolve(sys, sys->AcRowId, sys->AcColId, sys->Acval, *leng_Ac, v0caJ, *leng_v0c, y0c); // Numerically find D_sig,0^-1*(V0ca^T*[D_eps*V0d*(D_eps,0\(V0da^T*J)) - I]*J)
            //cout << " Finished hypreSolve() for y0c on port" << sourcePort + 1 << endl;
        }
        free(v0caJ); v0caJ = NULL;
        /* End of V0c solving */


        /* Get third vector term with V0c*y0c */
        double *yc = (double*)calloc(sys->N_edge, sizeof(double));
        double *yccp = (double*)malloc(sys->N_edge * sizeof(double));
        alpha = 1.0;
        beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, *V0ct, descr, y0c, beta, yc); // Complete third vector term: V0c*(D_sig,0\(V0ca^T*[D_eps*V0d*(D_eps,0\(V0da^T*J)) - I]*J))
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply yc (status = " << s << ")" << endl;
        }
        free(y0c); y0c = NULL;
        for (myint indi = 0; indi < sys->N_edge; indi++) {
            yccp[indi] = -yc[indi] * sys->stackEpsn[(indi + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0; // Modified "copy" of yc includes flipped bracket term: D_eps*V0c*(D_sig,0\(V0ca^T*[I - D_eps*V0d*(D_eps,0\(V0da^T*J))]*J))
        }

        /* Compute second dielectric region right hand side */
        double *dRhs2 = (double*)calloc(*leng_v0d1, sizeof(double));
        alpha = 1.0;
        beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, *V0dat, descr, yccp, beta, dRhs2); // V0da^T*D_eps*V0c*(D_sig,0\(V0ca^T*[I - D_eps*V0d*(D_eps,0\(V0da^T*J))]*J))
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply dRhs2 (status = " << s << ")" << endl;
        }
        free(yccp); yccp = NULL;

        /* Solve y0d2 system */
        double *y0d2 = (double*)calloc(*leng_v0d1, sizeof(double));
        //cout << " About to call hypreSolve() for y0d2 on port" << sourcePort + 1 << endl;
        status = hypreSolve(sys, sys->AdRowId, sys->AdColId, sys->Adval, *leng_Ad, dRhs2, *leng_v0d1, y0d2); // Numerically find D_eps,0^-1*V0da^T*D_eps*V0c*(D_sig,0\(V0ca^T*[I - D_eps*V0d*(D_eps,0\(V0da^T*J))]*J))
        //cout << " Finished hypreSolve() for y0d2 on port" << sourcePort + 1 << endl;
        free(dRhs2); dRhs2 = NULL;
        /* End of y0d2 solving*/

        /* Get second vector term with V0d*yd2 */
        double *yd2 = (double*)calloc(sys->N_edge, sizeof(double));
        alpha = 1.0;
        beta = 0.0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, *V0dt, descr, y0d2, beta, yd2); // Complete second vector term: V0d*(D_eps,0^-1\V0da^T*D_eps*V0c*(D_sig,0\(V0ca^T*[I - D_eps*V0d*(D_eps,0\(V0da^T*J))]*J)))
        if (s != 0) {
            cerr << "Could not sparse matrix-vector multiply yd2 (status = " << s << ")" << endl;
        }
        free(y0d2); y0d2 = NULL;

        //nn = 0;
        //nna = 0;
        //for (indi = 0; indi < sys->N_edge; indi++) {
        //    nn += (yd2[indi] + yc[indi]) * (yd2[indi] + yc[indi]);
        //    nna += (yd2a[indi] + yca[indi]) * (yd2a[indi] + yca[indi]);
        //}
        //nn = sqrt(nn);
        //nna = sqrt(nna);
        //for (indi = sys->N_edge_s; indi < sys->N_edge - sys->N_edge_s; indi++) {
        //    u0[sys->N_edge - 2 * sys->N_edge_s + indi - sys->N_edge_s].real = (yd2[indi] + yc[indi]) / nn;    // u0c is the other vector in u0
        //    u0a[sys->N_edge - 2 * sys->N_edge_s + indi - sys->N_edge_s].real = (yd2a[indi] + yca[indi]) / nna;    // u0ca
        //}
        //free(yca); yca = NULL;
        //free(yd2a); yd2a = NULL;


        /* Final complex vector solution for nullspace layout response of this port */
        complex<double> *ynullsol = (complex<double>*)malloc(sys->N_edge * sizeof(complex<double>));
        for (myint indi = 0; indi < sys->N_edge; indi++) {
            ynullsol[indi] = yd2[indi] - (1i)*(yd1[indi]); // yd = V0d*(D_eps,0^-1\V0da^T*D_eps*V0c*(D_sig,0\(V0ca^T*[I - D_eps*V0d*(D_eps,0\(V0da^T*J))]*J))) - 1j/omega * V0d*(D_eps,0\(-V0da^T*J)) = 1st term + 2nd term
            ynullsol[indi] += yc[indi]; // yd = 1st term + 2nd term + 3rd term
        }

        free(yd2); yd2 = NULL;
        free(yd1); yd1 = NULL;
        free(yc); yc = NULL;

        /* Construct Z-parameters for V0 part */
        sys->Construct_Z_V0(ynullsol, sourcePort);
#endif

        /* Calculate the Vh part */
#ifndef SKIP_VH
        // Find the Vh eigenmodes
        int step = 1000;
        //sys->find_Vh(sourcePort, step, bdl, bdu);
        status = find_Vh(sys, u0, u0a, sourcePort);
        double freq; // Frequency

        for (indi = 0; indi < sys->nfreq; indi++) {
            // this point's frequency

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

            //status = V0_reference(sys, sourcePort, freq);

            // Vh = Vh - u0*((u0a'*A*u0)\(u0a'*A*Vh))
            lapack_complex_double* V_re2 = (lapack_complex_double*)malloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * sys->leng_Vh * sizeof(lapack_complex_double));
            for (myint inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {    // A*Vh
                for (myint inde2 = 0; inde2 < sys->leng_Vh; inde2++) {
                    if (sys->markEdge[inde + bdl * sys->N_edge_s] != 0) {    // if this edge is inside the conductor
                        V_re2[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real = sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
                            - 2 * M_PI * freq * sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag * SIGMA;
                        V_re2[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag = sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real * (freq * 2 * M_PI) * SIGMA
                            - sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
                    }
                    else {
                        V_re2[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real = sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
                        V_re2[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag = -sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
                    }
                }
            }

            lapack_complex_double *y_re = (lapack_complex_double*)calloc(2 * sys->leng_Vh, sizeof(lapack_complex_double));    // u0a'*A*Vh
            status = matrix_multi('T', u0a, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), 2, V_re2, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), sys->leng_Vh, y_re);
            free(V_re2); V_re2 = NULL;

            lapack_complex_double *tmp3 = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * 2, sizeof(lapack_complex_double));
            for (myint inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {    // A*u0
                for (myint inde2 = 0; inde2 < 2; inde2++) {
                    if (sys->markEdge[inde + bdl * sys->N_edge_s] != 0) {
                        tmp3[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real = u0[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
                            - 2 * M_PI * freq * u0[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag * SIGMA;
                        tmp3[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag = u0[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real * (freq * 2 * M_PI) * SIGMA
                            - u0[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
                    }
                    else {
                        tmp3[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real = u0[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real * (-freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
                        tmp3[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag = -u0[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag * (freq * freq * 4 * pow(M_PI, 2)) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0;
                    }
                }
            }

            lapack_complex_double *tmp4 = (lapack_complex_double*)calloc(2 * 2, sizeof(lapack_complex_double));    // u0a'*A*u0
            status = matrix_multi('T', u0a, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), 2, tmp3, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), 2, tmp4);    // u0a'*A*u0
            ipiv = (lapack_int*)malloc(2 * sizeof(lapack_int));
            lapack_complex_double *y_new = (lapack_complex_double*)calloc(2 * sys->leng_Vh, sizeof(lapack_complex_double));
            int *info = LAPACKE_zcgesv(LAPACK_COL_MAJOR, 2, sys->leng_Vh, tmp4, 2, ipiv, y_re, 2, y_new, 2, iter);    // (u0a'*A*u0)\(u0a'*A*V_re)
            free(ipiv); ipiv = NULL;
            free(y_re); y_re = NULL;
            free(tmp3); tmp3 = NULL;
            free(tmp4); tmp4 = NULL;

            lapack_complex_double *m_new = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * sys->leng_Vh, sizeof(lapack_complex_double));
            status = matrix_multi('N', u0, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), 2, y_new, 2, sys->leng_Vh, m_new);    // u0*((u0a'*A*u0)\(u0a'*A*V_re))
            free(y_new); y_new = NULL;
            lapack_complex_double* Vh = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * (sys->leng_Vh), sizeof(lapack_complex_double));
            for (myint inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {
                for (myint inde2 = 0; inde2 < sys->leng_Vh; inde2++) {
                    Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real = sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real - m_new[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real;
                    Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag = sys->Vh[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag - m_new[inde2 * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag;
                }
            }
            free(m_new); m_new = NULL;

            // Vh'*(A+C)*Vh
            int inde;
            myint start;
            lapack_complex_double *tmp = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * sys->leng_Vh, sizeof(lapack_complex_double));
            for (myint indj = 0; indj < sys->leng_Vh; indj++) {    // calculate (A+C)*V_re1
                inde = 0;
                while (inde < sys->leng_S) {
                    start = sys->SRowId[inde];
                    while (inde < sys->leng_S && sys->SRowId[inde] == start) {
                        tmp[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + sys->SRowId[inde]].real += sys->Sval[inde] * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + sys->SColId[inde]].real;
                        tmp[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + sys->SRowId[inde]].imag += sys->Sval[inde] * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + sys->SColId[inde]].imag;
                        inde++;
                    }
                }

                for (inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {
                    if (sys->markEdge[inde + sys->N_edge_s] != 0) {
                        tmp[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real += -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real
                            - freq * 2 * M_PI * SIGMA * sys->sig[inde + bdl * sys->N_edge_s] * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag;
                        tmp[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag += -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag
                            + freq * 2 * M_PI * SIGMA * sys->sig[inde + bdl * sys->N_edge_s] * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real;
                    }
                    else {
                        tmp[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real += -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].real;
                        tmp[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag += -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + bdl * sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * Vh[indj * (sys->N_edge - (bdl + bdu) * sys->N_edge_s) + inde].imag;
                    }
                }
            }


            lapack_complex_double *m_h = (lapack_complex_double*)calloc(sys->leng_Vh * sys->leng_Vh, sizeof(lapack_complex_double));
            status = matrix_multi('T', Vh, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), sys->leng_Vh, tmp, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), sys->leng_Vh, m_h);    // V_re1'*(A+C)*V_re1
            cout << "Vh left matrix is generated!" << endl;

            lapack_complex_double *rhs_h = (lapack_complex_double*)calloc(sys->leng_Vh * 1, sizeof(lapack_complex_double));
            lapack_complex_double *J = (lapack_complex_double*)calloc(sys->N_edge - (bdl + bdu) * sys->N_edge_s, sizeof(lapack_complex_double));
            for (inde = bdl * sys->N_edge_s; inde < sys->N_edge - bdu * sys->N_edge_s; inde++) {
                J[inde - bdl * sys->N_edge_s].imag = -Jport[inde] * freq * 2 * M_PI;

            }
            status = matrix_multi('T', Vh, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), sys->leng_Vh, J, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), 1, rhs_h);    // -1i*omega*V_re1'*J
            free(J); J = NULL;
            cout << "Vh right hand side is generated!" << endl;

            /* V_re1'*A*u */
            //free(tmp);
            //tmp = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s), sizeof(lapack_complex_double));
            //for (inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {
            //    if (sys->markEdge[inde + sys->N_edge_s] != 0) {
            //        tmp[inde].real = -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * yd[inde + sys->N_edge_s].real() - freq * 2 * M_PI * SIGMA * yd[inde + sys->N_edge_s].imag() * sys->freqStart * sys->freqUnit / freq;
            //        tmp[inde].imag = freq * 2 * M_PI * SIGMA * yd[inde + sys->N_edge_s].real() - pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * yd[inde + sys->N_edge_s].imag() * sys->freqStart * sys->freqUnit / freq;
            //    }
            //    else {
            //        tmp[inde].real = -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * yd[inde + sys->N_edge_s].real();
            //        tmp[inde].imag = -pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inde + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * yd[inde + sys->N_edge_s].imag() * sys->freqStart * sys->freqUnit / freq;
            //    }
            //}
            //lapack_complex_double *rhs_h0 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
            //status = matrix_multi('T', Vh, sys->N_edge - (bdl + bdu) * sys->N_edge_s, sys->leng_Vh, tmp, sys->N_edge - (bdl + bdu) * sys->N_edge_s, 1, rhs_h0);    // V_re1'*A*u
            //for (inde = 0; inde < sys->leng_Vh; inde++) {
            //    rhs_h[inde].real = rhs_h[inde].real - rhs_h0[inde].real;
            //    rhs_h[inde].imag = rhs_h[inde].imag - rhs_h0[inde].imag;
            //}
            //free(rhs_h0); rhs_h0 = NULL;
            //free(tmp); tmp = NULL;


            lapack_int *ipiv = (lapack_int*)malloc(sys->leng_Vh * sizeof(lapack_int));
            lapack_int info1 = LAPACKE_zgesv(LAPACK_COL_MAJOR, sys->leng_Vh, 1, m_h, sys->leng_Vh, ipiv, rhs_h, sys->leng_Vh);// , y_h, sys->leng_Vh, iter);    // yh is generated
            free(ipiv); ipiv = NULL;
            free(m_h); m_h = NULL;

            lapack_complex_double *y_h = (lapack_complex_double*)calloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s), sizeof(lapack_complex_double));
            status = matrix_multi('N', Vh, (sys->N_edge - (bdl + bdu) * sys->N_edge_s), sys->leng_Vh, rhs_h, sys->leng_Vh, 1, y_h);

            complex<double> *final_x = (complex<double>*)malloc((sys->N_edge - (bdl + bdu) * sys->N_edge_s) * sizeof(complex<double>));
            for (inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {
                final_x[inde] = yd[inde + bdl * sys->N_edge_s].real() + y_h[inde].real + 1indi * (yd[inde + bdu * sys->N_edge_s].imag() *sys->freqStart * sys->freqUnit / freq + y_h[inde].imag);
            }

            free(y_h); y_h = NULL;
            free(rhs_h); rhs_h = NULL;
            cout << "final_x is generated!" << endl;
            /*xr = (complex<double>*)calloc(sys->N_edge - (bdl + bdu) * sys->N_edge_s, sizeof(complex<double>));
            status = reference1(sys, freq, sourcePort, sys->SRowId, sys->SColId, sys->Sval, xr);
            cout << "xr is generated!\n";*/
            // Construct Z parameters
            sys->Construct_Z_V0_Vh(final_x, indi, sourcePort);

            /* Compare with the solution from (-w^2*D_eps+iw*D_sig+S)xr=-iwJ */
            /*double err = 0, total_norm = 0, err0 = 0;
            for (inde = 0; inde < sys->N_edge - (bdl + bdu) * sys->N_edge_s; inde++) {
                err += sqrt((xr[inde].real() - final_x[inde].real()) * (xr[inde].real() - final_x[inde].real()) + (xr[inde].imag() - final_x[inde].imag()) * (xr[inde].imag() - final_x[inde].imag()));
                err0 += sqrt((xr[inde].real() - yd[inde + sys->N_edge_s].real()) * (xr[inde].real() - yd[inde + sys->N_edge_s].real()) + (xr[inde].imag() - yd[inde + sys->N_edge_s].imag() *sys->freqStart * sys->freqUnit / freq) * (xr[inde].imag() - yd[inde + sys->N_edge_s].imag() *sys->freqStart * sys->freqUnit / freq));
                total_norm += sqrt((xr[inde].real()) * (xr[inde].real()) + (xr[inde].imag()) * (xr[inde].imag()));
            }


            cout << "Freq " << freq << " the total error is " << err / total_norm << endl;
            cout << "Freq " << freq << " the y0 total error is " << err0 / total_norm << endl;*/


            free(final_x); final_x = NULL;
            free(Vh); Vh = NULL;
        }
#endif

#ifdef GENERATE_V0_SOLUTION
        free(ynullsol); ynullsol = NULL;
        free(Jport); Jport = NULL;
#endif

        // Solve system for x in (-omega^2 * D_eps + j * omega * D_sigma + S) * x = -j * omega * J
        /*free(u0); u0 = NULL;*/
        /*free(u0a); u0a = NULL;*/
        (*xcol)++;
    });
    PG11Start.work([=]() {
        cout << endl << "Begin to solve for network parameters!" << endl;
        MPI_Init(argc, argv);
        free(mapd);
        free(mapc);
        free(sys->markNode); sys->markNode = NULL;
        for (myint indi = 0; indi < sys->numCdt; indi++) {
            free(sys->conductor[indi].node); sys->conductor[indi].node = NULL;
        }
        free(sys->conductor); sys->conductor = NULL;
    });
    PG11Start.name("Synchronize for Parallel HYPRE Solves");
    PG11End.work([=]() {
        cout << "Finished solve for network parameters!" << endl;
        MPI_Finalize();
        mkl_sparse_destroy(*V0dt);
        mkl_sparse_destroy(*V0dat);
        mkl_sparse_destroy(*V0ct);
        mkl_sparse_destroy(*V0cat);

        free(leng_v0d1);
        free(v0d1num);
        free(leng_v0d1a);
        free(v0d1anum);
        free(leng_v0d2);
        free(leng_v0d2a);
        free(v0d2num);
        free(v0d2anum);
        free(leng_Ad);
        free(V0dt);
        free(V0dat);
        free(leng_v0c);
        free(v0cnum);
        free(leng_v0ca);
        free(v0canum);
        free(leng_Ac);
        free(xcol);
        free(V0cat);
        free(V0ct);
    });
    PG11End.name("Finalize Parallel HYPRE Solves");
    taskPG3.precede(PG11Start);
    taskPG4.precede(PG11Start);
    taskPG8.precede(PG11Start);
    taskPG9.precede(PG11Start);
    taskPG10.precede(PG11Start);

    tf::Task taskPG12 = subflow.emplace([=]() {
#ifndef SKIP_STIFF_REFERENCE
        /*  Generate the reference results and S parameters in .citi file for different frequencies with multiple right hand side */
        double freq; // Frequency

        for (myint indi = 0; indi < sys->nfreq; indi++) {

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
            cout << "Z-parameters matrix at frequency " << freq << " Hz is shown below:" << endl;
            status = reference(sys, indi, sys->SRowId, sys->SColId, sys->Sval, bdl, bdu);

            /*for (int indj = 0; indj < sys->numPorts; indj++) {
            for (int indk = 0; indk < sys->numPorts; indk++) {
            cout << sys->x[indi * (sys->numPorts * sys->numPorts) + indj * sys->numPorts + indk] << " ";
            }
            cout << endl;
            }*/
        }
#endif
    });
    taskPG12.name("Find Reference Results");
    PG11End.precede(taskPG12);

    /* Report the Z-parameters and prepare them for export */
    tf::Task taskPG13 = subflow.emplace([=]() {
#ifdef PRINT_V0_Z_PARAM
        sys->print_z_V0();
#endif

#ifdef PRINT_V0_Vh_Z_PARAM
        sys->print_z_V0_Vh();
#endif

        sys->cindex.clear();
        sys->acu_cnno.clear();

        free(sys->AdColId); sys->AdColId = NULL;
        free(sys->Adval); sys->Adval = NULL;
        free(sys->AdRowId); sys->AdRowId = NULL;
        free(sys->v0d1RowId); sys->v0d1RowId = NULL;
        //free(sys->v0d1ColId); sys->v0d1ColId = NULL;
        //free(sys->v0d1ColIdo); sys->v0d1ColIdo = NULL;
        //free(sys->v0d1val); sys->v0d1val = NULL;
        free(sys->v0d1valo); sys->v0d1valo = NULL;
        //free(sys->v0d1aRowId); sys->v0d1aRowId = NULL;
        //free(sys->v0d1aColId); sys->v0d1aColId = NULL;
        //free(sys->v0d1aColIdo); sys->v0d1aColIdo = NULL;
        //free(sys->v0d1aval); sys->v0d1aval = NULL;
        free(sys->v0d1avalo); sys->v0d1avalo = NULL;
        free(sys->v0cRowId);  sys->v0cRowId = NULL;
        free(sys->v0cColId);  sys->v0cColId = NULL;
        //free(sys->v0cColIdo); sys->v0cColIdo = NULL;
        //free(sys->v0cval); sys->v0cval = NULL;
        free(sys->v0cvalo); sys->v0cvalo = NULL;
        //free(sys->v0caRowId); sys->v0caRowId = NULL;
        //free(sys->v0caColId); sys->v0caColId = NULL;
        //free(sys->v0caColIdo); sys->v0caColIdo = NULL;
        //free(sys->v0caval); sys->v0caval = NULL;
        free(sys->v0cavalo); sys->v0cavalo = NULL;
        free(sys->AcRowId); sys->AcRowId = NULL;
        free(sys->AcColId); sys->AcColId = NULL;
        free(sys->Acval); sys->Acval = NULL;
        free(sys->xn); sys->xn = NULL;
        free(sys->yn); sys->yn = NULL;
        free(sys->zn); sys->zn = NULL;
        sys->stackEpsn.clear();
        //cout << "Here freeing pointers in sys" << endl;
        cout << "paraGenerator dynamic execution time is " << (clock() - t0) * 1.0 / CLOCKS_PER_SEC << " s" << endl << endl;
    });
    taskPG13.name("Report / Export Z-param");
    taskPG12.precede(taskPG13);
    return 0;
}

#ifndef SKIP_PARDISO
int pardisoSolve_c(fdtdMesh *sys, double *rhs, double *solution, int nodestart, int nodeend, int indstart, int indend) {

    /* Solve Ac system block by block*/
    int leng = nodeend - nodestart + 1;
    int status;
    double *a = (double*)malloc((indend - indstart + 1) * sizeof(double));
    int *ia = (int*)malloc((indend - indstart + 1) * sizeof(int));
    int *ja = (int*)malloc((indend - indstart + 1) * sizeof(int));

    for (int indi = 0; indi < (indend - indstart + 1); indi++) {
        a[indi] = sys->Acval[indstart + indi];
        ia[indi] = sys->AcRowId[indstart + indi] - sys->AcRowId[indstart];
        ja[indi] = sys->AcColId[indstart + indi] - sys->AcColId[indstart];
    }

    /* Use pardiso to solve */
    {
        int *ia1 = (int*)malloc((leng + 1) * sizeof(int));
        /* status = COO2CSR_malloc(ia, ja, a, indstart - indend + 1, leng, ia1);
         if (status != 0)
             return status;*/
        int count = 0;
        int indi = 0;
        int indk = 0;
        int start;
        ia1[indk] = 0;
        indk++;
        while (indi < (indend - indstart + 1)) {
            start = ia[indi];
            while (indi < (indend - indstart + 1) && ia[indi] == start) {
                count++;
                indi++;
            }
            ia1[indk] = (count);
            indk++;
        }

        void *pt[64];
        int mtype;
        int iparm[64];
        double dparm[64];
        int maxfct, mnum, phase, error, solver;
        int num_process;   //number of processors
        int v0csin;
        int perm;
        int nrhs = 1;
        int msglvl = 0;    //print statistical information

        mtype = 11;    // real and not symmetric
        solver = 0;
        error = 0;
        maxfct = 1;    //maximum number of numerical factorizations
        mnum = 1;    //which factorization to use
        phase = 13;    //analysis

        pardisoinit(pt, &mtype, iparm);
        nrhs = 1;
        iparm[38] = 1;
        iparm[34] = 1;    //0-based indexing

        pardiso(pt, &maxfct, &mnum, &mtype, &phase, &leng, a, ia1, ja, &perm, &nrhs, iparm, &msglvl, rhs, solution, &error);
        free(ia1); ia1 = NULL;
    }

    /* use HYPRE to solve */
    /*{
        status = hypreSolve(sys, ia, ja, a, indend - indstart + 1, rhs, leng, solution);
    }*/

    free(a); a = NULL;
    free(ia); ia = NULL;
    free(ja); ja = NULL;

    return 0;
}
#endif

//int mklMatrixMulti(fdtdMesh *sys, int &leng_A, int *aRowId, int *aColId, double *aval, int arow, int acol, int *bRowId, int *bColId, double *bval, int mark) {
//    // ArowId, AcolId, and Aval should be in the COO format
//    sparse_status_t s0;
//    sparse_matrix_t a, a_csr;
//    sparse_index_base_t indexing1 = SPARSE_INDEX_BASE_ZERO;
//    int row, col;
//    int *cols, *cole, *rowi;
//    double *val;
//    MKL_INT *AcolId;
//    double *Aval;
//    int indk;
//
//    s0 = mkl_sparse_d_create_csc(&a, SPARSE_INDEX_BASE_ZERO, arow, acol, &aColId[0], &aColId[1], aRowId, aval);
//
//
//    sparse_matrix_t b, b_csr;
//    s0 = mkl_sparse_d_create_csc(&b, SPARSE_INDEX_BASE_ZERO, arow, acol, &bColId[0], &bColId[1], bRowId, bval);
//    s0 = mkl_sparse_convert_csr(a, SPARSE_OPERATION_NON_TRANSPOSE, &a_csr);
//    s0 = mkl_sparse_convert_csr(b, SPARSE_OPERATION_NON_TRANSPOSE, &b_csr);
//
//    sparse_matrix_t A;
//    s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, a_csr, b_csr, &A);
//
//    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
//    int ARows, ACols;
//    MKL_INT *ArowStart, *ArowEnd;
//    s0 = mkl_sparse_d_export_csr(A, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
//
//    leng_A = ArowEnd[ARows - 1];
//
//    if (mark == 1) {    // dielectric
//        sys->AdRowId = (int*)malloc(leng_A * sizeof(int));
//        sys->AdColId = (int*)malloc(leng_A * sizeof(int));
//        sys->Adval = (double*)malloc(leng_A * sizeof(double));
//        int count, num, indj;
//
//        indj = 0;
//        for (int indi = 0; indi < ARows; indi++) {
//            num = ArowEnd[indi] - ArowStart[indi];
//            count = 0;
//            while (count < num) {
//                sys->AdRowId[indj] = indi;
//                sys->AdColId[indj] = AcolId[indj];
//                sys->Adval[indj] = Aval[indj];
//                indj++;
//                count++;
//            }
//        }
//    }
//    else if (mark == 2) {    // conductor
//        sys->AcRowId = (int*)malloc(leng_A * sizeof(int));
//        sys->AcColId = (int*)malloc(leng_A * sizeof(int));
//        sys->Acval = (double*)malloc(leng_A * sizeof(double));
//        int count, num, indj;
//
//        indj = 0;
//        indk = 1;
//        sys->cindex[indk] = sys->cindex[indk - 1];
//        for (int indi = 0; indi < ARows; indi++) {
//            num = ArowEnd[indi] - ArowStart[indi];
//            count = 0;
//            while (count < num) {
//                sys->AcRowId[indj] = indi;
//                sys->AcColId[indj] = AcolId[indj];
//                sys->Acval[indj] = Aval[indj];
//                if (indk - 1 <= indi) {
//                    sys->cindex[indk]++;
//                }
//                else {
//                    indk++;
//                    sys->cindex[indk] = sys->cindex[indk - 1];
//                    sys->cindex[indk]++;
//                }
//                indj++;
//                count++;
//            }
//        }
//    }
//
//    mkl_sparse_destroy(a);
//    mkl_sparse_destroy(b);
//    mkl_sparse_destroy(A);
//
//    return 0;
//}

int matrixMul(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> bRowId, vector<int> bColId, vector<double> bval, vector<int> &cRowId, vector<int> &cColId, vector<double> &cval) {
    // The first matrix is row by row, the second matrix is column by column
    int indi = 0, indj = 0;
    int flaga = aRowId[0];
    int flagb = bColId[0];
    int starta = 0;
    double sum = 0.;
    while (indi < aRowId.size()) {
        while (indj < bColId.size() && bColId[indj] == flagb && indi < aRowId.size() && aRowId[indi] == flaga) {
            if (aColId[indi] == bRowId[indj]) {
                sum += aval[indi] * bval[indj];
                indj++;
                indi++;
            }
            else if (aColId[indi] < bRowId[indj]) {
                indi++;
            }
            else if (aColId[indi] > bRowId[indj]) {
                indj++;
            }
        }
        if (sum != 0) {
            cRowId.push_back(flaga);
            cColId.push_back(flagb);
            cval.push_back(sum);
            sum = 0;
        }
        if (indi == aRowId.size()) {
            if (indj == bColId.size())
                break;
            else {
                indi = starta;
                while (bColId[indj] == flagb) {
                    indj++;
                    if (indj == bColId.size()) {
                        while (indi < aRowId.size() && aRowId[indi] == flaga) {
                            indi++;
                        }
                        starta = indi;
                        if (indi == aRowId.size())    //run all of the datas
                            break;
                        flaga = aRowId[indi];
                        indj = 0;
                        break;
                    }
                }
                flagb = bColId[indj];
                continue;
            }
        }
        if (indj == bColId.size()) {
            while (indi < aRowId.size() && aRowId[indi] == flaga) {
                indi++;
            }
            starta = indi;
            if (indi == aRowId.size())    //run all of the datas
                break;
            flaga = aRowId[indi];
            indj = 0;
        }
        else {
            if (bColId[indj] != flagb && aRowId[indi] != flaga) {
                flagb = bColId[indj];
                indi = starta;
            }
            else if (bColId[indj] != flagb) {
                flagb = bColId[indj];
                indi = starta;
            }
            else if (aRowId[indi] != flaga) {
                indi = starta;
                while (bColId[indj] == flagb) {
                    indj++;
                    if (indj == bColId.size()) {
                        while (indi < aRowId.size() && aRowId[indi] == flaga) {
                            indi++;
                        }
                        starta = indi;
                        if (indi == aRowId.size())    //run all of the datas
                            break;
                        flaga = aRowId[indi];
                        indj = 0;
                        break;
                    }
                }
                flagb = bColId[indj];
            }
        }
    }

    return 0;
}

int COO2CSR(vector<int> &rowId, vector<int> &ColId, vector<double> &val) {
    int indi = 0;
    vector<int> rowId2;
    int count = 0, start = 0;

    rowId2.push_back(0);
    while (indi < rowId.size()) {
        start = rowId[indi];
        while (indi < rowId.size() && rowId[indi] == start) {
            count++;
            indi++;
        }
        rowId2.push_back(count);
    }

    rowId.clear();
    rowId = rowId2;
    return 0;
}

#ifndef SKIP_PARDISO
int solveV0dSystem(fdtdMesh *sys, double *dRhs, double *y0d, int leng_v0d1) {

    clock_t t1 = clock();
    /* A\b1 */
    double *d = &(sys->Adval[0]);
    int *id = &(sys->AdRowId1[0]);
    int *jd = &(sys->AdColId[0]);

    void *ptd[64];
    int mtyped;
    int iparmd[64];
    double dparmd[64];
    int maxfctd, mnumd, phased, errord, solverd;
    int num_processd;   //number of processors
    int v0csin;
    int permd;
    int nrhs = 1;
    int msglvld = 0;    //print statistical information

    mtyped = 11;    // real and not symmetric
    solverd = 0;
    errord = 0;
    maxfctd = 1;    //maximum number of numerical factorizations
    mnumd = 1;    //which factorization to use
    phased = 13;    //analysis

    pardisoinit(ptd, &mtyped, iparmd);
    nrhs = 1;
    iparmd[38] = 1;
    iparmd[34] = 1;    //0-based indexing
    pardiso(ptd, &maxfctd, &mnumd, &mtyped, &phased, &leng_v0d1, d, id, jd, &permd, &nrhs, iparmd, &msglvld, dRhs, y0d, &errord);
    cout << "Time to this point: " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << " s" << endl;

    return 0;
}
#endif

int setsideLen(int node, double sideLen, int *markLayerNode, int *markProSide, fdtdMesh *sys) {
    queue<int> q;
    int *visited;
    int indx = 0, indy = 0;
    int mark;

    if (sideLen == 0) {
        return 0;
    }

    visited = (int*)calloc(sys->N_node_s, sizeof(int));
    q.push(node);
    visited[node] = 1;
    double startx = sys->xn[node / (sys->N_cell_y + 1)];
    double starty = sys->yn[node % (sys->N_cell_y + 1)];
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
