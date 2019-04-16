//#include "stdafx.h"
#include <ctime>
#include "fdtd.h"
#include "hypreSolverh.h"


static bool comp(pair<double, int> a, pair<double, int> b)
{
    return a.first <= b.first;
};

int setHYPREMatrix(myint *ARowId, myint *AColId, double *Aval, myint leng_v0, HYPRE_IJMatrix &a, HYPRE_ParCSRMatrix &parcsr_a){
    HYPRE_Int i;
    int myid, num_procs;
    HYPRE_Int N, n;

    HYPRE_Int ilower, iupper;
    HYPRE_Int local_size, extra;

    int solver_id;
    HYPRE_Int vis, print_system;

    double h, h2;

    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;

    /* Initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    N = leng_v0;
    local_size = N / num_procs;
    extra = N - local_size * num_procs;

    ilower = local_size * myid;
    ilower += hypre_min(myid, extra);

    iupper = local_size * (myid + 1);
    iupper += hypre_min(myid + 1, extra);
    iupper = iupper - 1;

    /* How many rows do I have? */
    local_size = iupper - ilower + 1;

    /* Create the matrix.
    Note that this is a square matrix, so we indicate the row partition
    size twice (since number of rows = number of cols) */
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

    /* Choose a parallel csr format storage (see the User's Manual) */
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    /* Initialize before setting coefficients */
    HYPRE_IJMatrixInitialize(A);

    {
        HYPRE_Int nnz;
        vector<double> values;
        vector<HYPRE_Int> cols;
        int index = 0;

        for (i = ilower; i <= iupper; i++)
        {
            nnz = 0;   // Number of non-zeros on row i

            while (ARowId[index] == i){
                cols.push_back(AColId[index]);
                values.push_back(Aval[index]);
                nnz++;
                index++;
            }

            /* Set the values for row i */
            HYPRE_IJMatrixSetValues(A, 1, &nnz, &i, &cols[0], &values[0]);

            cols.clear();
            values.clear();
        }
    }

    /* Assemble after setting the coefficients */
    HYPRE_IJMatrixAssemble(A);

    /* Get the parcsr matrix object to use */
    HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);

    a = A;
    parcsr_a = parcsr_A;

    return 0;
}

int paraGenerator(fdtdMesh *sys, unordered_map<double, int> xi, unordered_map<double, int> yi, unordered_map<double, int> zi){
    myint i, j, mark, k, l, n;
    int status;
    int count;
    int xcol;
    ofstream outfile1;
    vector<int> rowId;
    vector<int> colId;
    vector<double> val;
    vector<int> temp2;
    int inx, iny, inz;

    cout << "Number of edges in one layer: " << sys->N_edge_s << endl;
    cout << "Number of edges: " << sys->N_edge << endl;
    cout << "Number of nodes in one layer: " << sys->N_node_s << endl;
    cout << "Number of nodes: " << sys->N_node << endl;

    /* Construct V0d with row id, col id and its val */
    myint leng_v0d1 = 0, v0d1num = 0;    // store the num of v0d1 vectors, which are nodes outside the conductors
    myint leng_v0d1a = 0, v0d1anum = 0;
    myint leng_Ad = 0;
    myint *map = (myint*)calloc(sys->N_node, (sizeof(myint)));
    double block1_x, block1_y, block2_x, block2_y, block3_x, block3_y, sideLen;
    unordered_map<myint, unordered_map<myint, double>> Ad1;

    block1_x = 0;    // 
    block1_y = 0;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 10;
    block2_x = 0;// (sys->xn[sys->nx - 1] - sys->xn[0]) / 50;
    block2_y = 0;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 50;
    block3_x = 0;
    block3_y = 0;
    cout << "V0d's block1_x and block1_y are " << block1_x << " " << block1_y << endl;
    cout << "V0d's block2_x and block2_y are " << block2_x << " " << block2_y << endl;
    cout << "V0d's block3_x and block3_y are " << block3_x << " " << block3_y << endl;
    sideLen = 0;    // around the conductor 10um is considered with rapid potential change


    clock_t t1 = clock();
    status = merge_v0d1(sys, block1_x, block1_y, block2_x, block2_y, block3_x, block3_y, v0d1num, leng_v0d1, v0d1anum, leng_v0d1a, map, sideLen);
    cout << "Merge V0d1 time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
    
    /*outfile1.open("map.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_node; i++){
        outfile1 << map[i] << endl;
    }
    outfile1.close();*/
    
    for (i = 0; i < v0d1anum; i++){    // the upper and lower planes are PEC
        if (map[sys->edgelink[sys->v0d1aRowId[i] * 2]] != sys->v0d1aColId[i] + 1 && map[sys->edgelink[sys->v0d1aRowId[i] * 2]] != 0){
            Ad1[sys->v0d1aColId[i]][map[sys->edgelink[sys->v0d1aRowId[i] * 2]] - 1] += sys->v0d1aval[i] * 1 / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->eps[sys->v0d1aRowId[i]];
            Ad1[sys->v0d1aColId[i]][sys->v0d1aColId[i]] += sys->v0d1aval[i] * (-1) / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->eps[sys->v0d1aRowId[i]];
        }
        else if (map[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1]] != sys->v0d1aColId[i] + 1 && map[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1]] != 0){
            Ad1[sys->v0d1aColId[i]][map[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1]] - 1] += sys->v0d1aval[i] * (-1) / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->eps[sys->v0d1aRowId[i]];
            Ad1[sys->v0d1aColId[i]][sys->v0d1aColId[i]] += sys->v0d1aval[i] * 1 / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->eps[sys->v0d1aRowId[i]];
        }
        else {//if (map[sys->edgelink[sys->v0d1aRowId[i] * 2]] == 0 || map[sys->edgelink[sys->v0d1aRowId[i] * 2] + 1] == 0){
            Ad1[sys->v0d1aColId[i]][sys->v0d1aColId[i]] += abs(sys->v0d1aval[i] * 1 / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0d1aRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->eps[sys->v0d1aRowId[i]]);
        }
    }
    for (i = 0; i < leng_v0d1; i++){
        leng_Ad += Ad1[i].size();
    }

    sys->v0d1valo = (double*)malloc(v0d1num * sizeof(double));
    for (i = 0; i < v0d1num; i++)
        sys->v0d1valo[i] = sys->v0d1val[i];
    for (i = 0; i < v0d1num; i++){    // compute sqrt(D_eps)*V0d1
        sys->v0d1val[i] = sys->v0d1val[i] * sqrt(sys->eps[sys->v0d1RowId[i]]);
    }
    sys->v0d1avalo = (double*)malloc(v0d1anum * sizeof(double));
    for (i = 0; i < v0d1anum; i++)
        sys->v0d1avalo[i] = sys->v0d1aval[i];
    for (i = 0; i < v0d1anum; i++){
        sys->v0d1aval[i] = sys->v0d1aval[i] * sqrt(sys->eps[sys->v0d1aRowId[i]]);
    }
    

    myint leng_v0d = leng_v0d1;
    sys->v0d1ColIdo = (myint*)malloc(v0d1num * sizeof(myint));
    for (i = 0; i < v0d1num; i++)
        sys->v0d1ColIdo[i] = sys->v0d1ColId[i];
    free(sys->v0d1ColId); sys->v0d1ColId = (myint*)malloc((leng_v0d1 + 1) * sizeof(myint));
    status = COO2CSR_malloc(sys->v0d1ColIdo, sys->v0d1RowId, sys->v0d1val, v0d1num, leng_v0d1, sys->v0d1ColId);
    if (status != 0)
        return status;
    sys->v0d1aColIdo = (myint*)malloc(v0d1anum * sizeof(myint));
    for (i = 0; i < v0d1anum; i++)
        sys->v0d1aColIdo[i] = sys->v0d1aColId[i];
    free(sys->v0d1aColId); sys->v0d1aColId = (myint*)malloc((leng_v0d1a + 1) * sizeof(myint));
    status = COO2CSR_malloc(sys->v0d1aColIdo, sys->v0d1aRowId, sys->v0d1aval, v0d1anum, leng_v0d1a, sys->v0d1aColId);
    if (status != 0)
        return status;
    
    cout << "Length of V0d1 is " << leng_v0d1 << endl;
    cout << "Number of NNZ in V0d1 is " << v0d1num << endl;
    
    sparse_status_t s;

    /* V0d^T's csr form handle for MKL */
    sparse_matrix_t V0dt;
    s = mkl_sparse_d_create_csr(&V0dt, SPARSE_INDEX_BASE_ZERO, leng_v0d1, sys->N_edge, &sys->v0d1ColId[0], &sys->v0d1ColId[1], sys->v0d1RowId, sys->v0d1valo);

    /* V0da^T's csr form handle for MKL */
    sparse_matrix_t V0dat;
    s = mkl_sparse_d_create_csr(&V0dat, SPARSE_INDEX_BASE_ZERO, leng_v0d1, sys->N_edge, &sys->v0d1aColId[0], &sys->v0d1aColId[1], sys->v0d1aRowId, sys->v0d1avalo);

    ///******************************************************************/
    ///* use MKL to do matrix multiplication */
    ///* Note that the result for each row, the col # is not in the increased order */
    /*leng_Ad = 0;
    clock_t t2 = clock();
    status = mklMatrixMulti(sys, leng_Ad, sys->v0d1aRowId, sys->v0d1aColId, sys->v0d1aval, sys->N_edge, leng_v0d1, sys->v0d1RowId, sys->v0d1ColId, sys->v0d1val, 1);
    cout << "Matrix mutiplication time is " << (clock() - t2) * 1.0 / CLOCKS_PER_SEC << endl;*/
    
    ///*****************************************************************/
    sys->AdRowId = (myint*)calloc(leng_Ad, sizeof(myint));
    sys->AdColId = (myint*)calloc(leng_Ad, sizeof(myint));
    sys->Adval = (double*)calloc(leng_Ad, sizeof(double));
    j = 0;
    for (i = 0; i < leng_v0d1; i++){
        vector<pair<myint, double>> v(Ad1[i].begin(), Ad1[i].end());
        sort(v.begin(), v.end());
        for (auto adi : v){
            if (abs(adi.second) > 1e-8){
                sys->AdRowId[j] = i;
                sys->AdColId[j] = adi.first;
                sys->Adval[j] = adi.second;
                j++;
            }
        }
        v.clear();
    }
    Ad1.clear();

    int *argc;
    char ***argv;
    /*  trial of first set HYPRE matrix Ad */
    HYPRE_IJMatrix ad;
    HYPRE_ParCSRMatrix parcsr_ad;
    MPI_Init(argc, argv);
    status = setHYPREMatrix(sys->AdRowId, sys->AdColId, sys->Adval, leng_v0d1, ad, parcsr_ad);
   
    /* End */
    
    //sys->AdRowId1 = (int*)malloc((leng_v0d1 + 1) * sizeof(int));
    //status = COO2CSR_malloc(sys->AdRowId, sys->AdColId, sys->Adval, leng_Ad, leng_v0d1, sys->AdRowId1);
    //if (status != 0)
    //    return status;

    cout << "The number of nonzeros in Ad is " << leng_Ad << endl;
    /*outfile1.open("Ad.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < leng_Ad; i++){
        outfile1 << sys->AdRowId[i] + 1 << " " << sys->AdColId[i] + 1 << " " << sys->Adval[i] << endl;
    }
    outfile1.close();*/

    /* Construct V0c with row id, col id and its val */
    myint leng_v0c = 0, v0cnum = 0;
    myint leng_v0ca = 0, v0canum = 0;

    int numPortCdt = 0;
    myint leng_Ac = 0;
    count = 0;


    sys->cindex = (int*)calloc((sys->numCdt + 1), (sizeof(int)));
    sys->acu_cnno = (int*)calloc((sys->numCdt + 1), sizeof(int));
    sys->cindex[0] = -1;    // the last index in the sparse form for each conductor in V0c, i denotes ith conductor (starting from 1)
    sys->acu_cnno[0] = 0;    // how many v0c are in each conductor (accumulative), i denotes the ith conductor (starting from 1)
    


    unordered_map<myint, unordered_map<myint, double>> Ac;
    count = 0;
    free(map);
    map = (myint*)calloc(sys->N_node, sizeof(myint));
    block1_x = 0;// (sys->xn[sys->nx - 1] - sys->xn[0]) / 10;
    block1_y = 0;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 10;
    block2_x = 0;// (sys->xn[sys->nx - 1] - sys->xn[0]) / 300;
    block2_y = 0;// (sys->yn[sys->ny - 1] - sys->yn[0]) / 300;
    //cout << "V0c's block1_x and block1_y are " << block1_x << " " << block1_y << endl;
    //cout << "V0c's block2_x and block2_y are " << block2_x << " " << block2_y << endl;
    status = merge_v0c(sys, block1_x, block1_y, block2_x, block2_y, v0cnum, leng_v0c, v0canum, leng_v0ca, map);
    j = 0;

    /*for (i = 0; i < sys->numCdt + 1; i++){
        cout << sys->acu_cnno[i] << " ";
    }
    cout << endl;*/

    for (i = 0; i < v0canum; i++){    // the upper and lower planes are PEC
        if (map[sys->edgelink[sys->v0caRowId[i] * 2]] != sys->v0caColId[i] + 1 && map[sys->edgelink[sys->v0caRowId[i] * 2]] != 0){
            Ac[sys->v0caColId[i]][map[sys->edgelink[sys->v0caRowId[i] * 2]] - 1] += sys->v0caval[i] * 1 / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->sig[sys->v0caRowId[i]];
            Ac[sys->v0caColId[i]][sys->v0caColId[i]] += sys->v0caval[i] * (-1) / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->sig[sys->v0caRowId[i]];
        }
        else if (map[sys->edgelink[sys->v0caRowId[i] * 2 + 1]] != sys->v0caColId[i] + 1 && map[sys->edgelink[sys->v0caRowId[i] * 2 + 1]] != 0){
            Ac[sys->v0caColId[i]][map[sys->edgelink[sys->v0caRowId[i] * 2 + 1]] - 1] += sys->v0caval[i] * (-1) / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->sig[sys->v0caRowId[i]];
            Ac[sys->v0caColId[i]][sys->v0caColId[i]] += sys->v0caval[i] * 1 / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->sig[sys->v0caRowId[i]];
        }
        else {
            Ac[sys->v0caColId[i]][sys->v0caColId[i]] += abs(sys->v0caval[i] * 1 / sqrt(pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 1], 2)
                + pow(sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->v0caRowId[i] * 2 + 1] * 3 + 2], 2)) * sys->sig[sys->v0caRowId[i]]);
        }
    }
    sys->cindex[j + 1] = v0canum - 1;
    for (i = 0; i < leng_v0c; i++){
        leng_Ac += Ac[i].size();
    }
    sys->AcRowId = (myint*)calloc(leng_Ac, sizeof(myint));
    sys->AcColId = (myint*)calloc(leng_Ac, sizeof(myint));
    sys->Acval = (double*)calloc(leng_Ac, sizeof(double));
    j = 0;
    k = 1;
    for (i = 0; i < leng_v0c; i++){
        vector<pair<myint, double>> v(Ac[i].begin(), Ac[i].end());
        sort(v.begin(), v.end());
        for (auto aci : v){
            if (abs(aci.second) > 1e5){
                sys->AcRowId[j] = i;
                sys->AcColId[j] = aci.first;
                sys->Acval[j] = aci.second;
                if (sys->AcRowId[j] >= sys->acu_cnno[k]){
                    sys->cindex[k] = j - 1;
                    k++;
                }
                j++;
            }
        }
        v.clear();
    }
    sys->cindex[k] = j - 1;
    
    
    Ac.clear();

    /*  trial of first set HYPRE matrix Ac */
    HYPRE_IJMatrix ac;
    HYPRE_ParCSRMatrix parcsr_ac;
    status = setHYPREMatrix(sys->AcRowId, sys->AcColId, sys->Acval, leng_v0c, ac, parcsr_ac);
    /* End */


    /*outfile1.open("map.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_node; i++){
        outfile1 << map[i] << endl;
    }
    outfile1.close();*/
    /*cout << "error " << endl;
    for (i = 0; i < v0cnum; i++){
        if (map[sys->edgelink[sys->v0cRowId[i] * 2]] == map[sys->edgelink[sys->v0cRowId[i] * 2 + 1]] && map[sys->edgelink[sys->v0cRowId[i] * 2 + 1]] != 0){
            cout << sys->edgelink[sys->v0cRowId[i] * 2] << " " << sys->edgelink[sys->v0cRowId[i] * 2 + 1] << endl;
        }
    }*/
    
    
    
    sys->v0cvalo = (double*)malloc(v0cnum * sizeof(double));
    for (i = 0; i < v0cnum; i++)
        sys->v0cvalo[i] = sys->v0cval[i];    // v0cvalo is the v0c values without D_sig
    for (i = 0; i < v0cnum; i++){
        sys->v0cval[i] = sys->v0cval[i] * sqrt(sys->sig[sys->v0cRowId[i]]);       // Compute the sparse form of D_sig*V0c
    }
    sys->v0cavalo = (double*)malloc(v0canum * sizeof(double));
    for (i = 0; i < v0canum; i++)
        sys->v0cavalo[i] = sys->v0caval[i];
    for (i = 0; i < v0canum; i++){
        sys->v0caval[i] = sys->v0caval[i] * sqrt(sys->sig[sys->v0caRowId[i]]);
    }


    sys->v0cColIdo = (myint*)malloc(v0cnum * sizeof(myint));
    for (i = 0; i < v0cnum; i++)
        sys->v0cColIdo[i] = sys->v0cColId[i];
    free(sys->v0cColId); sys->v0cColId = (myint*)malloc((leng_v0c + 1) * sizeof(myint));
    status = COO2CSR_malloc(sys->v0cColIdo, sys->v0cRowId, sys->v0cval, v0cnum, leng_v0c, sys->v0cColId);
    if (status != 0)
        return status;
    sys->v0caColIdo = (myint*)malloc(v0canum * sizeof(myint));
    for (i = 0; i < v0canum; i++)
        sys->v0caColIdo[i] = sys->v0caColId[i];
    free(sys->v0caColId); sys->v0caColId = (myint*)malloc((leng_v0ca + 1)*sizeof(myint));
    status = COO2CSR_malloc(sys->v0caColIdo, sys->v0caRowId, sys->v0caval, v0canum, leng_v0ca, sys->v0caColId);
    if (status != 0)
        return status;

    //status = mklMatrixMulti(sys, leng_Ac, sys->v0caRowId, sys->v0caColId, sys->v0caval, sys->N_edge, leng_v0c, sys->v0cRowId, sys->v0cColId, sys->v0cval, 2);

    cout << "leng v0c " << leng_v0c << endl;
    cout << "The number of nonzeros in Ac is " << leng_Ac << endl;
    
    /*for (i = 0; i < sys->numCdt + 1; i++){
        cout << sys->acu_cnno[i] << " ";
    }
    cout << endl;
    for (i = 0; i < sys->numCdt + 1; i++){
        cout << sys->cindex[i] << " ";
    }
    cout << endl;*/
    
    
    /*outfile1.open("Ac2.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < leng_Ac; i++){
        outfile1 << sys->AcRowId[i] + 1 << " " << sys->AcColId[i] + 1 << " " << sys->Acval[i] << endl;
    }
    outfile1.close();*/


    //outfile1.open("v0c.txt", std::ofstream::out | std::ofstream::trunc);
    //for (i = 0; i < v0cnum; i++){
    //outfile1 << sys->v0cRowId[i] + 1 << " " << sys->v0cColId[i] + 1 << " " << sys->v0cvalo[i] << endl;
    //}
    //outfile1.close();
    /*outfile1.open("v0ca.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0canum; i++){
    outfile1 << sys->v0caRowId[i] + 1 << " " << sys->v0caColIdo[i] + 1 << " " << sys->v0cavalo[i] << endl;
    }
    outfile1.close();*/

    /* Compute the matrix V0c'*D_sig*V0c */
    /*status = matrixMulti(sys->v0caColId, sys->v0caRowId, sys->v0caval, sys->v0cRowId, sys->v0cColId, sys->v0cval, sys->AcRowId, sys->AcColId, sys->Acval);
    if (status != 0)
    return status;*/


    /*outfile1.open("v0d1.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d1num; i++){
        outfile1 << sys->v0d1RowId[i] + 1 << " " << sys->v0d1ColId[i] + 1 << " " << sys->v0d1val[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0d1a.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d1anum; i++){
        outfile1 << sys->v0d1aRowId[i] + 1 << " " << sys->v0d1aColId[i] + 1 << " " << sys->v0d1aval[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0d2.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d2num; i++){
        outfile1 << sys->v0d2RowId[i] + 1 << " " << sys->v0d2ColId[i] + 1 << " " << sys->v0d2val[i] << endl;
    }
    outfile1.close();
    outfile1.open("v0d2a.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < v0d2anum; i++){
        outfile1 << sys->v0d2aRowId[i] + 1 << " " << sys->v0d2aColId[i] + 1 << " " << sys->v0d2aval[i] << endl;
    }
    outfile1.close();*/
    

    

   /* double *d = &(sys->Adval[0]);
    int *id = &(sys->AdRowId1[0]);
    int *jd = &(sys->AdColId[0]);*/

    double *a;
    int *ia, *ja;
    

   
    /* Pick up a y0c2 cooresponding to one source port */
    ofstream outfile;
    complex<double> Zresult;
    double *bd1, *bd2;
    double *bdc1, *bdc2;
    double *temp;
    double mdone;
    int ione;
    int *ipiv;
    int info;
    double *workspace;
    double *xd2;
    double *temp1;
    int startCol;
    startCol = 0;
    sys->Y = (complex<double>*)calloc(sys->numPorts * sys->numPorts * sys->nfreq, sizeof(complex<double>));
    sys->x = (complex<double>*) calloc(sys->numPorts * sys->numPorts, sizeof(complex<double>));
    for (i = 0; i < sys->numPorts*sys->numPorts; i++){
        sys->x[i] = complex<double>(0., 0.); // Complex double constructor from real and imaginary
    }
    
    
    double *b, *xc;    // the array of the right hand side
    xcol = 0;

    vector<int> v0csRowId, v0csColId;
    vector<double> v0csval;
    double *epsy, *sigy;
    double *Jr, *Ji;
    int port, sourcePort;    // show which port it is using
    int node;
    double v0csedgeleng;
    complex<double> current;
    int port_N_node_s;
    double *crhs;
    double *yc_eps;
    

    

    char transa;
    int m;
    double alpha = 1, beta = 0;
    char matdescra[6];
    matdescra[0] = 'G'; matdescra[3] = 'C';    // general matrix multi, 0-based indexing
    cout << "Begin!\n";
    double *ydcp;
    double *y0c, *y0cs;
    double *yc;
    double *yccp;
    double *dRhs, *dRhs2, *crhss;
    complex<double> *yd;
    double *yd1, *yd2;
    double *v0caJ, *v0daJ, *y0d, *ydt, *y0d2;
    double leng;
    double *u0d, *u0c;
    double nn;
    
    sourcePort = 0;
    struct matrix_descr descr;


    /* V0ca^T's csr form handle for MKL */
    sparse_matrix_t V0cat;
    s = mkl_sparse_d_create_csr(&V0cat, SPARSE_INDEX_BASE_ZERO, leng_v0c, sys->N_edge, &sys->v0caColId[0], &sys->v0caColId[1], sys->v0caRowId, sys->v0cavalo);

    /* V0c^T's csr form handle for MKL */
    sparse_matrix_t V0ct;
    s = mkl_sparse_d_create_csr(&V0ct, SPARSE_INDEX_BASE_ZERO, leng_v0c, sys->N_edge, &sys->v0cColId[0], &sys->v0cColId[1], sys->v0cRowId, sys->v0cvalo);

    

    
    while (sourcePort < sys->numPorts){

        sys->J = (double*)calloc(sys->N_edge, sizeof(double));
        for (i = 0; i < sys->portEdge[sourcePort].size(); i++){
            sys->J[sys->portEdge[sourcePort][i]] = sys->portCoor[sourcePort].portDirection;
            //cout << sys->portEdge[sourcePort][i] - sys->N_edge_s << " ";
        }
        //cout << endl;
        

        dRhs = (double*)malloc(sys->N_edge*sizeof(double));
        for (i = 0; i < sys->N_edge; i++){
            dRhs[i] = -sys->J[i];
        }

        v0daJ = (double*)calloc(leng_v0d1, sizeof(double));
        y0d = (double*)calloc(leng_v0d1, sizeof(double));
        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, dRhs, beta, v0daJ);

        /*outfile.open("rhs.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < leng_v0d1; i++){
            outfile << v0daJ[i] << endl;
        }
        outfile.close();*/

        /* solve V0d system */
        t1 = clock();
        status = hypreSolve(sys, ad, parcsr_ad, leng_Ad, v0daJ, leng_v0d1, y0d);
        /*status = hypreSolve(sys, sys->AdRowId, sys->AdColId, sys->Adval, leng_Ad, v0daJ, leng_v0d1, y0d);*/
        cout << "HYPRE solve time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
        /* End of solving */

        /*t1 = clock();
        status = solveV0dSystem(sys, v0daJ, y0d, leng_v0d1);
        cout << "Pardiso solve time " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;*/
        /*outfile1.open("y0d.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < leng_v0d1; i++){
            outfile1 << y0d[i] << " ";
        }
        outfile1.close();*/

        for (i = 0; i < leng_v0d1; i++){
            y0d[i] = y0d[i] / (2 * PI * sys->freqStart * sys->freqUnit);    // y0d is imaginary
        }
        ydt = (double*)calloc(sys->N_edge, sizeof(double));
        yd1 = (double*)malloc(sys->N_edge * sizeof(double));
        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, y0d, beta, ydt);
        
        u0d = (double*)malloc(sys->N_edge * sizeof(double));
        nn = 0;
        for (i = 0; i < sys->N_edge; i++){
            nn += ydt[i] * ydt[i];
        }
        nn = sqrt(nn);
        for (i = 0; i < sys->N_edge; i++){
            u0d[i] = ydt[i] / nn;
        }

        /* Compute C right hand side */
        y0c = (double*)calloc(leng_v0c, sizeof(double));
        v0caJ = (double*)calloc(leng_v0c, sizeof(double));
        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, sys->J, beta, v0caJ);
        for (i = 0; i < leng_v0c; i++){
            v0caJ[i] = -v0caJ[i];
        }
        crhs = (double*)calloc(leng_v0c, sizeof(double));
        for (i = 0; i < sys->N_edge; i++){
            yd1[i] = ydt[i];
            ydt[i] = -ydt[i] * (2 * PI*sys->freqStart * sys->freqUnit)*sys->eps[i];
        }

        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0cat, descr, ydt, beta, crhs);

        for (i = 0; i < leng_v0c; i++){
            v0caJ[i] = (v0caJ[i] + crhs[i]);
            
        }
        
        

        /*solve c system block by block*/
        t1 = clock();
        for (i = 1; i <= 1; i++){

            //pardisoSolve_c(sys, &v0caJ[sys->acu_cnno[i - 1]], &y0c[sys->acu_cnno[i - 1]], sys->acu_cnno[i - 1], sys->acu_cnno[i] - 1, sys->cindex[i - 1] + 1, sys->cindex[i]);

            //a = (double*)malloc((sys->cindex[i] - sys->cindex[i - 1]) * sizeof(double));
            //ia = (int*)malloc((sys->cindex[i] - sys->cindex[i - 1]) * sizeof(int));
            //ja = (int*)malloc((sys->cindex[i] - sys->cindex[i - 1]) * sizeof(int));
            //crhss = (double*)malloc((sys->acu_cnno[i] - sys->acu_cnno[i - 1]) * sizeof(double));
            //y0cs = (double*)calloc((sys->acu_cnno[i] - sys->acu_cnno[i - 1]), sizeof(double));

            //for (j = sys->acu_cnno[i - 1]; j <= sys->acu_cnno[i] - 1; j++){
            //    crhss[j - sys->acu_cnno[i - 1]] = v0caJ[j];
            //}
            //for (j = sys->cindex[i - 1] + 1; j <= sys->cindex[i]; j++){
            //    a[j - sys->cindex[i - 1] - 1] = sys->Acval[j];
            //    ia[j - sys->cindex[i - 1] - 1] = sys->AcRowId[j] - sys->AcRowId[sys->cindex[i - 1] + 1];
            //    ja[j - sys->cindex[i - 1] - 1] = sys->AcColId[j] - sys->AcColId[sys->cindex[i - 1] + 1];
            //    //cout << a[j - sys->cindex[i - 1] - 1] << " " << ia[j - sys->cindex[i - 1] - 1] << " " << ja[j - sys->cindex[i - 1] - 1] << endl;
            //}

            //status = hypreSolve(sys, ia, ja, a, (sys->cindex[i] - sys->cindex[i - 1]), crhss, (sys->acu_cnno[i] - sys->acu_cnno[i - 1]), y0cs);

            //for (j = sys->acu_cnno[i - 1]; j <= sys->acu_cnno[i] - 1; j++){
            //    y0c[j] = y0cs[j - sys->acu_cnno[i - 1]];
            //}

            //free(a); a = NULL;
            //free(ia); ia = NULL;
            //free(ja); ja = NULL;
            //free(crhss); crhss = NULL;
            //free(y0cs); y0cs = NULL;

            status = hypreSolve(sys, ac, parcsr_ac, leng_Ac, v0caJ, leng_v0c, y0c);
            //status = hypreSolve(sys, sys->AcRowId, sys->AcColId, sys->Acval, leng_Ac, v0caJ, leng_v0c, y0c);

        }
        cout << "HYPRE solve time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
        


        /* V0cy0c */
        yc = (double*)calloc(sys->N_edge, sizeof(double));
        yccp = (double*)malloc(sys->N_edge * sizeof(double));
        dRhs2 = (double*)calloc(leng_v0d1, sizeof(double));
        y0d2 = (double*)calloc(leng_v0d1, sizeof(double));
        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0ct, descr, y0c, beta, yc);

        /*outfile1.open("yc.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < leng_v0c; i++){
            outfile1 << y0c[i] << endl;
        }
        outfile1.close();*/

        for (i = 0; i < sys->N_edge; i++){
            yccp[i] = -yc[i] * sys->eps[i];

        }
        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, V0dat, descr, yccp, beta, dRhs2);

        t1 = clock();
        status = hypreSolve(sys, ad, parcsr_ad, leng_Ad, dRhs2, leng_v0d1, y0d2);
        //status = hypreSolve(sys, sys->AdRowId, sys->AdColId, sys->Adval, leng_Ad, dRhs2, leng_v0d1, y0d2);
        cout << "HYPRE solve time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;


        /*outfile1.open("y0d.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < leng_v0d1; i++){
            outfile1 << sqrt(pow(y0d[i], 2) + pow(y0d2[i], 2)) << endl;
        }
        outfile1.close();*/

        yd2 = (double*)calloc(sys->N_edge, sizeof(double));
        
        alpha = 1;
        beta = 0;
        descr.type = SPARSE_MATRIX_TYPE_GENERAL;
        s = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE, alpha, V0dt, descr, y0d2, beta, yd2);

        u0c = (double*)malloc(sys->N_edge * sizeof(double));
        nn = 0;
        for (i = 0; i < sys->N_edge; i++){
            nn += (yd2[i] + yc[i]) * (yd2[i] + yc[i]);
        }
        nn = sqrt(nn);
        for (i = 0; i < sys->N_edge; i++){
            u0c[i] = (yd2[i] + yc[i]) / nn;
        }

        
        yd = (complex<double>*)malloc(sys->N_edge * sizeof(complex<double>));
        for (int id = 0; id < sys->N_edge; id++){
            yd[id] = yd2[id] - (1i)*(yd1[id]);
        }


        sys->y = (complex<double>*)malloc(sys->N_edge*sizeof(complex<double>));
        for (i = 0; i < sys->N_edge; i++){
            sys->y[i] = yd[i].real() + yc[i] + (1i) * yd[i].imag();
        }
        /*outfile1.open("y.txt", std::ofstream::out | std::ofstream::trunc);
        for (i = 0; i < sys->N_edge; i++){
        outfile1 << sys->y[i] << endl;
        }
        outfile1.close();*/

        for (i = 0; i < sys->numPorts; i++){
            for (j = 0; j < sys->portEdge[i].size(); j++){
                leng = pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3]), 2);
                leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3 + 1]), 2);
                leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3 + 2]), 2);
                leng = sqrt(leng);
                sys->x[i + sys->numPorts*xcol] = (sys->x[i + sys->numPorts*xcol].real() + sys->y[sys->portEdge[i][j]].real() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection))) + (1i)*(sys->y[sys->portEdge[i][j]].imag() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection)) + sys->x[i + sys->numPorts*xcol].imag());

            }
        }

        //status = find_Vh(sys, u0d, u0c, sourcePort);
        

        free(ydt); ydt = NULL;
        free(y0c); y0c = NULL;
        free(yd); yd = NULL;
        free(ydcp); ydcp = NULL;
        free(yc); yc = NULL;
        free(yccp); yccp = NULL;
        free(v0caJ); v0caJ = NULL;
        free(sys->y); sys->y = NULL;
        free(dRhs); dRhs = NULL;
        free(sys->J); sys->J = NULL;
        free(crhs); crhs = NULL;
        free(u0d); u0d = NULL;
        free(u0c); u0c = NULL;

        sourcePort++;
        xcol++;
        
    }
    MPI_Finalize();
    if (sys->nfreq > 1){
        for (int id = 0; id < sys->nfreq; id++){
            cout << "Z parameter at frequency " << (sys->freqStart + id * (sys->freqEnd - sys->freqStart) / (sys->nfreq - 1)) * sys->freqUnit << " is " << endl;
            for (i = 0; i < sys->numPorts; i++){
                for (j = 0; j < sys->numPorts; j++){
                    Zresult = sys->x[j + i*sys->numPorts].real() + (1i) * sys->x[j + i*sys->numPorts].imag() * sys->freqStart / (sys->freqStart + id * (sys->freqEnd - sys->freqStart) / (sys->nfreq - 1));
                    cout << Zresult << "\n";
                }
            }
        }
    }
    else{
        cout << "Z parameter at frequency " << (sys->freqStart) * sys->freqUnit << " is " << endl;
        for (i = 0; i < sys->numPorts; i++){
            for (j = 0; j < sys->numPorts; j++){
                Zresult = sys->x[j + i*sys->numPorts].real() + (1i) * sys->x[j + i*sys->numPorts].imag();
                cout << Zresult << "\n";
            }
        }
    }

    free(sys->cindex); sys->cindex = NULL;
    free(sys->acu_cnno); sys->acu_cnno = NULL;
    free(sys->x); sys->x = NULL;

    free(sys->AdColId); sys->AdColId = NULL;
    free(sys->Adval); sys->Adval = NULL;
    free(sys->AdRowId); sys->AdRowId = NULL;
    free(sys->v0d1RowId); sys->v0d1RowId = NULL;
    free(sys->v0d1ColId); sys->v0d1ColId = NULL;
    free(sys->v0d1ColIdo); sys->v0d1ColIdo = NULL;
    free(sys->v0d1val); sys->v0d1val = NULL;
    free(sys->v0d1valo); sys->v0d1valo = NULL;
    free(sys->v0d1aRowId); sys->v0d1aRowId = NULL;
    free(sys->v0d1aColId); sys->v0d1aColId = NULL;
    free(sys->v0d1aColIdo); sys->v0d1aColIdo = NULL;
    free(sys->v0d1aval); sys->v0d1aval = NULL;
    free(sys->v0d1avalo); sys->v0d1avalo = NULL;
    
    free(sys->v0cRowId);  sys->v0cRowId = NULL;
    free(sys->v0cColId);  sys->v0cColId = NULL;
    free(sys->v0cColIdo); sys->v0cColIdo = NULL;
    free(sys->v0cval); sys->v0cval = NULL;
    free(sys->v0cvalo); sys->v0cvalo = NULL;
    free(sys->v0caRowId); sys->v0caRowId = NULL;
    free(sys->v0caColId); sys->v0caColId = NULL;
    free(sys->v0caColIdo); sys->v0caColIdo = NULL;
    free(sys->v0caval); sys->v0caval = NULL;
    free(sys->v0cavalo); sys->v0cavalo = NULL;
    free(sys->AcRowId); sys->AcRowId = NULL;
    free(sys->AcColId); sys->AcColId = NULL;
    free(sys->Acval); sys->Acval = NULL;

    mkl_sparse_destroy(V0dt);
    mkl_sparse_destroy(V0dat);
    mkl_sparse_destroy(V0ct);
    mkl_sparse_destroy(V0cat);

    HYPRE_IJMatrixDestroy(ad);
    HYPRE_IJMatrixDestroy(ac);

    return 0;
}
//int pardisoSolve_c(fdtdMesh *sys, double *rhs, double *solution, int nodestart, int nodeend, int indstart, int indend){
//
//    /*solve Ac system block by block*/
//    int leng = nodeend - nodestart + 1;
//    int status;
//    double *a = (double*)malloc((indend - indstart + 1) * sizeof(double));
//    int *ia = (int*)malloc((indend - indstart + 1) * sizeof(int));
//    int *ja = (int*)malloc((indend - indstart + 1) * sizeof(int));
//
//    for (int i = 0; i < (indend - indstart + 1); i++){
//        a[i] = sys->Acval[indstart + i];
//        ia[i] = sys->AcRowId[indstart + i] - sys->AcRowId[indstart];
//        ja[i] = sys->AcColId[indstart + i] - sys->AcColId[indstart];
//    }
//
//    /* use pardiso to solve */
//    {
//         int *ia1 = (int*)malloc((leng + 1) * sizeof(int));
//        /* status = COO2CSR_malloc(ia, ja, a, indstart - indend + 1, leng, ia1);
//         if (status != 0)
//             return status;*/
//         int count = 0;
//         int i = 0;
//         int k = 0;
//         int start;
//         ia1[k] = 0;
//         k++;
//         while (i < (indend - indstart + 1)){
//             start = ia[i];
//             while (i < (indend - indstart + 1) && ia[i] == start) {
//                 count++;
//                 i++;
//             }
//             ia1[k] = (count);
//             k++;
//         }
//
//         void *pt[64];
//         int mtype;
//         int iparm[64];
//         double dparm[64];
//         int maxfct, mnum, phase, error, solver;
//         int num_process;   //number of processors
//         int v0csin;
//         int perm;
//         int nrhs = 1;
//         int msglvl = 0;    //print statistical information
//
//         mtype = 11;    // real and not symmetric
//         solver = 0;
//         error = 0;
//         maxfct = 1;    //maximum number of numerical factorizations
//         mnum = 1;    //which factorization to use
//         phase = 13;    //analysis
//
//         pardisoinit(pt, &mtype, iparm);
//         nrhs = 1;
//         iparm[38] = 1;
//         iparm[34] = 1;    //0-based indexing
//         
//         pardiso(pt, &maxfct, &mnum, &mtype, &phase, &leng, a, ia1, ja, &perm, &nrhs, iparm, &msglvl, rhs, solution, &error);
//         free(ia1); ia1 = NULL;
//    }
//
//    /* use HYPRE to solve */
//    /*{
//        status = hypreSolve(sys, ia, ja, a, indend - indstart + 1, rhs, leng, solution);
//    }*/
//    
//    
//    free(a); a = NULL;
//    free(ia); ia = NULL;
//    free(ja); ja = NULL;
//    
//    
//
//    return 0;
//}



//int interativeSolver(int N, int nrhs, double *rhs, int *ia, int *ja, double *a, int *ib, int *jb, double *b, double *solution, fdtdMesh *sys){
//    // ia, ja, a are CSR form with one-based indexing and it is only upper triangular elements (symmetric)
//    double mdone = -1;
//    int ione = 1;
//    double euclidean_norm;
//    int RCI_request, itercount;
//    int *ipar;
//    double *dpar;
//    double *tmp;
//    char U = 'U';
//    double *temp;
//    vector<char> transa;
//    int m;    // number of rows
//    int k;    // number of columns
//    double alpha = 1, beta = 0;
//    char matdescra[6];
//    int *pntrb1, *pntre1;
//    int *pntrb2, *pntre2;
//    int i, j;
//    double *y;
//    int length = 128;
//    ipar = (int*)malloc((length + 2 * nrhs) * sizeof(int));
//    dpar = (double*)malloc((length + 2 * nrhs) * sizeof(double));
//    tmp = (double*)malloc(N*(3 + nrhs) * sizeof(double));
//    y = (double*)malloc(sys->N_edge * sizeof(double));
//    temp = (double*)malloc(N * sizeof(double));
//    transa.push_back('T'); transa.push_back('N');
//    m = N;
//    k = sys->N_edge;
//    matdescra[0] = 'G'; matdescra[3] = 'C';
//    pntrb1 = ia;
//    pntre1 = &ia[1];
//    pntrb2 = ib;
//    pntre2 = &ib[1];
//    alpha = 1; beta = 0;
//
//    dfgmres_init(&N, solution, rhs, &RCI_request, ipar, dpar, tmp);
//    ipar[14] = 2;
//    //ipar[7] = 0;
//    ipar[10] = 0;
//    dpar[0] = 5e-3;
//    ipar[4] = 1000;   // set the maximum iteration number
//    dfgmres_check(&N, solution, rhs, &RCI_request, ipar, dpar, tmp);
//
//    while (true){
//        dfgmres(&N, solution, rhs, &RCI_request, ipar, dpar, tmp);
//        //cout << RCI_request << " ";
//        if (RCI_request == 0){
//            break;
//        }
//        if (RCI_request == 1){
//            mkl_dcsrmv(&transa[0], &m, &k, &alpha, matdescra, a, ja, pntrb1, pntre1, tmp, &beta, y);
//            mkl_dcsrmv(&transa[1], &m, &k, &alpha, matdescra, b, jb, pntrb2, pntre2, y, &beta, &tmp[N]);
//            continue;
//        }
//        if (RCI_request == 2){
//            for (j = 0; j < nrhs; j++){
//                mkl_dcsrmv(&transa[0], &m, &k, &alpha, matdescra, a, ja, pntrb1, pntre1, solution, &beta, y);
//                mkl_dcsrmv(&transa[1], &m, &k, &alpha, matdescra, b, jb, pntrb2, pntre2, y, &beta, temp);
//            }
//            cblas_daxpy(N, mdone, rhs, ione, temp, ione);
//            euclidean_norm = cblas_dnrm2(N, temp, ione) / cblas_dnrm2(N, rhs, ione);
//            //cout << euclidean_norm << "\n";
//            if (euclidean_norm > 1.e-3)
//                continue;
//            else
//                break;
//        }
//        else{
//            if (dpar[6] < 1e-12) break;
//            else continue;
//        }
//    }
//
//    dfgmres_get(&N, solution, rhs, &RCI_request, ipar, dpar, tmp, &itercount);
//    //cout << itercount << "\n";
//    //cout << euclidean_norm << "\n";
//    free(tmp);
//    free(temp);
//    free(ipar);
//    free(dpar);
//    free(y);
//
//    return 0;
//}

//int mklMatrixMulti(fdtdMesh *sys, int &leng_A, int *aRowId, int *aColId, double *aval, int arow, int acol, int *bRowId, int *bColId, double *bval, int mark){
//    // ArowId, AcolId, and Aval should be in the COO format
//    sparse_status_t s0;
//    sparse_matrix_t a, a_csr;
//    sparse_index_base_t indexing1 = SPARSE_INDEX_BASE_ZERO;
//    int row, col;
//    int *cols, *cole, *rowi;
//    double *val;
//    MKL_INT *AcolId;
//    double *Aval;
//    int k;
//
//    s0 = mkl_sparse_d_create_csc(&a, SPARSE_INDEX_BASE_ZERO, arow, acol, &aColId[0], &aColId[1], aRowId, aval);
//
//
//    sparse_matrix_t b, b_csr;
//    s0 = mkl_sparse_d_create_csc(&b, SPARSE_INDEX_BASE_ZERO, arow, acol, &bColId[0], &bColId[1], bRowId, bval);
//    
//    s0 = mkl_sparse_convert_csr(a, SPARSE_OPERATION_NON_TRANSPOSE, &a_csr);
//
//    s0 = mkl_sparse_convert_csr(b, SPARSE_OPERATION_NON_TRANSPOSE, &b_csr);
//
//    sparse_matrix_t A;
//    s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, a_csr, b_csr, &A);
//    
//
//    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
//    int ARows, ACols;
//    MKL_INT *ArowStart, *ArowEnd;
//    s0 = mkl_sparse_d_export_csr(A, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
//    
//    leng_A = ArowEnd[ARows - 1];
//    
//    if (mark == 1){    // dielectric
//        sys->AdRowId = (int*)malloc(leng_A * sizeof(int));
//        sys->AdColId = (int*)malloc(leng_A * sizeof(int));
//        sys->Adval = (double*)malloc(leng_A * sizeof(double));
//        int count, num, j;
//
//        j = 0;
//        for (int i = 0; i < ARows; i++){
//            num = ArowEnd[i] - ArowStart[i];
//            count = 0;
//            while (count < num){
//                sys->AdRowId[j] = i;
//                sys->AdColId[j] = AcolId[j];
//                sys->Adval[j] = Aval[j];
//                j++;
//                count++;
//            }
//        }
//    }
//    else if (mark == 2){    // conductor
//        sys->AcRowId = (int*)malloc(leng_A * sizeof(int));
//        sys->AcColId = (int*)malloc(leng_A * sizeof(int));
//        sys->Acval = (double*)malloc(leng_A * sizeof(double));
//        int count, num, j;
//
//        j = 0;
//        k = 1;
//        sys->cindex[k] = sys->cindex[k - 1];
//        for (int i = 0; i < ARows; i++){
//            num = ArowEnd[i] - ArowStart[i];
//            count = 0;
//            while (count < num){
//                sys->AcRowId[j] = i;
//                sys->AcColId[j] = AcolId[j];
//                sys->Acval[j] = Aval[j];
//                if (k - 1 <= i){
//                    sys->cindex[k]++;
//                }
//                else{
//                    k++;
//                    sys->cindex[k] = sys->cindex[k - 1];
//                    sys->cindex[k]++;
//                }
//                j++;
//                count++;
//            }
//        }
//    }
//    
//    mkl_sparse_destroy(a);
//    mkl_sparse_destroy(b);
//    mkl_sparse_destroy(A);
//    
//   
//    return 0;
//}

//int matrixMulti(int *aRowId, int *aColId, double *aval, int anum, int *bRowId, int *bColId, double *bval, int bnum, fdtdMesh *sys, int &leng, int mark){
//    //the first matrix is row by row, the second matrix is column by column
//    // leng is used to record the number of elements in the derived matrix
//
//    int i = 0, j = 0;
//    int flaga, flagb, k;
//    int starta;
//    double sum;
//    /*vector<int> cRowId;
//    vector<int> cColId;
//    vector<double> cval;*/
//    
//
//    flaga = aRowId[0];
//    flagb = bColId[0];
//    starta = 0;
//    sum = 0;
//    while (i < anum){
//        while (j < bnum && bColId[j] == flagb && i < anum && aRowId[i] == flaga){
//            if (aColId[i] == bRowId[j]){
//                sum += aval[i] * bval[j];
//                j++;
//                i++;
//            }
//            else if (aColId[i] < bRowId[j]){
//                i++;
//            }
//            else if (aColId[i] > bRowId[j]){
//                j++;
//            }
//        }
//        if (sum != 0){
//            /*cRowId.push_back(flaga);
//            cColId.push_back(flagb);
//            cval.push_back(sum);*/
//            sum = 0;
//            leng++;
//        }
//        if (i == anum){
//            if (j == bnum)
//                break;
//            else{
//                i = starta;
//                while (bColId[j] == flagb){
//                    j++;
//                    if (j == bnum){
//                        while (i < anum && aRowId[i] == flaga){
//                            i++;
//                        }
//                        starta = i;
//                        if (i == anum)    //run all of the datas
//                            break;
//                        flaga = aRowId[i];
//                        j = 0;
//                        break;
//                    }
//                }
//                flagb = bColId[j];
//                continue;
//            }
//        }
//        if (j == bnum){
//            while (i < anum && aRowId[i] == flaga){
//                i++;
//            }
//            starta = i;
//            if (i == anum)    //run all of the datas
//                break;
//            flaga = aRowId[i];
//            j = 0;
//        }
//        else{
//            if (bColId[j] != flagb && aRowId[i] != flaga){
//                flagb = bColId[j];
//                i = starta;
//            }
//            else if (bColId[j] != flagb){
//                flagb = bColId[j];
//                i = starta;
//            }
//            else if (aRowId[i] != flaga){
//                i = starta;
//                while (bColId[j] == flagb){
//                    j++;
//                    if (j == bnum){
//                        while (i < anum && aRowId[i] == flaga){
//                            i++;
//                        }
//                        starta = i;
//                        if (i == anum)    //run all of the datas
//                            break;
//                        flaga = aRowId[i];
//                        j = 0;
//                        break;
//                    }
//                }
//                flagb = bColId[j];
//            }
//
//        }
//    }
//    if (mark == 1){
//        sys->AdRowId = (int*)malloc(leng * sizeof(int));
//        sys->AdColId = (int*)malloc(leng * sizeof(int));
//        sys->Adval = (double*)malloc(leng * sizeof(double));
//    }
//    else if (mark == 2) {
//        sys->AcRowId = (int*)malloc(leng * sizeof(int));
//        sys->AcColId = (int*)malloc(leng * sizeof(int));
//        sys->Acval = (double*)malloc(leng * sizeof(double));
//    }
//    leng = 0;
//
//    flaga = aRowId[0];
//    flagb = bColId[0];
//    starta = 0;
//    sum = 0;
//    i = 0;
//    j = 0;
//    int count = 1;
//    sys->cindex[count] = sys->cindex[count - 1];
//    while (i < anum){
//        while (j < bnum && bColId[j] == flagb && i < anum && aRowId[i] == flaga){
//            if (aColId[i] == bRowId[j]){
//                sum += aval[i] * bval[j];
//                j++;
//                i++;
//            }
//            else if (aColId[i] < bRowId[j]){
//                i++;
//            }
//            else if (aColId[i] > bRowId[j]){
//                j++;
//            }
//        }
//        if (sum != 0){
//            if (mark == 1){
//                sys->AdRowId[leng] = (flaga);
//                sys->AdColId[leng] = (flagb);
//                sys->Adval[leng] = (sum);
//            }
//            else if (mark == 2){
//                sys->AcRowId[leng] = (flaga);
//                sys->AcColId[leng] = (flagb);
//                sys->Acval[leng] = (sum);
//                if (flaga < sys->acu_cnno[count]){
//                    sys->cindex[count]++;
//                }
//                else{
//                    count++;
//                    sys->cindex[count] = sys->cindex[count - 1];
//                    sys->cindex[count]++;
//                }
//            }
//            sum = 0;
//            leng++;
//        }
//        if (i == anum){
//            if (j == bnum)
//                break;
//            else{
//                i = starta;
//                while (bColId[j] == flagb){
//                    j++;
//                    if (j == bnum){
//                        while (i < anum && aRowId[i] == flaga){
//                            i++;
//                        }
//                        starta = i;
//                        if (i == anum)    //run all of the datas
//                            break;
//                        flaga = aRowId[i];
//                        j = 0;
//                        break;
//                    }
//                }
//                flagb = bColId[j];
//                continue;
//            }
//        }
//        if (j == bnum){
//            while (i < anum && aRowId[i] == flaga){
//                i++;
//            }
//            starta = i;
//            if (i == anum)    //run all of the datas
//                break;
//            flaga = aRowId[i];
//            j = 0;
//        }
//        else{
//            if (bColId[j] != flagb && aRowId[i] != flaga){
//                flagb = bColId[j];
//                i = starta;
//            }
//            else if (bColId[j] != flagb){
//                flagb = bColId[j];
//                i = starta;
//            }
//            else if (aRowId[i] != flaga){
//                i = starta;
//                while (bColId[j] == flagb){
//                    j++;
//                    if (j == bnum){
//                        while (i < anum && aRowId[i] == flaga){
//                            i++;
//                        }
//                        starta = i;
//                        if (i == anum)    //run all of the datas
//                            break;
//                        flaga = aRowId[i];
//                        j = 0;
//                        break;
//                    }
//                }
//                flagb = bColId[j];
//            }
//        }
//    }
//    
//    return 0;
//}

int matrixMul(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> bRowId, vector<int> bColId, vector<double> bval, vector<int> &cRowId, vector<int> &cColId, vector<double> &cval){
    //the first matrix is row by row, the second matrix is column by column

    int i = 0, j = 0;
    int flaga, flagb, k;
    int starta;
    double sum;

    flaga = aRowId[0];
    flagb = bColId[0];
    starta = 0;
    sum = 0;
    while (i < aRowId.size()){
        while (j < bColId.size() && bColId[j] == flagb && i < aRowId.size() && aRowId[i] == flaga){
            if (aColId[i] == bRowId[j]){
                sum += aval[i] * bval[j];
                j++;
                i++;
            }
            else if (aColId[i] < bRowId[j]){
                i++;
            }
            else if (aColId[i] > bRowId[j]){
                j++;
            }
        }
        if (sum != 0){
            cRowId.push_back(flaga);
            cColId.push_back(flagb);
            cval.push_back(sum);
            sum = 0;
        }
        if (i == aRowId.size()){
            if (j == bColId.size())
                break;
            else{
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bColId.size()){
                        while (i < aRowId.size() && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == aRowId.size())    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
                continue;
            }
        }
        if (j == bColId.size()){
            while (i < aRowId.size() && aRowId[i] == flaga){
                i++;
            }
            starta = i;
            if (i == aRowId.size())    //run all of the datas
                break;
            flaga = aRowId[i];
            j = 0;
        }
        else{
            if (bColId[j] != flagb && aRowId[i] != flaga){
                flagb = bColId[j];
                i = starta;
            }
            else if (bColId[j] != flagb){
                flagb = bColId[j];
                i = starta;
            }
            else if (aRowId[i] != flaga){
                i = starta;
                while (bColId[j] == flagb){
                    j++;
                    if (j == bColId.size()){
                        while (i < aRowId.size() && aRowId[i] == flaga){
                            i++;
                        }
                        starta = i;
                        if (i == aRowId.size())    //run all of the datas
                            break;
                        flaga = aRowId[i];
                        j = 0;
                        break;
                    }
                }
                flagb = bColId[j];
            }

        }
    }

    return 0;
}


int COO2CSR(vector<int> &rowId, vector<int> &ColId, vector<double> &val){
    int i;
    vector<int> rowId2;
    int count, start;

    rowId2.push_back(0);
    count = 0;
    i = 0;
    while (i < rowId.size()){
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

int COO2CSR_malloc(myint *rowId, myint *ColId, double *val, myint totalnum, myint leng, myint *rowId1){    // totalnum is the total number of entries, leng is the row number
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
    while (i < totalnum){
        start = rowId[i];
        while (i < totalnum && rowId[i] == start) {
            count++;
            i++;
        }
        rowId2[k] = (count);
        k++;
    }

    for (i = 0; i <= leng; i++){
        rowId1[i] = rowId2[i];
    }

    free(rowId2); rowId2 = NULL;
    return 0;
}

int mvMulti(vector<int> aRowId, vector<int> aColId, vector<double> aval, vector<int> &bRowId, vector<int> &bColId, vector<double> &bval, double *index_val, int size){
    //the same sequence in aColId and index
    double *v;
    int i;

    i = 0;
    v = (double*)calloc(size, sizeof(double));
    while (i < aColId.size()){
        v[aRowId[i]] += index_val[aColId[i]] * aval[i];
        i++;
    }
    for (i = 0; i < size; i++){
        if (abs(v[i]) > 1.e-1){
            bRowId.push_back(i);
            bColId.push_back(0);
            bval.push_back(v[i]);
        }
    }

    return 0;
}

int nodeAdd_count(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng){    // size is the size of index, total_size is the size of the vector

    int i, j;
    double *v;

    v = (double*)calloc(total_size, sizeof(double));
    for (i = 0; i < size; i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }
    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            num++;
        }
    }
    leng++;

    free(v);
    v = NULL;
    return 0;

}

int nodeAdd(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int mark){    // size is the size of index, total_size is the size of the vector

    int i, j;
    double *v;

    v = (double*)calloc(total_size, sizeof(double));
    for (i = 0; i < size; i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }
    if (mark == 1){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d1RowId[num] = (i);
                sys->v0d1ColId[num] = (leng);
                sys->v0d1val[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 2){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d2RowId[num] = (i);
                sys->v0d2ColId[num] = (leng);
                sys->v0d2val[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 3){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0cRowId[num] = (i);
                sys->v0cColId[num] = (leng);
                sys->v0cval[num] = (v[i]);
                
                num++;
            }
        }
        leng++;
    }

    free(v);
    v = NULL;
    return 0;

}

int nodeAddLarger(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int *RowId, int *ColId, double *Val){    // size is the size of index, total_size is the size of the vector

    /* original code */
    /*int i, j;
    double *v;

    v = (double*)calloc(total_size, sizeof(double));
    for (i = 0; i < size; i++){
        for (j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            v[sys->nodeEdge[index[i]][j].first] += sys->nodeEdge[index[i]][j].second;
        }
    }

    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            RowId[num] = (i);
            ColId[num] = (leng);
            Val[num] = (v[i]);
            num++;
        }
    }
    leng++;

    free(v);
    v = NULL;
    return 0;*/

    /******************************************************************/
    unordered_map<int, double> v;
    for (int i = 0; i < size; i++){
        for (int j = 0; j < sys->nodeEdge[index[i]].size(); j++){
            if (v.find(sys->nodeEdge[index[i]][j].first) == v.end()){
                v[sys->nodeEdge[index[i]][j].first] = sys->nodeEdge[index[i]][j].second;
            }
            else{
                v.erase(sys->nodeEdge[index[i]][j].first);
            }
        }
    }
    for (auto vi : v){
        if (abs(vi.second) > 1e-5){
            RowId[num] = vi.first;
            ColId[num] = leng;
            Val[num] = vi.second;
            num++;
        }
    }
    leng++;

    return 0;
}


int nodeAddAvg_count(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng){    // Get the average V0d2 (around the conductor)
    int i, j;
    double *v;
    int inx, iny, inz;
    vector<int> rowId, colId;
    vector<double> val;

    int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < size; i++){
        nodeset[index[i]] = 1;
    }


    v = (double*)calloc(total_size, sizeof(double));
    inz = index[0] / sys->N_node_s;
    inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index[0] % (sys->N_cell_y + 1);
    if (iny == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    }
    else if (iny == sys->N_cell_y){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    }

    if (inx == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    }
    else if (inx == sys->N_cell_x){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
        colId.push_back(1);
        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    }

    if (inz == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    }
    else if (inz == sys->N_cell_z){
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    }

    for (i = 0; i < val.size(); i++){
        v[rowId[i]] = val[i];
    }

    while (!rowId.empty()){
        rowId.pop_back();
        colId.pop_back();
        val.pop_back();
    }

    int *visited;
    stack<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int mark;
    int count;
    st.push(index[0]);
    visited[index[0]] = 1;
    while (!st.empty()){
        mark = 0;
        for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
            if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
                /*cout << "h\n";*/
                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] = v[rowId[i]] + ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }

                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
                mark = 1;

                break;
            }
            else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] += ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }
                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
                mark = 1;
                break;
            }
        }
        if (mark == 0){
            st.pop();
        }
    }

    for (i = 0; i < total_size; i++){
        if (abs(v[i]) > 1e-5) {
            num++;
        }
    }
    leng++;
    free(visited);
    visited = NULL;
    free(v);
    v = NULL;
    free(nodeset);
    nodeset = NULL;

    return 0;
}

int nodeAddAvg(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int mark){    // Get the average V0d2 (around the conductor)
    int i, j;
    double *v;
    int inx, iny, inz;
    vector<int> rowId, colId;
    vector<double> val;

    int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < size; i++){
        nodeset[index[i]] = 1;
    }

    v = (double*)calloc(total_size, sizeof(double));
    inz = index[0] / sys->N_node_s;
    inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index[0] % (sys->N_cell_y + 1);
    if (iny == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    }
    else if (iny == sys->N_cell_y){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
        colId.push_back(1);
        val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    }

    if (inx == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    }
    else if (inx == sys->N_cell_x){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
        colId.push_back(1);
        val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    }

    if (inz == 0){
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    }
    else if (inz == sys->N_cell_z){
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    }
    else{
        rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
        rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
        colId.push_back(1);
        val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    }

    for (i = 0; i < val.size(); i++){
        v[rowId[i]] = val[i];
    }
    
    while (!rowId.empty()){
        rowId.pop_back();
        colId.pop_back();
        val.pop_back();
    }

    int *visited;
    stack<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int record;
    int count;
    st.push(index[0]);
    visited[index[0]] = 1;
    while (!st.empty()){
        record = 0;
        for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
            if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
                /*cout << "h\n";*/
                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] = v[rowId[i]] + ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }

                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
                record = 1;

                break;
            }
            else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

                if (iny == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (iny == sys->N_cell_y){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inx == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inx == sys->N_cell_x){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                if (inz == 0){
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else if (inz == sys->N_cell_z){
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
                    }
                }
                else{
                    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
                    colId.push_back(1);
                    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
                    if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
                        ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
                    }
                }

                for (i = 0; i < rowId.size(); i++){
                    v[rowId[i]] += ratio * val[i];
                }
                while (!rowId.empty()){
                    rowId.pop_back();
                    colId.pop_back();
                    val.pop_back();
                }
                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
                record = 1;
                break;
            }
        }
        if (record == 0){
            st.pop();
        }
    }

    // mark = 1 is for v0d1a, mark = 2 is for v0d2a, mark = 3 is for v0c
    if (mark == 1){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d1aRowId[num] = (i);
                sys->v0d1aColId[num] = (leng);
                sys->v0d1aval[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 2){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0d2aRowId[num] = (i);
                sys->v0d2aColId[num] = (leng);
                sys->v0d2aval[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }
    else if (mark == 3){
        for (i = 0; i < total_size; i++){
            if (abs(v[i]) > 1e-5) {
                sys->v0caRowId[num] = (i);
                sys->v0caColId[num] = (leng);
                sys->v0caval[num] = (v[i]);
                num++;
            }
        }
        leng++;
    }

    free(visited);
    visited = NULL;
    free(v);
    v = NULL;
    free(nodeset);
    nodeset = NULL;

    return 0;
}

int nodeAddAvgLarger(int *index, int size, int total_size, fdtdMesh *sys, int &num, int &leng, int* RowId, int *ColId, double *Val){    // Get the average V0d2 (around the conductor)
    /* orginal code */
    //int i, j;
    //double *v;
    //int inx, iny, inz;
    //vector<int> rowId, colId;
    //vector<double> val;

    //int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    //for (i = 0; i < size; i++){
    //    nodeset[index[i]] = 1;
    //}

    //v = (double*)calloc(total_size, sizeof(double));
    //inz = index[0] / sys->N_node_s;
    //inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    //iny = index[0] % (sys->N_cell_y + 1);
    //if (iny == 0){
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
    //    colId.push_back(1);
    //    val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    //}
    //else if (iny == sys->N_cell_y){
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //    colId.push_back(1);
    //    val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    //}
    //else{
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
    //    colId.push_back(1);
    //    val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //    colId.push_back(1);
    //    val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //}

    //if (inx == 0){
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    //}
    //else if (inx == sys->N_cell_x){
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    //}
    //else{
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
    //    colId.push_back(1);
    //    val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //}

    //if (inz == 0){
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    //}
    //else if (inz == sys->N_cell_z){
    //    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    //}
    //else{
    //    rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //    rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //    colId.push_back(1);
    //    val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //}

    //for (i = 0; i < val.size(); i++){
    //    v[rowId[i]] = val[i];
    //}

    //while (!rowId.empty()){
    //    rowId.pop_back();
    //    colId.pop_back();
    //    val.pop_back();
    //}

    //int *visited;
    //stack<int> st;
    //double ratio;
    //visited = (int*)calloc(sys->N_node, sizeof(int));
    //int record;
    //int count;
    //st.push(index[0]);
    //visited[index[0]] = 1;
    //while (!st.empty()){
    //    record = 0;
    //    for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
    //        if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1)){
    //            visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

    //            inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
    //            inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    //            iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
    //            /*cout << "h\n";*/
    //            if (iny == 0){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
    //                colId.push_back(1);
    //                val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //                colId.push_back(1);
    //                val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else{
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
    //                colId.push_back(1);
    //                val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //                colId.push_back(1);
    //                val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }

    //            if (inx == 0){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else if (inx == sys->N_cell_x){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else{
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
    //                colId.push_back(1);
    //                val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }

    //            if (inz == 0){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else if (inz == sys->N_cell_z){
    //                rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    //                if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else{
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //                rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //                if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }

    //            for (i = 0; i < rowId.size(); i++){
    //                v[rowId[i]] = v[rowId[i]] + ratio * val[i];
    //            }
    //            while (!rowId.empty()){
    //                rowId.pop_back();
    //                colId.pop_back();
    //                val.pop_back();
    //            }

    //            st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
    //            record = 1;

    //            break;
    //        }
    //        else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1)){
    //            visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

    //            inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
    //            inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    //            iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

    //            if (iny == 0){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny);
    //                colId.push_back(1);
    //                val.push_back(1 / (sys->yn[iny + 1] - sys->yn[iny]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->yn[iny + 1] - sys->yn[iny]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else if (iny == sys->N_cell_y){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //                colId.push_back(1);
    //                val.push_back(-1 / (sys->yn[iny] - sys->yn[iny - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*(sys->N_cell_y) + iny - 1 == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->yn[iny] - sys->yn[iny - 1]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else{
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny);
    //                colId.push_back(1);
    //                val.push_back(2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1);
    //                colId.push_back(1);
    //                val.push_back(-2 / (sys->yn[iny + 1] - sys->yn[iny - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + inx*sys->N_cell_y + iny - 1 == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->yn[iny + 1] - sys->yn[iny - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }

    //            if (inx == 0){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(1 / (sys->xn[inx + 1] - sys->xn[inx]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->xn[inx + 1] - sys->xn[inx]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else if (inx == sys->N_cell_x){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-1 / (sys->xn[inx] - sys->xn[inx - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->xn[inx] - sys->xn[inx - 1]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else{
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx - 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny);
    //                colId.push_back(1);
    //                val.push_back(2 / (sys->xn[inx + 1] - sys->xn[inx - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_y)*(sys->N_cell_x + 1) + (sys->N_cell_y + 1)*(inx)+iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->xn[inx + 1] - sys->xn[inx - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }

    //            if (inz == 0){
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(1 / (sys->zn[inz + 1] - sys->zn[inz]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->zn[inz + 1] - sys->zn[inz]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else if (inz == sys->N_cell_z){
    //                rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-1 / (sys->zn[inz] - sys->zn[inz - 1]));
    //                if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->zn[inz] - sys->zn[inz - 1]) * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }
    //            else{
    //                rowId.push_back(inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //                if (inz*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = -(sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //                rowId.push_back((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny);
    //                colId.push_back(1);
    //                val.push_back(-2 / (sys->zn[inz + 1] - sys->zn[inz - 1]));
    //                if ((inz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + inx*(sys->N_cell_y + 1) + iny == sys->nodeEdge[st.top()][j].first){
    //                    ratio = (sys->zn[inz + 1] - sys->zn[inz - 1]) / 2 * v[sys->nodeEdge[st.top()][j].first];
    //                }
    //            }

    //            for (i = 0; i < rowId.size(); i++){
    //                v[rowId[i]] += ratio * val[i];
    //            }
    //            while (!rowId.empty()){
    //                rowId.pop_back();
    //                colId.pop_back();
    //                val.pop_back();
    //            }
    //            st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
    //            record = 1;
    //            break;
    //        }
    //    }
    //    if (record == 0){
    //        st.pop();
    //    }
    //}


    //for (i = 0; i < total_size; i++){
    //    if (abs(v[i]) > 1e-5) {
    //        RowId[num] = (i);
    //        ColId[num] = (leng);
    //        Val[num] = (v[i]);
    //        num++;
    //    }
    //}
    //leng++;



    //free(visited);
    //visited = NULL;
    //free(v);
    //v = NULL;
    //free(nodeset);
    //nodeset = NULL;

    //return 0;

    /***************************************************************************************/
    int i, j;
    unordered_map<int, double> v;
    int inx, iny, inz;
    vector<int> rowId, colId;
    vector<double> val;

    int *nodeset = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < size; i++){
        nodeset[index[i]] = 1;
    }

    inz = index[0] / sys->N_node_s;
    inx = (index[0] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
    iny = index[0] % (sys->N_cell_y + 1);
    for (i = 0; i < sys->nodeEdgea[index[0]].size(); i++){
        v[sys->nodeEdgea[index[0]][i].first] = sys->nodeEdgea[index[0]][i].second;
    }


    int *visited;
    stack<int> st;
    double ratio;
    visited = (int*)calloc(sys->N_node, sizeof(int));
    int record;
    int count;
    st.push(index[0]);
    visited[index[0]] = 1;
    while (!st.empty()){
        record = 0;
        for (j = 0; j < sys->nodeEdge[st.top()].size(); j++){
            if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2] % (sys->N_cell_y + 1);
                
                for (i = 0; i < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]].size(); i++){
                    if (sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]][i].first == sys->nodeEdge[st.top()][j].first){
                        ratio = -1 / sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]][i].second * v[sys->nodeEdge[st.top()][j].first];
                        break;
                    }
                }

                for (i = 0; i < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]].size(); i++){
                    if (v.find(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]][i].first) == v.end()){
                        v[sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]][i].first] = ratio * sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]][i].second;
                    }
                    else{
                        v.erase(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]][i].first);
                    }
                }


                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2]);
                record = 1;

                break;
            }
            else if ((sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] != st.top() && visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 0) && (nodeset[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] == 1)){
                visited[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]] = 1;

                inz = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] / sys->N_node_s;
                inx = (sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] - inz * sys->N_node_s) / (sys->N_cell_y + 1);
                iny = sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1] % (sys->N_cell_y + 1);

                for (i = 0; i < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]].size(); i++){
                    if (sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]][i].first == sys->nodeEdge[st.top()][j].first){
                        ratio = -1 / sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]][i].second * v[sys->nodeEdge[st.top()][j].first];
                        break;
                    }
                }

                for (i = 0; i < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]].size(); i++){
                    if (v.find(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]][i].first) == v.end()){
                        v[sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]][i].first] = ratio * sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]][i].second;
                    }
                    else{
                        v.erase(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]][i].first);
                    }
                }


                
                st.push(sys->edgelink[sys->nodeEdge[st.top()][j].first * 2 + 1]);
                record = 1;
                break;
            }
        }
        if (record == 0){
            st.pop();
        }
    }



    for (auto vi : v){
        if (abs(vi.second) > 1e-5){
            RowId[num] = vi.first;
            ColId[num] = leng;
            Val[num] = vi.second;
            num++;
        }
    }
    leng++;


    free(visited);
    visited = NULL;
    
    free(nodeset);
    nodeset = NULL;

    return 0;

}


//int solveV0dSystem(fdtdMesh *sys, double *dRhs, double *y0d, int leng_v0d1){
//    
//    clock_t t1 = clock();
//    /* A\b1 */
//    double *d = &(sys->Adval[0]);
//    int *id = &(sys->AdRowId1[0]);
//    int *jd = &(sys->AdColId[0]);
//    
//    void *ptd[64];
//    int mtyped;
//    int iparmd[64];
//    double dparmd[64];
//    int maxfctd, mnumd, phased, errord, solverd;
//    int num_processd;   //number of processors
//    int v0csin;
//    int permd;
//    int nrhs = 1;
//    int msglvld = 0;    //print statistical information
//
//    mtyped = 11;    // real and not symmetric
//    solverd = 0;
//    errord = 0;
//    maxfctd = 1;    //maximum number of numerical factorizations
//    mnumd = 1;    //which factorization to use
//    phased = 13;    //analysis
//
//    pardisoinit(ptd, &mtyped, iparmd);
//    nrhs = 1;
//    iparmd[38] = 1;
//    iparmd[34] = 1;    //0-based indexing
//    pardiso(ptd, &maxfctd, &mnumd, &mtyped, &phased, &leng_v0d1, d, id, jd, &permd, &nrhs, iparmd, &msglvld, dRhs, y0d, &errord);
//    
//    cout << "Time to this point: " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
//
//    return 0;
//
//}

int merge_v0d1(fdtdMesh *sys, double block1_x, double block1_y, double block2_x, double block2_y, double block3_x, double block3_y, myint &v0d1num, myint &leng_v0d1, myint &v0d1anum, myint &leng_v0d1a, myint *map, double sideLen){
    int *visited;
    clock_t t1;
    double t, ta;
    double ratio;
    double startx, starty;    // the start coordinates of each block
    vector<int> st;    // dfs stack
    //vector<int> ind;
    int indsize, i;
    int indx, indy;
    int mark;
    int status;
    int indnum;
    int *markLayerNode = (int*)calloc(sys->N_node_s, sizeof(int));
    //int *markProSide = (int*)calloc(sys->N_node_s, sizeof(int));


    for (int i = 0; i < sys->numPorts; i++){
        for (int j = 0; j < sys->cdtNumNode[sys->portCoor[i].portCnd-1]; j++){
            markLayerNode[sys->conductor[sys->portCoor[i].portCnd - 1].node[j] % (sys->N_node_s)] = 1;
        }
    }
    
    /*for (int i = 0; i < sys->N_node_s; i++){
        indx = i / (sys->N_cell_y + 1);
        indy = i % (sys->N_cell_y + 1);
        if (sys->xn[indx] >= 2256.2e-6 && sys->xn[indx] <= 2814.4e-6){
            if (sys->yn[indy] >= 3079.8e-6 && sys->yn[indy] <= 3861.3e-6){
                markLayerNode[i] = 1;
            }
        }
    }*/

    
    //t1 = clock();
    //for (auto ni : layerNode){
    //    indx = ni / (sys->N_cell_y + 1);
    //    indy = ni % (sys->N_cell_y + 1);
    //    if ((indx - 1) * (sys->N_cell_y + 1) + indy >= 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (sys->xn[indx] - sys->xn[indx - 1]) < sideLen){    // this node is on the boundary of the projection
    //        status = setsideLen((indx) * (sys->N_cell_y + 1) + indy, sideLen, markLayerNode, markProSide, sys);
    //        continue;
    //    }
    //    else if ((indx + 1)* (sys->N_cell_y + 1) + indy < sys->N_node_s && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (sys->xn[indx + 1] - sys->xn[indx]) < sideLen){
    //        status = setsideLen((indx) * (sys->N_cell_y + 1) + indy, sideLen, markLayerNode, markProSide, sys);
    //        continue;
    //    }
    //    else if (indx * (sys->N_cell_y + 1) + indy - 1 >= 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy - 1] == 0 && (sys->yn[indy] - sys->yn[indy - 1]) < sideLen){
    //        status = setsideLen(indx * (sys->N_cell_y + 1) + indy, sideLen, markLayerNode, markProSide, sys);
    //        continue;
    //    }
    //    else if (indx * (sys->N_cell_y + 1) + indy + 1 < sys->N_node_s && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (sys->yn[indy + 1] - sys->yn[indy]) < sideLen){
    //        status = setsideLen(indx * (sys->N_cell_y + 1) + indy, sideLen, markLayerNode, markProSide, sys);
    //        continue;
    //    }
    //}
    //cout << "Time for finding the side nodes is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;


    leng_v0d1 = 0;
    leng_v0d1a = 0;
    v0d1num = 0;
    v0d1anum = 0;
    /* first asign a larger number of storage, don't need to calculate the entries twice */
    myint *v0d1RowId = (myint*)malloc(sys->N_edge * 2 * sizeof(myint));
    myint *v0d1ColId = (myint*)malloc(sys->N_edge * 2 * sizeof(myint));
    double *v0d1val = (double*)malloc(sys->N_edge * 2 * sizeof(double));
    myint *v0d1aRowId = (myint*)malloc(sys->N_edge * 2 * sizeof(myint));
    myint *v0d1aColId = (myint*)malloc(sys->N_edge * 2 * sizeof(myint));
    double *v0d1aval = (double*)malloc(sys->N_edge * 2 * sizeof(double));
    int count = 1;    /* count which box it is */
    clock_t t2 = clock();
    unordered_map<myint, double> va, v;
    t = 0; ta = 0;
    for (int iz = 1; iz < sys->nz - 1; iz++){    // merge on each layer, not in the conductor
        visited = (int*)calloc(sys->nx * sys->ny, sizeof(int));
        for (int ix = 0; ix < sys->nx; ix++){
            for (int iy = 0; iy < sys->ny; iy++){
                if (visited[ix * (sys->N_cell_y + 1) + iy] == 0 && sys->markNode[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] == 0){
                    va.clear();
                    v.clear();
                    if (markLayerNode[ix * (sys->N_cell_y + 1) + iy] == 0 && sys->markProSide[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] == 0){    // this point is not visited and it is outside the conductor, not in the projection of the excited conductor
                        //if (!ind.empty())
                        //    ind.clear();
                        startx = sys->xn[ix];
                        starty = sys->yn[iy];
                        
                        map[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] = count;
                        //ind.push_back(iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy);
                        st.push_back(ix * (sys->N_cell_y + 1) + iy);
                        visited[ix * (sys->N_cell_y + 1) + iy] = 1;
                        for (int i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            va[sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first] = sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                        }
                        for (int i = 0; i < sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            v[sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first] = sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                        }
                        while (!st.empty()){
                            mark = 0;
                            indx = (st.back()) / (sys->N_cell_y + 1);
                            indy = st.back() % (sys->N_cell_y + 1);
                            if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                                if (sys->markNode[iz * sys->N_node_s + st.back() + sys->N_cell_y + 1] == 0 && visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && sys->markProSide[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block1_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block1_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + indx * (sys->N_cell_y + 1) + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                        st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                        visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                        map[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] = count;

                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indx != 0){    // it must have a left x edge, thus left x node
                                if (sys->markNode[iz * sys->N_node_s + st.back() - sys->N_cell_y - 1] == 0 && visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && sys->markProSide[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block1_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block1_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (indx - 1) * (sys->N_cell_y + 1) + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                        st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                        visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                        map[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] = count;

                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                                if (sys->markNode[iz * sys->N_node_s + st.back() + 1] == 0 && visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && sys->markProSide[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block1_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block1_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                        st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                        visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                        map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = count;

                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indy != 0){    // it must have a closer y edge, thus closer y node
                                if (sys->markNode[iz * sys->N_node_s + st.back() - 1] == 0 && visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && sys->markProSide[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block1_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block1_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy - 1){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                        st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                        visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                        map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = count;

                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (mark == 0){
                                st.pop_back();
                            }
                        }
                        indnum = v.size();
                        
                        if (indnum != 0){
                            
                            for (auto vi : v){
                                v0d1RowId[v0d1num] = vi.first;
                                v0d1ColId[v0d1num] = leng_v0d1;
                                v0d1val[v0d1num] = vi.second;
                                v0d1num++;
                            }
                            leng_v0d1++;
                            
                            //status = nodeAddAvgLarger(&ind[0], indnum, sys->N_edge, sys, v0d1anum, leng_v0d1a, v0d1aRowId, v0d1aColId, v0d1aval);    // used to find v0d1anum and leng_v0d1a
                            //if (status != 0)
                            //    return status;
                            for (auto vai : va){
                                v0d1aRowId[v0d1anum] = vai.first;
                                v0d1aColId[v0d1anum] = leng_v0d1a;
                                v0d1aval[v0d1anum] = vai.second;
                                v0d1anum++;
                            }
                            leng_v0d1a++;
                            count++;
                        }  

                    }

                    else if (markLayerNode[ix * (sys->N_cell_y + 1) + iy] == 1 && sys->markProSide[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] == 0){//&& sys->exciteCdtLayer[iz] == 1){    // this point is not visited and it is outside the conductor, in the projection of the excited conductor
                        //while (!ind.empty())
                        //    ind.clear();
                        startx = sys->xn[ix];
                        starty = sys->yn[iy];
                        //ind.push_back(iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy);
                        map[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] = count;
                        st.push_back(ix * (sys->N_cell_y + 1) + iy);
                        visited[ix * (sys->N_cell_y + 1) + iy] = 1;
                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            va[sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first] = sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                        }
                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            v[sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first] = sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                        }
                        while (!st.empty()){
                            mark = 0;
                            indx = (st.back()) / (sys->N_cell_y + 1);
                            indy = st.back() % (sys->N_cell_y + 1);

                            if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                                if (sys->markNode[iz * sys->N_node_s + st.back() + sys->N_cell_y + 1] == 0 && visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 1 && sys->markProSide[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + indx * (sys->N_cell_y + 1) + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                        st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                        visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                        map[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indx != 0){    // it must have a left x edge, thus left x node
                                if (sys->markNode[iz * sys->N_node_s + st.back() - sys->N_cell_y - 1] == 0 && visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 1 && sys->markProSide[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (indx - 1) * (sys->N_cell_y + 1) + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                        st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                        visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                        map[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                                if (sys->markNode[iz * sys->N_node_s + st.back() + 1] == 0 && visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 1 && sys->markProSide[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block2_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                        st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                        visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                        map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indy != 0){    // it must have a closer y edge, thus closer y node
                                if (sys->markNode[iz * sys->N_node_s + st.back() - 1] == 0 && visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 1 && sys->markProSide[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block2_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy - 1){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                        st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                        visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                        map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (mark == 0){
                                st.pop_back();
                            }
                        }
                        indnum = v.size();
                        
                        if (indnum != 0){
                            
                            for (auto vi : v){
                                v0d1RowId[v0d1num] = vi.first;
                                v0d1ColId[v0d1num] = leng_v0d1;
                                v0d1val[v0d1num] = vi.second;
                                v0d1num++;
                            }
                            leng_v0d1++;

                            for (auto vai : va){
                                v0d1aRowId[v0d1anum] = vai.first;
                                v0d1aColId[v0d1anum] = leng_v0d1a;
                                v0d1aval[v0d1anum] = vai.second;
                                v0d1anum++;
                            }
                            leng_v0d1a++;
                            count++;
                        }
                    }
                    else{
                        /*for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            v0d1RowId[v0d1num] = sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first;
                            v0d1ColId[v0d1num] = leng_v0d1;
                            v0d1val[v0d1num] = sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                            v0d1num++;
                        }
                        map[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] = count;
                        count++;
                        leng_v0d1++;
                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            v0d1aRowId[v0d1anum] = sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first;
                            v0d1aColId[v0d1anum] = leng_v0d1a;
                            v0d1aval[v0d1anum] = sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                            v0d1anum++;
                        }
                        leng_v0d1a++;*/


                       
                        
                        
                        startx = sys->xn[ix];
                        starty = sys->yn[iy];

                        map[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy] = count;
                        st.push_back(ix * (sys->N_cell_y + 1) + iy);
                        visited[ix * (sys->N_cell_y + 1) + iy] = 1;
                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            va[sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first] = sys->nodeEdgea[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                        }
                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy].size(); i++){
                            v[sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].first] = sys->nodeEdge[iz * sys->N_node_s + ix * (sys->N_cell_y + 1) + iy][i].second;
                        }
                        while (!st.empty()){
                            mark = 0;
                            indx = (st.back()) / (sys->N_cell_y + 1);
                            indy = st.back() % (sys->N_cell_y + 1);

                            if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                                if (sys->markNode[iz * sys->N_node_s + st.back() + sys->N_cell_y + 1] == 0 && visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && sys->markProSide[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 1){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block3_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block3_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + indx * (sys->N_cell_y + 1) + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                        st.push_back((indx + 1)*(sys->N_cell_y + 1) + indy);
                                        visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                        map[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indx != 0){    // it must have a left x edge, thus left x node
                                if (sys->markNode[iz * sys->N_node_s + st.back() - sys->N_cell_y - 1] == 0 && visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && sys->markProSide[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 1){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block3_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block3_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + (sys->N_cell_x + 1) * sys->N_cell_y + (indx - 1) * (sys->N_cell_y + 1) + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                        st.push_back((indx - 1)*(sys->N_cell_y + 1) + indy);
                                        visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                        map[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                                if (sys->markNode[iz * sys->N_node_s + st.back() + 1] == 0 && visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && sys->markProSide[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 1){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block3_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block3_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                        st.push_back((indx)*(sys->N_cell_y + 1) + indy + 1);
                                        visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                        map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (indy != 0){    // it must have a closer y edge, thus closer y node
                                if (sys->markNode[iz * sys->N_node_s + st.back() - 1] == 0 && visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && sys->markProSide[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 1){    // this node is in dielectric and this node is not visited
                                    if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block3_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block3_y){    // this node is within the block area
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy - 1){
                                                ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first];
                                                break;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == va.end()){
                                                va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                            else{
                                                va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                                //va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] += ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                        }
                                        for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                            if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == v.end()){
                                                v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                            }
                                            else{
                                                v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                            }
                                        }
                                        //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                        st.push_back((indx)*(sys->N_cell_y + 1) + indy - 1);
                                        visited[(indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                        map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = count;
                                        mark = 1;
                                        continue;
                                    }
                                }
                            }
                            if (mark == 0){
                                st.pop_back();
                            }
                        }
                        indnum = v.size();

                        if (indnum != 0){

                            for (auto vi : v){
                                v0d1RowId[v0d1num] = vi.first;
                                v0d1ColId[v0d1num] = leng_v0d1;
                                v0d1val[v0d1num] = vi.second;
                                v0d1num++;
                            }
                            leng_v0d1++;

                            for (auto vai : va){
                                v0d1aRowId[v0d1anum] = vai.first;
                                v0d1aColId[v0d1anum] = leng_v0d1a;
                                v0d1aval[v0d1anum] = vai.second;
                                v0d1anum++;
                            }
                            leng_v0d1a++;
                            count++;
                        }
                    }
                }

            }
        }
        free(visited); visited = NULL;
    }

    

    /* V0d2 generation */
    int j;
    
    //int indz;
    //if (sys->numCdt > 0){

    //    int *ind2;
    //    
    //    for (i = 0; i < sys->numCdt; i++){
    //        if (sys->conductor[i].markPort != -1){
    //            ind2 = (int*)malloc(sys->cdtNumNode[i] * sizeof(int));
    //            indnum = 0;
    //            for (j = 0; j < sys->cdtNumNode[i]; j++){
    //                ind2[indnum] = (sys->conductor[i].node[j]);
    //                indnum++;
    //            }
    //            status = nodeAddLarger(ind2, indnum, sys->N_edge, sys, v0d1num, leng_v0d1, v0d1RowId, v0d1ColId, v0d1val);    // used to generate v0d2
    //            if (status != 0)
    //                return status;
    //            status = nodeAddAvgLarger(ind2, indnum, sys->N_edge, sys, v0d1anum, leng_v0d1a, v0d1aRowId, v0d1aColId, v0d1aval);    // used to generate v0d2a
    //            if (status != 0)
    //                return status;
    //            free(ind2); ind2 = NULL;
    //        }
    //    }
    //}
    visited = (int*)calloc(sys->N_node, sizeof(int));
    for (i = 0; i < sys->numCdt; i++){
        if (sys->conductor[i].markPort == -1){
            continue;
        }
        else{
            v.clear();
            va.clear();
            st.push_back(sys->conductor[i].node[0]);
            visited[sys->conductor[i].node[0]] = 1;
            map[sys->conductor[i].node[0]] = count;
            
            for (indx = 0; indx < sys->nodeEdge[st.back()].size(); indx++){
                v[sys->nodeEdge[st.back()][indx].first] = sys->nodeEdge[st.back()][indx].second;
            }
            for (indx = 0; indx < sys->nodeEdgea[st.back()].size(); indx++){
                va[sys->nodeEdgea[st.back()][indx].first] = sys->nodeEdgea[st.back()][indx].second;
            }
            while (!st.empty()){
                mark = 0;
                for (j = 0; j < sys->nodeEdge[st.back()].size(); j++){
                    if (sys->markEdge[sys->nodeEdge[st.back()][j].first] == i + 1){
                        if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] == 0)){
                            map[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] = count;
                            
                            for (indx = 0; indx < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]].size(); indx++){
                                if (sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first == sys->nodeEdge[st.back()][j].first){
                                    ratio = -1 / sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].second * va[sys->nodeEdge[st.back()][j].first];
                                    break;
                                }
                            }
                            for (indx = 0; indx < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]].size(); indx++){
                                if (va.find(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first) == va.end()){
                                    va[sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first] = ratio * sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].second;
                                }
                                else{
                                    va.erase(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first);
                                }
                            }
                            for (indx = 0; indx < sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]].size(); indx++){
                                if (v.find(sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first) == v.end()){
                                    v[sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first] = sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].second;
                                }
                                else{
                                    v.erase(sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]][indx].first);
                                }
                            }
                            
                            visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]] = 1;
                            st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2]);
                            mark = 1;
                            break;
                        }
                        else if ((sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1] != st.back() && visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] == 0)){
                            map[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] = count;
                           
                            for (indx = 0; indx < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]].size(); indx++){
                                if (sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first == sys->nodeEdge[st.back()][j].first){
                                    ratio = -1 / sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].second * va[sys->nodeEdge[st.back()][j].first];
                                    break;
                                }
                            }
                            for (indx = 0; indx < sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]].size(); indx++){
                                if (va.find(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first) == va.end()){
                                    va[sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first] = ratio * sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].second;
                                }
                                else{
                                    va.erase(sys->nodeEdgea[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first);
                                }
                            }
                            for (indx = 0; indx < sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]].size(); indx++){
                                if (v.find(sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first) == v.end()){
                                    v[sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first] = sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].second;
                                }
                                else{
                                    v.erase(sys->nodeEdge[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]][indx].first);
                                }
                            }
                            visited[sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]] = 1;
                            st.push_back(sys->edgelink[sys->nodeEdge[st.back()][j].first * 2 + 1]);
                            mark = 1;
                            break;
                        }
                    }
                }
                if (mark == 0){
                    st.pop_back();
                }
            }
            indnum = v.size();
            if (indnum != 0){
                for (auto vi : v){
                    v0d1RowId[v0d1num] = vi.first;
                    v0d1ColId[v0d1num] = leng_v0d1;
                    v0d1val[v0d1num] = vi.second;
                    v0d1num++;
                }
                leng_v0d1++;
                for (auto vai : va){
                    v0d1aRowId[v0d1anum] = vai.first;
                    v0d1aColId[v0d1anum] = leng_v0d1a;
                    v0d1aval[v0d1anum] = vai.second;
                    v0d1anum++;
                }
                leng_v0d1a++;
                count++;
            }
        }
    }
    

    
    
    sys->v0d1RowId = (myint*)malloc(v0d1num * sizeof(myint));
    sys->v0d1ColId = (myint*)malloc(v0d1num * sizeof(myint));
    sys->v0d1val = (double*)malloc(v0d1num * sizeof(double));
    sys->v0d1aRowId = (myint*)malloc(v0d1anum * sizeof(myint));
    sys->v0d1aColId = (myint*)malloc(v0d1anum * sizeof(myint));
    sys->v0d1aval = (double*)malloc(v0d1anum * sizeof(double));
    

    for (myint i = 0; i < v0d1num; i++){
        sys->v0d1RowId[i] = v0d1RowId[i];
        sys->v0d1ColId[i] = v0d1ColId[i];
        sys->v0d1val[i] = v0d1val[i];
    }
    for (myint i = 0; i < v0d1anum; i++){
        sys->v0d1aRowId[i] = v0d1aRowId[i];
        sys->v0d1aColId[i] = v0d1aColId[i];
        sys->v0d1aval[i] = v0d1aval[i];
    }
    

    free(v0d1RowId); v0d1RowId = NULL;
    free(v0d1ColId); v0d1ColId = NULL;
    free(v0d1val); v0d1val = NULL;
    free(v0d1aRowId); v0d1aRowId = NULL;
    free(v0d1aColId); v0d1aColId = NULL;
    free(v0d1aval); v0d1aval = NULL;

    

    /*ofstream outfile;
    outfile.open("color.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_node; i++){
        outfile << coloring[i] << endl;
    }
    outfile.close();*/

    return 1;
}

int setsideLen(int node, double sideLen, int *markLayerNode, int *markProSide, fdtdMesh *sys){
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
    while (!q.empty()){
        indx = q.front() / (sys->N_cell_y + 1);
        indy = q.front() % (sys->N_cell_y + 1);

        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node

            if (visited[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx + 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx + 1] - startx), 2) + pow((sys->yn[indy] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx + 1)*(sys->N_cell_y + 1) + indy);
                    visited[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                    markProSide[(indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                }
            }
        }

        if (indx != 0){    // it must have a left x edge, thus left x node
            if (visited[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && markLayerNode[(indx - 1) * (sys->N_cell_y + 1) + indy] == 0){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx - 1] - startx), 2) + pow((sys->yn[indy] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx - 1)*(sys->N_cell_y + 1) + indy);
                    visited[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                    markProSide[(indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                }
            }
        }

        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
            if (visited[indx * (sys->N_cell_y + 1) + indy + 1] == 0 && markLayerNode[indx * (sys->N_cell_y + 1) + indy + 1] == 0){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx] - startx), 2) + pow((sys->yn[indy + 1] - starty), 2)) <= sideLen){    // this node is within the block area
                    q.push((indx)*(sys->N_cell_y + 1) + indy + 1);
                    visited[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                    markProSide[(indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                }
            }
        }

        if (indy != 0){    // it must have a closer y edge, thus closer y node
            if (visited[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && markLayerNode[(indx)* (sys->N_cell_y + 1) + indy - 1] == 0){    // this node is in dielectric and this node is not visited
                if (sqrt(pow((sys->xn[indx] - startx), 2) + pow((sys->yn[indy - 1] - starty), 2)) <= sideLen){    // this node is within the block area
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

int merge_v0c(fdtdMesh *sys, double block_x, double block_y, double block2_x, double block2_y, myint &v0cnum, myint &leng_v0c, myint &v0canum, myint &leng_v0ca, myint *map){
    int *visited;
    double ratio;
    double startx, starty;    // the start coordinates of each block
    queue<int> st;    // dfs stack
    //vector<int> ind;
    int indsize;
    int indx, indy;
    int mark;
    int status;
    int markcond;
    int count = 0;
    leng_v0c = 0;
    leng_v0ca = 0;
    v0cnum = 0;
    v0canum = 0;
    int indnum;
    int ix, iy, iz;
    int n;
    int i, j;
    int map_count = 1;
    
    myint *v0cRowId = (myint*)malloc(2 * sys->N_edge * sizeof(myint));
    myint *v0cColId = (myint*)malloc(2 * sys->N_edge * sizeof(myint));
    double *v0cval = (double*)malloc(2 * sys->N_edge * sizeof(double));
    myint *v0caRowId = (myint*)malloc(2 * sys->N_edge * sizeof(myint));
    myint *v0caColId = (myint*)malloc(2 * sys->N_edge * sizeof(myint));
    double *v0caval = (double*)malloc(2 * sys->N_edge * sizeof(double));
    unordered_map<myint, double> v, va;
    visited = (int*)calloc(sys->N_node, sizeof(int));


    for (int ic = 0; ic < sys->numCdt; ic++){
        
        if (sys->conductor[ic].markPort <= 0){    // not excited conductors
            markcond = ic + 1;
            //visited = (int*)calloc(sys->N_node, sizeof(int));
            n = sys->cdtNumNode[ic] - 1;
            if (sys->conductor[ic].markPort == -1)
                n = sys->cdtNumNode[ic];
            for (int jc = 0; jc < n; jc++){
                if (visited[sys->conductor[ic].node[jc]] == 0 && sys->conductor[ic].node[jc] >= sys->N_node_s && sys->conductor[ic].node[jc] < sys->N_node - sys->N_node_s){
                    v.clear();
                    va.clear();
                    //if (!ind.empty())
                    //    ind.clear();
                    iz = sys->conductor[ic].node[jc] / (sys->N_node_s);
                    ix = (sys->conductor[ic].node[jc] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iy = sys->conductor[ic].node[jc] % (sys->N_cell_y + 1);
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    //ind.push_back(sys->conductor[ic].node[jc]);
                    st.push(ix * (sys->N_cell_y + 1) + iy);
                    visited[sys->conductor[ic].node[jc]] = 1;
                    map[sys->conductor[ic].node[jc]] = map_count;
                    
                    for (i = 0; i < sys->nodeEdge[st.front() + iz * sys->N_node_s].size(); i++){
                        v[sys->nodeEdge[st.front() + iz * sys->N_node_s][i].first] = sys->nodeEdge[st.front() + iz * sys->N_node_s][i].second;
                    }
                    for (i = 0; i < sys->nodeEdgea[st.front() + iz * sys->N_node_s].size(); i++){
                        va[sys->nodeEdgea[st.front() + iz * sys->N_node_s][i].first] = sys->nodeEdgea[st.front() + iz * sys->N_node_s][i].second;
                    }
                    while (!st.empty()){

                        mark = 0;
                        indx = (st.front()) / (sys->N_cell_y + 1);
                        indy = st.front() % (sys->N_cell_y + 1);

                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * indx + indy] == markcond && visited[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + indx * (sys->N_cell_y + 1) + indy){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first];
                                            break;
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    map[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1) * (sys->N_cell_y + 1) + indy){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first];
                                            break;
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    map[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy] == markcond && visited[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first];
                                            break;
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy - 1] == markcond && visited[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy - 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy - 1){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first];
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        //if (mark == 0){
                            st.pop();
                        //}
                    }

                    indnum = v.size();


                    //if (indnum != 0){
                    //    status = nodeAddLarger(&ind[0], indnum, sys->N_edge, sys, v0cnum, leng_v0c, v0cRowId, v0cColId, v0cval);    // used to find v0d1num and leng_v0d1
                    //    if (status != 0)
                    //        return status;

                    //    status = nodeAddAvgLarger(&ind[0], indnum, sys->N_edge, sys, v0canum, leng_v0ca, v0caRowId, v0caColId, v0caval);    // used to find v0d1anum and leng_v0d1a
                    //    if (status != 0)
                    //        return status;
                    //}
                    if (indnum != 0){
                        for (auto vi : v){
                            v0cRowId[v0cnum] = vi.first;
                            v0cColId[v0cnum] = leng_v0c;
                            v0cval[v0cnum] = vi.second;
                            v0cnum++;
                        }
                        leng_v0c++;
                        for (auto vai : va){
                            v0caRowId[v0canum] = vai.first;
                            v0caColId[v0canum] = leng_v0ca;
                            v0caval[v0canum] = vai.second;
                            v0canum++;
                        }
                        leng_v0ca++;
                        map_count++;
                    }

                }
            }

            sys->acu_cnno[count + 1] = leng_v0c;
            count++;
            //free(visited); visited = NULL;
        }
        else{


            markcond = ic + 1;
            
            //visited = (int*)calloc(sys->N_node, sizeof(int));
            n = sys->cdtNumNode[ic] - 1;
            for (int jc = 0; jc < n; jc++){
                if (visited[sys->conductor[ic].node[jc]] == 0 && sys->conductor[ic].node[jc] >= sys->N_node_s && sys->conductor[ic].node[jc] < sys->N_node - sys->N_node_s){
                    v.clear();
                    va.clear();
                    //if (!ind.empty())
                    //    ind.clear();
                    iz = sys->conductor[ic].node[jc] / (sys->N_node_s);
                    ix = (sys->conductor[ic].node[jc] - iz * sys->N_node_s) / (sys->N_cell_y + 1);
                    iy = sys->conductor[ic].node[jc] % (sys->N_cell_y + 1);
                    startx = sys->xn[ix];
                    starty = sys->yn[iy];
                    //ind.push_back(sys->conductor[ic].node[jc]);
                    st.push(ix * (sys->N_cell_y + 1) + iy);
                    visited[sys->conductor[ic].node[jc]] = 1;
                    map[sys->conductor[ic].node[jc]] = map_count;
                    
                    for (i = 0; i < sys->nodeEdge[st.front() + iz * sys->N_node_s].size(); i++){
                        v[sys->nodeEdge[st.front() + iz * sys->N_node_s][i].first] = sys->nodeEdge[st.front() + iz * sys->N_node_s][i].second;
                    }
                    for (i = 0; i < sys->nodeEdgea[st.front() + iz * sys->N_node_s].size(); i++){
                        va[sys->nodeEdgea[st.front() + iz * sys->N_node_s][i].first] = sys->nodeEdgea[st.front() + iz * sys->N_node_s][i].second;
                    }
                    while (!st.empty()){

                        mark = 0;
                        indx = (st.front()) / (sys->N_cell_y + 1);
                        indy = st.front() % (sys->N_cell_y + 1);

                        if (indx != sys->nx - 1){    // it must have a right x edge, thus right x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * indx + indy] == markcond && visited[iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx + 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx + 1] - startx) >= 0 && (sys->xn[indx + 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + indx * (sys->N_cell_y + 1) + indy){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first];
                                            break;
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy);
                                    st.push((indx + 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    map[iz * sys->N_node_s + (indx + 1)*(sys->N_cell_y + 1) + indy] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        if (indx != 0){    // it must have a left x edge, thus left x node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_x + 1) * (indx - 1) + indy] == markcond && visited[iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy] == 0 && (iz * sys->N_node_s + (indx - 1) * (sys->N_cell_y + 1) + indy != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx - 1] - startx) >= 0 && (sys->xn[indx - 1] - startx) <= block2_x && (sys->yn[indy] - starty) >= 0 && (sys->yn[indy] - starty) <= block2_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1) * (sys->N_cell_y + 1) + indy){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first];
                                            break;
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy);
                                    st.push((indx - 1)*(sys->N_cell_y + 1) + indy);
                                    visited[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = 1;
                                    map[iz * sys->N_node_s + (indx - 1)*(sys->N_cell_y + 1) + indy] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        if (indy != sys->ny - 1){    // it must have a farther y edge, thus farther y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy] == markcond && visited[iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy + 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy + 1] - starty) >= 0 && (sys->yn[indy + 1] - starty) <= block2_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first];
                                            break;
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1);
                                    st.push((indx)*(sys->N_cell_y + 1) + indy + 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = 1;
                                    map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy + 1] = map_count;
                                    //mark = 1;
                                    
                                    //continue;
                                }
                            }
                        }
                        if (indy != 0){    // it must have a closer y edge, thus closer y node
                            if (sys->markEdge[iz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * indx + indy - 1] == markcond && visited[iz * sys->N_node_s + (indx)* (sys->N_cell_y + 1) + indy - 1] == 0 && (iz * sys->N_node_s + indx * (sys->N_cell_y + 1) + indy - 1 != sys->conductor[ic].node[sys->cdtNumNode[ic] - 1] || sys->conductor[ic].markPort == -1)){    // this node is in conductor and this node is not visited
                                if ((sys->xn[indx] - startx) >= 0 && (sys->xn[indx] - startx) <= block2_x && (sys->yn[indy - 1] - starty) >= 0 && (sys->yn[indy - 1] - starty) <= block2_y){    // this node is within the block area
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                        if (sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first == iz * (sys->N_edge_s + sys->N_edge_v) + indx * sys->N_cell_y + indy - 1){
                                            ratio = -1 / sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second * va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first];
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                        if (va.find(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == va.end()){
                                            va[sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = ratio * sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                        }
                                        else{
                                            va.erase(sys->nodeEdgea[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                        }
                                    }
                                    for (i = 0; i < sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1].size(); i++){
                                        if (v.find(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first) == v.end()){
                                            v[sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first] = sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].second;
                                        }
                                        else{
                                            v.erase(sys->nodeEdge[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1][i].first);
                                        }
                                    }
                                    //ind.push_back(iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1);
                                    st.push((indx)*(sys->N_cell_y + 1) + indy - 1);
                                    visited[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = 1;
                                    map[iz * sys->N_node_s + (indx)*(sys->N_cell_y + 1) + indy - 1] = map_count;
                                    //mark = 1;

                                    //continue;
                                }
                            }
                        }
                        //if (mark == 0){

                            st.pop();

                        //}
                    }

                    indnum = v.size();


                    //if (indnum != 0){

                    //    status = nodeAddLarger(&ind[0], indnum, sys->N_edge, sys, v0cnum, leng_v0c, v0cRowId, v0cColId, v0cval);    // used to find v0d1num and leng_v0d1
                    //    if (status != 0)
                    //        return status;

                    //    status = nodeAddAvgLarger(&ind[0], indnum, sys->N_edge, sys, v0canum, leng_v0ca, v0caRowId, v0caColId, v0caval);    // used to find v0d1anum and leng_v0d1a
                    //    if (status != 0)
                    //        return status;
                    //}
                    if (indnum != 0){
                        for (auto vi : v){
                            v0cRowId[v0cnum] = vi.first;
                            v0cColId[v0cnum] = leng_v0c;
                            v0cval[v0cnum] = vi.second;
                            v0cnum++;
                        }
                        leng_v0c++;
                        for (auto vai : va){
                            v0caRowId[v0canum] = vai.first;
                            v0caColId[v0canum] = leng_v0ca;
                            v0caval[v0canum] = vai.second;
                            v0canum++;
                        }
                        leng_v0ca++;
                        map_count++;
                    }
                }
            }

            sys->acu_cnno[count + 1] = leng_v0c;
            count++;
            //free(visited); visited = NULL;

        }
    }

    cout << endl;
    sys->v0cRowId = (myint*)malloc(v0cnum * sizeof(myint));
    sys->v0cColId = (myint*)malloc(v0cnum * sizeof(myint));
    sys->v0cval = (double*)malloc(v0cnum * sizeof(double));
    sys->v0caRowId = (myint*)malloc(v0canum * sizeof(myint));
    sys->v0caColId = (myint*)malloc(v0canum * sizeof(myint));
    sys->v0caval = (double*)malloc(v0canum * sizeof(double));
    
    for (int i = 0; i < v0cnum; i++){
        sys->v0cRowId[i] = v0cRowId[i];
        sys->v0cColId[i] = v0cColId[i];
        sys->v0cval[i] = v0cval[i];
    }
    for (int i = 0; i < v0canum; i++){
        sys->v0caRowId[i] = v0caRowId[i];
        sys->v0caColId[i] = v0caColId[i];
        sys->v0caval[i] = v0caval[i];
    }
    
    free(v0cRowId); v0cRowId = NULL;
    free(v0cColId); v0cColId = NULL;
    free(v0cval); v0cval = NULL;
    free(v0caRowId); v0caRowId = NULL;
    free(v0caColId); v0caColId = NULL;
    free(v0caval); v0caval = NULL;
    free(visited); visited = NULL;
    //ind.clear();

    return 1;
}
