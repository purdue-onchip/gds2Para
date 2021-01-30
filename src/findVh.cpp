#include "fdtd.hpp"

static bool comp(pair<complex<double>, myint> a, pair<complex<double>, myint> b)
{
    return (pow(a.first.real(), 2) + pow(a.first.imag(), 2)) < (pow(b.first.real(), 2) + pow(b.first.imag(), 2));
};

int find_Vh(fdtdMesh *sys, lapack_complex_double *u0, lapack_complex_double *u0a, int sourcePort){
    // upper and lower plane is PEC
    int t0n = 3;
    double dt = DT;
    double tau = 1.e-11;
    double t0 = t0n * tau;
    myint nt = 2 * t0 / dt;
    myint SG = 1000;
    double zero_value = 1e10;    // smaller than this value is discarded as null-space eigenvector
    double eps1 = 1.e-7;    // weight
    double eps2 = 1.e-2;    // wavelength error tolerance
    double *rsc = (double*)calloc((sys->N_edge - sys->bden), sizeof(double));
    double *xr = (double*)calloc((sys->N_edge - sys->bden) * 3, sizeof(double));
    myint start;
    myint index;
    int l = 0;
    
    vector<double> U0_base(sys->N_edge - sys->bden, 0);
    vector<vector<double>> U0;   // Since I need to put in new vectors into U0, I need dynamic storage
    double *U0c, *m0;
    vector<vector<double>> temp, temp1, temp2;
    double nn, nnr, nni;
    double *Cr, *D_sigr, *D_epsr;
    vector<vector<double>> Cr_p, D_sigr_p, D_epsr_p;
    vector<double> Cr_base;
    double *Ar, *Br, *BBr;
    vector<pair<complex<double>, myint>> dp, dp1;
    double *vl, *vr, *vl1;
    double *vr1;
    lapack_int info;
    lapack_int n, m;
    char jobvl, jobvr;
    lapack_int lda, ldb;
    lapack_int ldvl, ldvr;
    double *alphar;
    double *alphai;
    double *beta;
    double *tmp, *tmp1, *tmp2;
    int i_re1, i_nre1;
    lapack_complex_double *V_re, *V_nre, *V_re1, *V_nre1;
    V_re = (lapack_complex_double*)calloc(1, sizeof(lapack_complex_double));
    V_nre = (lapack_complex_double*)calloc(1, sizeof(lapack_complex_double));
    V_re1 = (lapack_complex_double*)calloc(1, sizeof(lapack_complex_double));
    V_nre1 = (lapack_complex_double*)calloc(1, sizeof(lapack_complex_double));
    int ind_i, ind_j;
    int i_re, i_nre;
    vector<int> re, nre;    // used to store the No. of eigenvectors in re and nre
    char job, side;
    lapack_int ilo, ihi;
    double *lscale, *rscale;
    double nor;
    lapack_complex_double *tmp3, *tmp4;
    int status;
    lapack_complex_double *y_re;
    lapack_complex_double *y_nre, *y_nre1;
    double x_re, x_nre;
    double maxvalue = 1.e+100;
    lapack_complex_double *tau_qr;
    int U0_i = 0;
    double zero_norm = 1.e-7;
    lapack_complex_double *m_new;
    lapack_int *ipiv;
    lapack_int iter;
    double *V_nre_norm, *tau1;
    lapack_complex_double *y_new, *V_new;
    lapack_complex_double *y;
    int *select;
    double *q;
    lapack_complex_double *qc;
    ofstream outfile;
    lapack_complex_double *m_re;
    double scale = 1.e+12;
    lapack_complex_double *RR;
    double eps = 1e-5;
    for (int ind = 1; ind <= nt; ind++){    // time marching to find the repeated eigenvectors
        for (int sourcePortSide = 0; sourcePortSide < sys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
            for (int inde = 0; inde < sys->portCoor[sourcePort].portEdge[sourcePortSide].size(); inde++){
                rsc[sys->mapEdge[sys->portCoor[sourcePort].portEdge[sourcePortSide][inde]]] = 2000 * exp(-pow((((dt * ind) - t0) / tau), 2)) + 2000 * (dt * ind - t0) * exp(-pow(((dt * ind - t0) / tau), 2)) * (-2 * (dt * ind - t0) / pow(tau, 2));

            }
        }

        // central difference
        index = 0;
        nn = 0;
        while (index < sys->leng_S){
            start = sys->SRowId[index];
            while (index < sys->leng_S && sys->SRowId[index] == start){
                xr[2 * (sys->N_edge - sys->bden) + sys->SRowId[index]] += sys->Sval[index] * xr[1 * (sys->N_edge - sys->bden) + sys->SColId[index]] * (-2) * pow(dt, 2);
                index++;
            }
            if (sys->markEdge[sys->mapEdgeR[start]] != 0){//SIGMA next line
                xr[2 * (sys->N_edge - sys->bden) + start] += -rsc[start] * 2 * pow(dt, 2) + dt * sys->stackSign[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * xr[start] - 2 * sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0 * xr[start] + 4 * sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0 * xr[1 * (sys->N_edge - sys->bden) + start];
                xr[2 * (sys->N_edge - sys->bden) + start] /= (2 * sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0 + dt * sys->stackSign[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)]);
            }
            else {
                xr[2 * (sys->N_edge - sys->bden) + start] += -rsc[start] * 2 * pow(dt, 2) - 2 * sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0 * xr[start] + 4 * sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0 * xr[1 * (sys->N_edge - sys->bden) + start];
                xr[2 * (sys->N_edge - sys->bden) + start] /= (2 * sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0);
            }
            nn += xr[1 * (sys->N_edge - sys->bden) + start] * xr[1 * (sys->N_edge - sys->bden) + start];
        }
        nn = sqrt(nn);
        //cout << "Step " << ind << "'s norm is " << nn << endl;
        if (ind == SG){
            l++;
            U0.push_back(U0_base);
            temp.push_back(U0_base);
            temp1.push_back(U0_base);
            temp2.push_back(U0_base);
            Cr_p.push_back(Cr_base);
            D_sigr_p.push_back(Cr_base);
            D_epsr_p.push_back(Cr_base);
            for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++){
                U0[U0_i][inde] = xr[1 * (sys->N_edge - sys->bden) + inde] / nn;
            }
            U0_i++;
            Cr = (double*)calloc(l * l, sizeof(double));
            D_sigr = (double*)calloc(l * l, sizeof(double));
            D_epsr = (double*)calloc(l * l, sizeof(double));
            index = 0;
            
            while (index < sys->leng_S){
                start = sys->SRowId[index];
                while (index < sys->leng_S && sys->SRowId[index] == start){
                    temp[0][sys->SRowId[index]] += sys->Sval[index] * U0[U0_i - 1][sys->SColId[index]];
                    index++;
                }
                if (sys->markEdge[sys->mapEdgeR[start]] != 0) {
                    temp1[0][start] = sqrt(sys->stackSign[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)]) * U0[U0_i - 1][start];//SIGMA
                    temp2[0][start] = sqrt(sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0) * U0[U0_i - 1][start];
                }
                else {
                    temp1[0][start] = 0;
                    temp2[0][start] = sqrt(sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0) * U0[U0_i - 1][start];
                }
            }
            for (myint inde = 0; inde < l; inde++){
                Cr_p.push_back(Cr_base);
                D_sigr_p.push_back(Cr_base);
                D_epsr_p.push_back(Cr_base);
                for (myint inde2 = 0; inde2 < l; inde2++){
                    Cr_p[inde].push_back(0);
                    D_sigr_p[inde].push_back(0);
                    D_epsr_p[inde].push_back(0);
                    for (myint inde3 = 0; inde3 < sys->N_edge - 2 * sys->N_edge_s; inde3++){
                        Cr[inde * l + inde2] += U0[inde][inde3] * temp[inde2][inde3];
                        D_sigr[inde * l + inde2] += temp1[inde][inde3] * temp1[inde2][inde3];
                        D_epsr[inde * l + inde2] += temp2[inde][inde3] * temp2[inde2][inde3];
                        Cr_p[inde][inde2] += U0[inde][inde3] * temp[inde2][inde3];
                        D_sigr_p[inde][inde2] += temp1[inde][inde3] * temp1[inde2][inde3];
                        D_epsr_p[inde][inde2] += temp2[inde][inde3] * temp2[inde2][inde3];
                    }
                }
            }
            Ar = (double*)malloc(4 * l * l * sizeof(double));
            Br = (double*)malloc(4 * l * l * sizeof(double));
            for (myint inde = 0; inde < l; inde++){
                for (myint inde2 = 0; inde2 < l; inde2++){
                    Ar[inde * 2 * l + inde2] = -Cr[inde * l + inde2];
                    Br[inde * 2 * l + inde2] = D_sigr[inde * l + inde2] * scale;
                }
            }
            for (myint inde = 0; inde < l; inde++){
                for (myint inde2 = l; inde2 < 2 * l; inde2++){
                    Ar[inde * 2 * l + inde2] = 0;
                    Br[inde * 2 * l + inde2] = D_epsr[inde * l + inde2 - l] * pow(scale, 2);
                }
            }
            for (myint inde = l; inde < 2 * l; inde++){
                for (myint inde2 = 0; inde2 < l; inde2++){
                    Ar[inde * 2 * l + inde2] = 0;
                    Br[inde * 2 * l + inde2] = D_epsr[(inde - l) * l + inde2] * pow(scale, 2);
                }
            }
            for (myint inde = l; inde < 2 * l; inde++){
                for (myint inde2 = l; inde2 < 2 * l; inde2++){
                    Ar[inde * 2 * l + inde2] = D_epsr[(inde - l) * l + inde2 - l] * pow(scale, 2);
                    Br[inde * 2 * l + inde2] = 0;
                }
            }

            n = 2 * l;
            jobvl = 'N'; jobvr = 'V';
            lda = 2 * l; ldb = 2 * l;
            ldvl = 2 * l; ldvr = 2 * l;
            alphar = (double*)malloc(2 * l * sizeof(double));
            alphai = (double*)malloc(2 * l * sizeof(double));
            beta = (double*)malloc(2 * l * sizeof(double));
            vl = (double*)malloc(4 * l * l * sizeof(double));
            vr = (double*)malloc(4 * l * l * sizeof(double));
            job = 'B';    // scale only
            lscale = (double*)malloc(n * sizeof(double));
            rscale = (double*)malloc(n * sizeof(double));
            
            info = LAPACKE_dggev(LAPACK_COL_MAJOR, jobvl, jobvr, n, Ar, lda, Br, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr);
            for (myint inde = 0; inde < 2 * l; inde++){
                //if (alphai[inde] == 0){
                //    dp.push_back(make_pair(alphar[inde] / beta[inde] * scale, inde));
                //}
                //else{
                if (abs(alphai[inde] / beta[inde]) > eps) {   // only those imaginary part larger than 0 eigenvalues are considered
                    dp.push_back(make_pair(alphar[inde] / beta[inde] * scale + 1i * alphai[inde] / beta[inde] * scale, inde));
                    inde++;
                    dp.push_back(make_pair(alphar[inde] / beta[inde] * scale + 1i * alphai[inde] / beta[inde] * scale, inde));
                }
            }
            
            sort(dp.begin(), dp.end(), comp);
            
            free(Cr); Cr = NULL;
            free(D_epsr); D_epsr = NULL;
            free(D_sigr); D_sigr = NULL;
            free(Ar); Ar = NULL;
            free(Br); Br = NULL;
            free(alphar); alphar = NULL;
            free(alphai); alphai = NULL;
            free(beta); beta = NULL;
        }
        else if (ind % SG == 0){
            tmp = (double*)malloc((sys->N_edge - sys->bden) * sizeof(double));
            tmp1 = (double*)calloc((sys->N_edge - sys->bden), sizeof(double));
            tmp2 = (double*)calloc(U0_i, sizeof(double));
            nn = 0;
            nor = 0;
            for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                tmp[inde] = xr[1 * (sys->N_edge - sys->bden) + inde];
                nor += xr[1 * (sys->N_edge - sys->bden) + inde] * xr[1 * (sys->N_edge - sys->bden) + inde];
            }
            nor = sqrt(nor);

            

            /* modified Gran Schmidt */
            clock_t t1 = clock();
            /*for (myint inde = 0; inde < U0_i; inde++){
                nn = 0;
                for (myint inde1 = 0; inde1 < sys->N_edge - 2 * sys->N_edge_s; inde1++){
                    nn += U0[inde][inde1] * U0[inde][inde1];
                }
                nn = sqrt(nn);
                for (myint inde1 = 0; inde1 < sys->N_edge - 2 * sys->N_edge_s; inde1++){
                    U0[inde][inde1] = U0[inde][inde1] / nn;
                }
                for (myint inde1 = inde + 1; inde1 < U0_i; inde1++){
                    nn = 0;
                    for (myint inde2 = 0; inde2 < sys->N_edge - 2 * sys->N_edge_s; inde2++){
                        nn += U0[inde][inde2] * U0[inde1][inde2];
                    }
                    q = (double*)malloc((sys->N_edge - 2 * sys->N_edge_s) * sizeof(double));
                    for (myint inde2 = 0; inde2 < sys->N_edge - 2 * sys->N_edge_s; inde2++){
                        q[inde2] = U0[inde][inde2] * nn;
                        U0[inde1][inde2] -= q[inde2];
                    }
                    free(q); q = NULL;
                }
                nn = 0;
                for (myint inde1 = 0; inde1 < sys->N_edge - 2 * sys->N_edge_s; inde1++){
                    nn += U0[inde][inde1] * tmp[inde1];
                }
                q = (double*)malloc((sys->N_edge - 2 * sys->N_edge_s) * sizeof(double));
                for (myint inde1 = 0; inde1 < sys->N_edge - 2 * sys->N_edge_s; inde1++){
                    q[inde1] = U0[inde][inde1] * nn;
                    tmp[inde1] = tmp[inde1] - q[inde1];
                }
                free(q); q = NULL;
            }
            nn = 0;
            for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
                nn += tmp[inde] * tmp[inde];
            }
            nn = sqrt(nn);*/
            
            for (myint inde = 0; inde < U0_i; inde++){
                nn = 0;
                for (myint inde1 = 0; inde1 < sys->N_edge - sys->bden; inde1++){
                    nn += tmp[inde1] * U0[inde][inde1];
                }
                for (myint inde1 = 0; inde1 < sys->N_edge - sys->bden; inde1++){
                    tmp[inde1] -= nn * U0[inde][inde1];
                }
            }
            nn = 0;
            for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++){
                nn += tmp[inde] * tmp[inde];
            }
            nn = sqrt(nn);

            //cout << "Modified Gram Schmidt takes " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
            //cout << "U0_i is " << U0_i << endl;
            /*ofstream out;
            out.open("U0.txt", std::ofstream::out | std::ofstream::trunc);
            for (myint inde = 0; inde < (sys->N_edge - 2 * sys->N_edge_s); inde++){
                for (myint inde2 = 0; inde2 < U0_i; inde2++){
                    out << U0[inde2][inde] << " ";
                }
                out << endl;
            }
            out.close();*/
            
            //cout << "New vector after deduction " << nn / nor << endl;
            if (nn / nor > zero_norm){     // this new vector contains new information
                U0.push_back(U0_base);
                temp.push_back(U0_base);
                temp1.push_back(U0_base);
                temp2.push_back(U0_base);
                for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                    U0[U0_i][inde] = tmp[inde] / nn;
                }
                U0_i++;
            }
            else{     // no new information go to the next loop
                for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++){
                    xr[inde] = xr[1 * (sys->N_edge - sys->bden) + inde];
                    xr[1 * (sys->N_edge - sys->bden) + inde] = xr[2 * (sys->N_edge - sys->bden) + inde];
                    xr[2 * (sys->N_edge - sys->bden) + inde] = 0;
                }
                free(tmp); tmp = NULL;
                free(tmp1); tmp1 = NULL;
                free(tmp2); tmp2 = NULL;
                continue;
            }
            /*tau_qr = (lapack_complex_double*)malloc(U0_i * sizeof(lapack_complex_double));

            info = LAPACKE_zgeqrf(LAPACK_COL_MAJOR, (sys->N_edge - 2 * sys->N_edge_s), U0_i, V_re, (sys->N_edge - 2 * sys->N_edge_s), tau_qr);

            info = LAPACKE_zungqr(LAPACK_COL_MAJOR, (sys->N_edge - 2 * sys->N_edge_s), i_re, i_re, V_re, (sys->N_edge - 2 * sys->N_edge_s), tau_qr);

            free(tau_qr); tau_qr = NULL;*/
            t1 = clock();
            free(tmp); tmp = NULL;
            free(tmp1); tmp1 = NULL;
            free(tmp2); tmp2 = NULL;
            l++;
            Cr = (double*)calloc(U0_i * U0_i, sizeof(double));
            D_sigr = (double*)calloc(U0_i * U0_i, sizeof(double));
            D_epsr = (double*)calloc(U0_i * U0_i, sizeof(double));
            index = 0;
            
            while (index < sys->leng_S){
                start = sys->SRowId[index];

                while (index < sys->leng_S && sys->SRowId[index] == start){
                    temp[(U0_i - 1)][sys->SRowId[index]] += sys->Sval[index] * U0[(U0_i - 1)][sys->SColId[index]];
                    index++;
                }
                if (sys->markEdge[sys->mapEdgeR[start]] != 0) {
                    temp1[(U0_i - 1)][start] = sqrt(sys->stackSign[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)]) * U0[U0_i - 1][start];//SIGMA
                    temp2[(U0_i - 1)][start] = sqrt(sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0) * U0[U0_i - 1][start];
                }
                else {
                    temp1[(U0_i - 1)][start] = 0;
                    temp2[(U0_i - 1)][start] = sqrt(sys->stackEpsn[(sys->mapEdgeR[start]) / (sys->N_edge_v + sys->N_edge_s)] * EPSILON0) * U0[U0_i - 1][start];
                }
            }
            index = U0_i - 1;
            Cr_p.push_back(Cr_base);
            D_sigr_p.push_back(Cr_base);
            D_epsr_p.push_back(Cr_base);
            for (myint inde2 = 0; inde2 < U0_i; inde2++){
                for (myint inde3 = 0; inde3 < (sys->N_edge - sys->bden); inde3++){
                    Cr[index * U0_i + inde2] += U0[inde2][inde3] * temp[index][inde3];
                    D_sigr[index * U0_i + inde2] += temp1[inde2][inde3] * temp1[index][inde3];
                    D_epsr[index * U0_i + inde2] += temp2[inde2][inde3] * temp2[index][inde3];
                }
                Cr_p[index].push_back(Cr[index * U0_i + inde2]);
                D_sigr_p[index].push_back(D_sigr[index * U0_i + inde2]);
                D_epsr_p[index].push_back(D_epsr[index * U0_i + inde2]);

            }

            index = U0_i - 1;
            for (myint inde2 = 0; inde2 < U0_i - 1; inde2++){
                for (myint inde3 = 0; inde3 < (sys->N_edge - sys->bden); inde3++){
                    Cr[inde2 * U0_i + index] += U0[index][inde3] * temp[inde2][inde3];
                    D_sigr[inde2 * U0_i + index] += temp1[index][inde3] * temp1[inde2][inde3];
                    D_epsr[inde2 * U0_i + index] += temp2[index][inde3] * temp2[inde2][inde3];
                    
                }
                Cr_p[inde2].push_back(Cr[inde2 * U0_i + index]);
                D_sigr_p[inde2].push_back(D_sigr[inde2 * U0_i + index]);
                D_epsr_p[inde2].push_back(D_epsr[inde2 * U0_i + index]);
            }
            for (myint inde = 0; inde < U0_i - 1; inde++){
                for (myint inde2 = 0; inde2 < U0_i - 1; inde2++){
                    Cr[inde * U0_i + inde2] += Cr_p[inde][inde2];
                    D_sigr[inde * U0_i + inde2] += D_sigr_p[inde][inde2];
                    D_epsr[inde * U0_i + inde2] += D_epsr_p[inde][inde2];
                }
            }
            
            Ar = (double*)malloc(4 * U0_i * U0_i * sizeof(double));
            Br = (double*)malloc(4 * U0_i * U0_i * sizeof(double));
            for (myint inde = 0; inde < U0_i; inde++){
                for (myint inde2 = 0; inde2 < U0_i; inde2++){
                    Ar[inde * 2 * U0_i + inde2] = -Cr[inde * U0_i + inde2];
                    Br[inde * 2 * U0_i + inde2] = D_sigr[inde * U0_i + inde2] * scale;
                }
            }
            for (myint inde = 0; inde < U0_i; inde++){
                for (myint inde2 = U0_i; inde2 < 2 * U0_i; inde2++){
                    Ar[inde * 2 * U0_i + inde2] = 0;
                    Br[inde * 2 * U0_i + inde2] = D_epsr[inde * U0_i + inde2 - U0_i] * pow(scale, 2);
                }
            }
            for (myint inde = U0_i; inde < 2 * U0_i; inde++){
                for (myint inde2 = 0; inde2 < U0_i; inde2++){
                    Ar[inde * 2 * U0_i + inde2] = 0;
                    Br[inde * 2 * U0_i + inde2] = D_epsr[(inde - U0_i) * U0_i + inde2] * pow(scale, 2);
                }
            }
            for (myint inde = U0_i; inde < 2 * U0_i; inde++){
                for (myint inde2 = U0_i; inde2 < 2 * U0_i; inde2++){
                    Ar[inde * 2 * U0_i + inde2] = D_epsr[(inde - U0_i) * U0_i + inde2 - U0_i] * pow(scale, 2);
                    Br[inde * 2 * U0_i + inde2] = 0;
                }
            }
            
            /*if (ind < 1000){
                cout << "Ar is " << endl;
                for (myint inde = 0; inde < 2 * U0_i; inde++){
                    for (myint inde2 = 0; inde2 < 2 * U0_i; inde2++){
                        cout << Ar[inde2 * 2 * U0_i + inde] << ", ";
                    }
                    cout << endl;
                }
                cout << endl;
                cout << "Br is " << endl;
                for (myint inde = 0; inde < 2 * U0_i; inde++){
                    for (myint inde2 = 0; inde2 < 2 * U0_i; inde2++){
                        cout << Br[inde2 * 2 * U0_i + inde] << ", ";
                    }
                    cout << endl;
                }
            }*/
            
            
            n = 2 * U0_i;
            jobvl = 'N'; jobvr = 'V';
            lda = 2 * U0_i; ldb = 2 * U0_i;
            ldvl = 2 * U0_i; ldvr = 2 * U0_i;
            alphar = (double*)malloc(2 * U0_i * sizeof(double));
            alphai = (double*)malloc(2 * U0_i * sizeof(double));
            beta = (double*)malloc(2 * U0_i * sizeof(double));
            vl1 = (double*)calloc(4 * U0_i * U0_i, sizeof(double));
            vr1 = (double *)calloc(4 * U0_i * U0_i, sizeof(double));
            //lscale = (double*)calloc(2 * U0_i, sizeof(double));
            //rscale = (double*)calloc(2 * U0_i, sizeof(double));
            //tau1 = (double*)calloc(2 * U0_i, sizeof(double));
            //RR = (double*)calloc(4 * U0_i * U0_i, sizeof(double));    // store the upper triangular matrix
            //QQ = (double*)calloc(4 * U0_i * U0_i, sizeof(double));
            //ZZ = (double*)calloc(4 * U0_i * U0_i, sizeof(double));
            //select = (int*)calloc(2 * U0_i, sizeof(int));
            dp1.clear();
            cout << "Generate eigenvalue problem time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
            t1 = clock();
            info = LAPACKE_dggev(LAPACK_COL_MAJOR, jobvl, jobvr, 2 * U0_i, Ar, lda, Br, ldb, alphar, alphai, beta, vl1, ldvl, vr1, ldvr);
            /*info = LAPACKE_dggbal(LAPACK_COL_MAJOR, 'B', 2 * U0_i, Ar, lda, Br, ldb, &ilo, &ihi, lscale, rscale);
            info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, 2 * U0_i, 2 * U0_i, Br, 2 * U0_i, tau1);
            for (myint inde = 0; inde < 2 * U0_i; inde++){
                for (myint inde2 = inde; inde2 < 2 * U0_i; inde2++){
                    RR[inde2 * 2 * U0_i + inde] = Br[inde2 * U0_i + inde];
                    Br[inde2 * 2 * U0_i + inde] = 0;
                }
            }
            info = LAPACKE_dormqr(LAPACK_COL_MAJOR, 'L', 'T', 2 * U0_i, 2 * U0_i, 2 * U0_i, Br, 2 * U0_i, tau1, Ar, 2 * U0_i);
            info = LAPACKE_dgghrd(LAPACK_COL_MAJOR, 'I', 'I', 2 * U0_i, ilo, ihi, Ar, 2 * U0_i, RR, 2 * U0_i, QQ, 2 * U0_i, ZZ, 2 * U0_i);
            info = LAPACKE_dhgeqz(LAPACK_COL_MAJOR, 'S', 'V', 'V', 2 * U0_i, ilo, ihi, Ar, 2 * U0_i, RR, 2 * U0_i, alphar, alphai, beta, QQ, 2 * U0_i, ZZ, 2 * U0_i);
            info = LAPACKE_dtgevc(LAPACK_COL_MAJOR, 'R', 'B', select, 2 * U0_i, Ar, 2 * U0_i, RR, 2 * U0_i, QQ, 2 * U0_i, ZZ, 2 * U0_i, 2 * U0_i, &m);
            info = LAPACKE_dggbak(LAPACK_COL_MAJOR, 'B', 'R', 2 * U0_i, ilo, ihi, lscale, rscale, m, vr1, 2 * U0_i);*/
            
            

            for (myint inde = 0; inde < 2 * U0_i; inde++){
                //if (alphai[inde] == 0){
                //    dp1.push_back(make_pair(alphar[inde] / beta[inde] * scale, inde));
                //}
                //else{
                if (abs(alphai[inde] / beta[inde]) > eps) {    // only those imaginary part larger than 0 eigenvalues are considered
                    dp1.push_back(make_pair(alphar[inde] / beta[inde] * scale + 1i * alphai[inde] / beta[inde] * scale, inde));
                    inde++;
                    dp1.push_back(make_pair(alphar[inde] / beta[inde] * scale + 1i * alphai[inde] / beta[inde] * scale, inde));
                }
            }
            /*if (ind < 1000){
                

                outfile.open("vr.txt", std::ofstream::out | std::ofstream::trunc);
                for (myint inde = 0; inde < 2 * U0_i; inde++){
                    for (myint inde1 = 0; inde1 < 2 * U0_i; inde1++){
                        outfile << vr1[inde1 * (2 * U0_i) + inde] << " ";
                    }
                    outfile << endl;
                }
                outfile.close();

                for (myint inde = 0; inde < dp1.size(); inde++){
                    cout << dp1[inde].first << " ";
                }
                cout << endl;
            }*/
            sort(dp1.begin(), dp1.end(), comp);
            
            ind_i = 0;
            ind_j = 0;
            i_re = 0;
            i_nre = 0;
            re.clear();
            nre.clear();
            
                /*for (myint inde = 0; inde < dp1.size(); inde++){
                    cout << dp1[inde].first << " ";
                }
                cout << endl;*/
                /* for (myint inde = 0; inde < 2 * U0_i; inde++){
                    for (myint inde2 = 0; inde2 < 2 * U0_i; inde2++){
                        cout << vr1[inde2 * 2 * U0_i + inde] << " ";
                    }
                    cout << endl;
                }*/
                /*if (ind == 3 * SG)
                    break;*/
            
            while (ind_i < dp.size() && abs(dp[ind_i].first) < zero_value){
                ind_i++;
            }
            while (ind_j < dp1.size() && abs(dp1[ind_j].first) < zero_value){
                ind_j++;
            }
            while (ind_i < dp.size() && ind_j < dp1.size()){
                if (abs(dp[ind_i].first) > maxvalue || abs(dp1[ind_j].first) > maxvalue){
                    break;
                }
                if (abs(abs(dp[ind_i].first) - abs(dp1[ind_j].first)) / abs(dp[ind_i].first) < eps2){
                    if (dp[ind_i].first.imag() == 0){    // if the eigenvalue is real, only corresponds to one eigenvector
                        re.push_back(ind_i);    // put the eigenvector index into re
                        ind_i++;
                        i_re++;
                    }
                    else if (abs(dp[ind_i].first.imag() + dp[ind_i + 1].first.imag()) / abs(dp[ind_i].first.imag()) < 1.e-3){    // if the eigenvalue is imaginary, corresponds to two eigenvectors
                        if (dp[ind_i].first.imag() > 0){    // first put the positive eigenvalue and then put the negative eigenvalue
                            re.push_back(ind_i);
                            re.push_back(ind_i + 1);
                        }
                        else if (dp[ind_i].first.imag() < 0){
                            re.push_back(ind_i + 1);
                            re.push_back(ind_i);
                        }
                        
                        ind_i += 2;    // two conjugate eigenvalues
                        i_re += 2;
                    }
                }
                
                else if (abs(abs(dp[ind_i].first) - abs(dp1[ind_j].first)) / abs(dp[ind_i].first) >= eps2 && abs(dp[ind_i].first) > abs(dp1[ind_j].first)){
                    ind_j++;
                }
                else if (abs(abs(dp[ind_i].first) - abs(dp1[ind_j].first)) / abs(dp[ind_i].first) >= eps2 && abs(dp[ind_i].first) < abs(dp1[ind_j].first)){
                    if (dp[ind_i].first.imag() == 0){
                        nre.push_back(ind_i);
                        ind_i++;
                        i_nre++;
                    }
                    else{
                        if (abs(dp[ind_i].first.imag() + dp[ind_i + 1].first.imag()) / abs(dp[ind_i].first.imag()) < 1.e-3){
                            if (dp[ind_i].first.imag() > 0){    // first put the positive eigenvalue and then put the negative eigenvalue
                                nre.push_back(ind_i);
                                nre.push_back(ind_i + 1);
                            }
                            else if (dp[ind_i].first.imag() < 0){
                                nre.push_back(ind_i + 1);
                                nre.push_back(ind_i);
                            }
                        }
                        ind_i += 2;    // two conjugate eigenvalues
                        i_nre += 2;
                    }
                }
            }
            
            while (ind_i < dp.size()){
                if (dp[ind_i].first.imag() == 0){
                    nre.push_back(ind_i);
                    ind_i++;
                    i_nre++;
                }
                else if (abs(dp[ind_i].first.imag() + dp[ind_i + 1].first.imag()) / abs(dp[ind_i].first.imag()) < 1.e-5){
                    if (dp[ind_i].first.imag() > 0){    // first put the positive eigenvalue and then put the negative eigenvalue
                        nre.push_back(ind_i);
                        nre.push_back(ind_i + 1);
                    }
                    else if (dp[ind_i].first.imag() < 0){
                        nre.push_back(ind_i + 1);
                        nre.push_back(ind_i);
                    }
                    ind_i += 2;    // two conjugate eigenvalues
                    i_nre += 2;
                }
                else if (abs(dp[ind_i].first) > maxvalue){
                    nre.push_back(ind_i);
                    ind_i++;
                    i_nre++;
                }

            }
            /*for (myint inde = 0; inde < dp.size(); inde++){
                cout << dp[inde].first << " ";
            }
            cout << endl;*/
            free(V_re);
            V_re = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * i_re, sizeof(lapack_complex_double));
            free(V_nre);
            V_nre = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * i_nre, sizeof(lapack_complex_double));

            i_re = 0;
            
            for (myint rei = 0; rei < re.size(); rei++){
                if (dp[re[rei]].first.imag() == 0){    // one eigenvector corresponding to real eigenvalue
                    for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                        for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                            V_re[i_re * (sys->N_edge - sys->bden) + inde].real += U0[inde3][inde] * vr[dp[re[rei]].second * 2 * (U0_i - 1) + inde3];    // only get the first half eigenvectors
                        }
                    }
                    i_re++;
                }
                else{    // two eigenvectors corresponding to imaginary eigenvalue
                    if (dp[re[rei]].first.imag() > 0){
                        for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                            for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                                V_re[i_re * (sys->N_edge - sys->bden) + inde].real = V_re[i_re * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[dp[re[rei]].second * 2 * (U0_i - 1) + inde3];
                                V_re[i_re * (sys->N_edge - sys->bden) + inde].imag = V_re[i_re * (sys->N_edge - sys->bden) + inde].imag + U0[inde3][inde] * vr[dp[re[rei + 1]].second * 2 * (U0_i - 1) + inde3];
                                V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].real = V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[(dp[re[rei]].second) * 2 * (U0_i - 1) + inde3];
                                V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].imag = V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].imag - U0[inde3][inde] * vr[dp[re[rei + 1]].second * 2 * (U0_i - 1) + inde3];
                            }
                        }
                        rei++;
                    }
                    else if (dp[re[rei]].first.imag() < 0){
                        for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                            for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                                V_re[i_re * (sys->N_edge - sys->bden) + inde].real = V_re[i_re * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[dp[re[rei + 1]].second * 2 * (U0_i - 1) + inde3];
                                V_re[i_re * (sys->N_edge - sys->bden) + inde].imag = V_re[i_re * (sys->N_edge - sys->bden) + inde].imag - U0[inde3][inde] * vr[dp[re[rei]].second * 2 * (U0_i - 1) + inde3];
                                V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].real = V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[dp[re[rei + 1]].second * 2 * (U0_i - 1) + inde3];
                                V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].imag = V_re[(i_re + 1) * (sys->N_edge - sys->bden) + inde].imag + U0[inde3][inde] * vr[dp[re[rei]].second * 2 * (U0_i - 1) + inde3];
                            }
                        }
                        rei++;
                    }
                    i_re += 2;
                }
            }

            
            i_nre = 0;
            for (myint nrei = 0; nrei < nre.size(); nrei++){
                if (abs(dp[nre[nrei]].first) > maxvalue){    // if this eigenvalue is infinity
                    for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                        for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                            V_nre[i_nre * (sys->N_edge - sys->bden) + inde].real += U0[inde3][inde] * vr[dp[nre[nrei]].second * 2 * (U0_i - 1) + inde3];
                        }
                    }
                    i_nre++;
                    continue;
                }
                if (dp[nre[nrei]].first.imag() == 0){    // one eigenvector corresponding to real eigenvalue
                    for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                        for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                            V_nre[i_nre * (sys->N_edge - sys->bden) + inde].real += U0[inde3][inde] * vr[dp[nre[nrei]].second * 2 * (U0_i - 1) + inde3];
                        }
                    }
                    i_nre++;
                }
                else{    // two eigenvectors corresponding to imaginary eigenvalue
                    if (dp[nre[nrei]].first.imag() > 0){
                        for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                            for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                                V_nre[i_nre * (sys->N_edge - sys->bden) + inde].real = V_nre[i_nre * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[dp[nre[nrei]].second * 2 * (U0_i - 1) + inde3];
                                V_nre[i_nre * (sys->N_edge - sys->bden) + inde].imag = V_nre[i_nre * (sys->N_edge - sys->bden) + inde].imag + U0[inde3][inde] * vr[dp[nre[nrei + 1]].second * 2 * (U0_i - 1) + inde3];
                                V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].real = V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[(dp[nre[nrei]].second) * 2 * (U0_i - 1) + inde3];
                                V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].imag = V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].imag - U0[inde3][inde] * vr[dp[nre[nrei + 1]].second * 2 * (U0_i - 1) + inde3];
                            }
                        }
                        nrei++;
                    }
                    else if (dp[nre[nrei]].first.imag() < 0){
                        for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                            for (myint inde3 = 0; inde3 < U0_i - 1; inde3++){
                                V_nre[i_nre * (sys->N_edge - sys->bden) + inde].real = V_nre[i_nre * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[dp[nre[nrei + 1]].second * 2 * (U0_i - 1) + inde3];
                                V_nre[i_nre * (sys->N_edge - sys->bden) + inde].imag = V_nre[i_nre * (sys->N_edge - sys->bden) + inde].imag - U0[inde3][inde] * vr[dp[nre[nrei]].second * 2 * (U0_i - 1) + inde3];
                                V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].real = V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].real + U0[inde3][inde] * vr[dp[nre[nrei + 1]].second * 2 * (U0_i - 1) + inde3];
                                V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].imag = V_nre[(i_nre + 1) * (sys->N_edge - sys->bden) + inde].imag + U0[inde3][inde] * vr[dp[nre[nrei]].second * 2 * (U0_i - 1) + inde3];
                            }
                        }
                        nrei++;
                    }
                    i_nre += 2;
                }
            }
            cout << "Solve generalized eigenvale problem time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
            /*cout << "The step No is " << ind << endl;
            cout << "The repeated eigenvalues are of No. " << i_re << endl;
            cout << "The non-repeated eigenvalues are of No. " << i_nre << endl;*/
            
            t1 = clock();
            // orthogonalized V_re
            if (i_re > 0){
                /*outfile.open("V_re.txt", std::ofstream::out | std::ofstream::trunc);
                for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
                    for (myint inde2 = 0; inde2 < i_re; inde2++){
                        outfile << V_re[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " " << V_re[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag << " ";
                    }
                    outfile << endl;
                }
                outfile.close();*/
                /*tau_qr = (lapack_complex_double*)malloc(i_re * sizeof(lapack_complex_double));
                RR = (lapack_complex_double*)calloc(i_re * i_re, sizeof(lapack_complex_double));

                info = LAPACKE_zgeqrf(LAPACK_COL_MAJOR, (sys->N_edge - 2 * sys->N_edge_s), i_re, V_re, (sys->N_edge - 2 * sys->N_edge_s), tau_qr);
                
                for (myint inde = 0; inde < i_re; inde++){
                    for (myint inde1 = 0; inde1 <= inde; inde1++){
                        RR[inde * i_re + inde1].real = V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].real;
                        RR[inde * i_re + inde1].imag = V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].imag;
                    }
                }
                info = LAPACKE_zungqr(LAPACK_COL_MAJOR, (sys->N_edge - 2 * sys->N_edge_s), i_re, i_re, V_re, (sys->N_edge - 2 * sys->N_edge_s), tau_qr);
                
                free(tau_qr); tau_qr = NULL;*/

               /* outfile.open("RR.txt", std::ofstream::out | std::ofstream::trunc);
                for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
                    for (myint inde2 = 0; inde2 < i_re; inde2++){
                        outfile << RR[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " " << RR[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag << " ";
                    }
                    outfile << endl;
                }
                outfile.close();*/
                

                /*i_re1 = 0;
                for (myint inde = 0; inde < i_re; inde++){
                    if (sqrt(pow(RR[inde * i_nre + inde].real, 2) + pow(RR[inde * i_nre + inde].imag, 2)) > 1.e-12){
                        i_re1++;
                    }
                }
                free(V_re1);
                V_re1 = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s) * i_re1, sizeof(lapack_complex_double));
                i_re1 = 0;
                for (myint inde = 0; inde < i_re; inde++){
                    if (sqrt(pow(RR[inde * i_nre + inde].real, 2) + pow(RR[inde * i_nre + inde].imag, 2)) > 1.e-12){
                        for (myint inde1 = 0; inde1 < sys->N_edge - 2 * sys->N_edge_s; inde1++){
                            V_re1[i_re1 * (sys->N_edge - 2 * sys->N_edge_s) + inde1].real = V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].real;
                            V_re1[i_re1 * (sys->N_edge - 2 * sys->N_edge_s) + inde1].imag = V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].imag;
                        }
                        i_re1++;
                    }
                }
                free(RR); RR = NULL;*/
                /*outfile.open("V_re1.txt", std::ofstream::out | std::ofstream::trunc);
                for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
                    for (myint inde2 = 0; inde2 < i_re; inde2++){
                        outfile << V_re[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " " << V_re[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag << " ";
                    }
                    outfile << endl;
                }
                outfile.close();*/


                for (myint inde = 0; inde < i_re; inde++){
                    nn = 0;
                    for (myint inde1 = 0; inde1 < sys->N_edge - sys->bden; inde1++){
                        nn += pow(V_re[inde * (sys->N_edge - sys->bden) + inde1].real, 2) + pow(V_re[inde * (sys->N_edge - sys->bden) + inde1].imag, 2);
                    }
                    nn = sqrt(nn);
                    for (myint inde1 = 0; inde1 < sys->N_edge - sys->bden; inde1++){
                        V_re[inde * (sys->N_edge - sys->bden) + inde1].real = V_re[inde * (sys->N_edge - sys->bden) + inde1].real / nn;
                        V_re[inde * (sys->N_edge - sys->bden) + inde1].imag = V_re[inde * (sys->N_edge - sys->bden) + inde1].imag / nn;
                    }
                    /*for (myint inde1 = inde + 1; inde1 < i_re; inde1++){
                        nnr = 0;
                        nni = 0;
                        for (myint inde2 = 0; inde2 < sys->N_edge - 2 * sys->N_edge_s; inde2++){
                            nnr += V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real * V_re[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real
                                + V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag * V_re[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag;
                            nni += -V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag * V_re[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real
                                +V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real * V_re[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag;
                        }
                        
                        qc = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s), sizeof(lapack_complex_double));
                        for (myint inde2 = 0; inde2 < sys->N_edge - 2 * sys->N_edge_s; inde2++){
                            qc[inde2].real = V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real * nnr - nni * V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag;
                            qc[inde2].imag = nnr * V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag + nni * V_re[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real;
                            
                           
                            V_re[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real -= qc[inde2].real;
                            V_re[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag -= qc[inde2].imag;

                            
                        }
                        free(qc); qc = NULL;
                    }*/
                    
                }
                
                

                if (i_nre > 0){
                    tmp3 = (lapack_complex_double*)calloc(i_re * i_nre, sizeof(lapack_complex_double));
                    tmp4 = (lapack_complex_double*)calloc((sys->N_edge - sys->bden) * i_nre, sizeof(lapack_complex_double));
                    m_re = (lapack_complex_double*)calloc(i_re * i_re, sizeof(lapack_complex_double));
                    
                    status = matrix_multi('T', V_re, (sys->N_edge - sys->bden), i_re, V_nre, (sys->N_edge - sys->bden), i_nre, tmp3);
                    status = matrix_multi('T', V_re, (sys->N_edge - sys->bden), i_re, V_re, (sys->N_edge - sys->bden), i_re, m_re);
                    ipiv = (lapack_int*)malloc((i_nre + i_re) * sizeof(lapack_int));
                    info = LAPACKE_zgesv(LAPACK_COL_MAJOR, i_re, i_nre, m_re, i_re, ipiv, tmp3, i_re);// , y_nre1, i_nre + i_re, &iter);
                    status = matrix_multi('N', V_re, (sys->N_edge - sys->bden), i_re, tmp3, i_re, i_nre, tmp4);
                    free(ipiv); ipiv = NULL;
                    free(m_re); m_re = NULL;

                    //V_nre_norm = (double*)calloc(i_nre, sizeof(double));
                    for (myint inde2 = 0; inde2 < i_nre; inde2++){
                        for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                            V_nre[inde2 * (sys->N_edge - sys->bden) + inde].real -= tmp4[inde2 * (sys->N_edge - sys->bden) + inde].real;
                            V_nre[inde2 * (sys->N_edge - sys->bden) + inde].imag -= tmp4[inde2 * (sys->N_edge - sys->bden) + inde].imag;
                            //V_nre_norm[inde2] += pow(V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real, 2) + pow(V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag, 2);
                        }
                        //V_nre_norm[inde2] = sqrt(V_nre_norm[inde2]);
                        //cout << V_nre_norm[inde2] << endl;
                    }

                    // orthogonalize V_nre by using Modified GS
                    for (myint inde = 0; inde < i_nre; inde++){
                        nn = 0;
                        for (myint inde1 = 0; inde1 < sys->N_edge - sys->bden; inde1++){
                            nn += pow(V_nre[inde * (sys->N_edge - sys->bden) + inde1].real, 2) + pow(V_nre[inde * (sys->N_edge - sys->bden) + inde1].imag, 2);
                        }
                        nn = sqrt(nn);
                        for (myint inde1 = 0; inde1 < sys->N_edge - sys->bden; inde1++){
                            V_nre[inde * (sys->N_edge - sys->bden) + inde1].real = V_nre[inde * (sys->N_edge - sys->bden) + inde1].real / nn;
                            V_nre[inde * (sys->N_edge - sys->bden) + inde1].imag = V_nre[inde * (sys->N_edge - sys->bden) + inde1].imag / nn;
                        }
                        /*for (myint inde1 = inde + 1; inde1 < i_nre; inde1++){
                            nnr = 0;
                            nni = 0;
                            for (myint inde2 = 0; inde2 < sys->N_edge - 2 * sys->N_edge_s; inde2++){
                                nnr += V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real * V_nre[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real
                                    + V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag * V_nre[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag;
                                nni += -V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag * V_nre[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real
                                    + V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real * V_nre[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag;
                            }

                            qc = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s), sizeof(lapack_complex_double));
                            for (myint inde2 = 0; inde2 < sys->N_edge - 2 * sys->N_edge_s; inde2++){
                                qc[inde2].real = V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real * nnr - nni * V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag;
                                qc[inde2].imag = nnr * V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag + nni * V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real;


                                V_nre[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].real -= qc[inde2].real;
                                V_nre[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde2].imag -= qc[inde2].imag;


                            }
                            free(qc); qc = NULL;
                        }*/

                    }

                   // make each column vector to be of norm 1
                   /*for (myint inde2 = 0; inde2 < i_nre; inde2++){

                       for (myint inde = 0; inde < (sys->N_edge - 2 * sys->N_edge_s); inde++){
                           V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real = V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real / V_nre_norm[inde2];
                           V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag = V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag / V_nre_norm[inde2];
                           
                       }
                   }
                   free(V_nre_norm); V_nre_norm = NULL;*/
                   
                   /*tau_qr = (lapack_complex_double*)malloc(i_nre * sizeof(lapack_complex_double));
                   RR = (lapack_complex_double*)calloc(i_nre * i_nre, sizeof(lapack_complex_double));
                   info = LAPACKE_zgeqrf(LAPACK_COL_MAJOR, (sys->N_edge - 2 * sys->N_edge_s), i_nre, V_nre, (sys->N_edge - 2 * sys->N_edge_s), tau_qr);
                   for (myint inde = 0; inde < i_nre; inde++){
                       for (myint inde1 = 0; inde1 <= inde; inde1++){
                           RR[inde * i_nre + inde1].real = V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].real;
                           RR[inde * i_nre + inde1].imag = V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].imag;
                       }
                   }
                   info = LAPACKE_zungqr(LAPACK_COL_MAJOR, (sys->N_edge - 2 * sys->N_edge_s), i_nre, i_nre, V_nre, (sys->N_edge - 2 * sys->N_edge_s), tau_qr);
                   free(tau_qr); tau_qr = NULL;
                   i_nre1 = 0;
                   for (myint inde = 0; inde < i_nre; inde++){
                       if (sqrt(pow(RR[inde * i_nre + inde].real, 2) + pow(RR[inde * i_nre + inde].imag, 2)) > 1.e-12){
                           i_nre1++;
                       }
                   }
                   free(V_nre1);
                   V_nre1 = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s) * i_nre1, sizeof(lapack_complex_double));
                   i_nre1 = 0;
                   for (myint inde = 0; inde < i_nre; inde++){
                       if (sqrt(pow(RR[inde * i_nre + inde].real, 2) + pow(RR[inde * i_nre + inde].imag, 2)) > 1.e-12){
                           for (myint inde1 = 0; inde1 < sys->N_edge - 2 * sys->N_edge_s; inde1++){
                               V_nre1[i_nre1 * (sys->N_edge - 2 * sys->N_edge_s) + inde1].real = V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].real;
                               V_nre1[i_nre1 * (sys->N_edge - 2 * sys->N_edge_s) + inde1].imag = V_nre[inde * (sys->N_edge - 2 * sys->N_edge_s) + inde1].imag;
                           }
                           i_nre1++;
                       }
                   }
                   free(RR); RR = NULL;*/
                        /*outfile.open("V_nre.txt", std::ofstream::out | std::ofstream::trunc);
                        for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
                            for (myint inde2 = 0; inde2 < i_nre; inde2++){
                                outfile << V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " " << V_nre[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag << " ";
                            }
                            outfile << endl;
                        }
                        outfile.close();*/


                    free(tmp3); tmp3 = NULL;
                    free(tmp4); tmp4 = NULL;
                    cout << "Repeated eigenvalues number is " << i_re << " and non-repeated eigenvalues number is " << i_nre << endl;
                }
            }
            cout << "The time to generate orthogonalized V_re and V_nre is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
            
            if (i_re == 0){
                dp.clear();
                for (myint inde = 0; inde < dp1.size(); inde++){
                    dp.push_back(dp1[inde]);
                }
                free(vr); vr = NULL;
                vr = (double*)malloc(4 * U0_i * U0_i * sizeof(double));
                for (myint inde = 0; inde < 2 * U0_i; inde++){
                    for (myint inde2 = 0; inde2 < 2 * U0_i; inde2++){
                        vr[inde * 2 * U0_i + inde2] = vr1[inde * 2 * U0_i + inde2];
                    }
                }
                free(vr1); vr1 = NULL;
                free(vl1); vl1 = NULL;
                //free(V_re); V_re = NULL;
                //free(V_nre); V_nre = NULL;
                free(Cr); Cr = NULL;
                free(D_epsr); D_epsr = NULL;
                free(D_sigr); D_sigr = NULL;
                free(Ar); Ar = NULL;
                free(Br); Br = NULL;
                free(alphar); alphar = NULL;
                free(alphai); alphai = NULL;
                free(beta); beta = NULL;
                for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                    xr[inde] = xr[1 * (sys->N_edge - sys->bden) + inde];
                    xr[1 * (sys->N_edge - sys->bden) + inde] = xr[2 * (sys->N_edge - sys->bden) + inde];
                    xr[2 * (sys->N_edge - sys->bden) + inde] = 0;
                }
                /*if (ind >= 5000)
                    break;*/
                cout << endl;
                continue;
            }
            else if (i_nre == 0){    // all the eigenvalues are repeated
                //free(V_nre); V_nre = NULL;
                //free(V_re); V_re = NULL;
                //break;
                dp.clear();
                for (myint inde = 0; inde < dp1.size(); inde++){
                    dp.push_back(dp1[inde]);
                }
                free(vr); vr = NULL;
                vr = (double*)malloc(4 * U0_i * U0_i * sizeof(double));
                for (myint inde = 0; inde < 2 * U0_i; inde++){
                    for (myint inde2 = 0; inde2 < 2 * U0_i; inde2++){
                        vr[inde * 2 * U0_i + inde2] = vr1[inde * 2 * U0_i + inde2];
                    }
                }
                free(vr1); vr1 = NULL;
                free(vl1); vl1 = NULL;
                //free(V_re); V_re = NULL;
                //free(V_nre); V_nre = NULL;
                free(Cr); Cr = NULL;
                free(D_epsr); D_epsr = NULL;
                free(D_sigr); D_sigr = NULL;
                free(Ar); Ar = NULL;
                free(Br); Br = NULL;
                free(alphar); alphar = NULL;
                free(alphai); alphai = NULL;
                free(beta); beta = NULL;
                for (myint inde = 0; inde < (sys->N_edge - sys->bden); inde++){
                    xr[inde] = xr[1 * (sys->N_edge - sys->bden) + inde];
                    xr[1 * (sys->N_edge - sys->bden) + inde] = xr[2 * (sys->N_edge - sys->bden) + inde];
                    xr[2 * (sys->N_edge - sys->bden) + inde] = 0;
                }
                /*if (ind >= 5000)
                    break;*/
                cout << endl;
                continue;
            }
            else{
                t1 = clock();
                y_new = (lapack_complex_double*)calloc((i_re + i_nre), sizeof(lapack_complex_double));
                for (myint inde = 0; inde < i_re; inde++){
                    for (myint inde2 = 0; inde2 < (sys->N_edge - sys->bden); inde2++){
                        y_new[inde].real = y_new[inde].real + (V_re[inde * (sys->N_edge - sys->bden) + inde2].real) * xr[1 * (sys->N_edge - sys->bden) + inde2];    //conjugate transpose
                        y_new[inde].imag = y_new[inde].imag - (V_re[inde * (sys->N_edge - sys->bden) + inde2].imag) * xr[1 * (sys->N_edge - sys->bden) + inde2];    //conjugate transpose

                    }
                }
                for (myint inde = 0; inde < i_nre; inde++){
                    for (myint inde2 = 0; inde2 < (sys->N_edge - sys->bden); inde2++){
                        y_new[inde + i_re].real = y_new[inde + i_re].real + V_nre[inde * (sys->N_edge - sys->bden) + inde2].real * xr[1 * (sys->N_edge - sys->bden) + inde2];    // conjugate transpose
                        y_new[inde + i_re].imag = y_new[inde + i_re].imag - V_nre[inde * (sys->N_edge - sys->bden) + inde2].imag * xr[1 * (sys->N_edge - sys->bden) + inde2];
                    }
                }
                
                
                V_new = (lapack_complex_double*)malloc((sys->N_edge - sys->bden) * (i_re + i_nre) * sizeof(lapack_complex_double));
                for (myint inde = 0; inde < i_re; inde++){
                    for (myint inde2 = 0; inde2 < sys->N_edge - sys->bden; inde2++){
                        V_new[inde * (sys->N_edge - sys->bden) + inde2].real = V_re[inde * (sys->N_edge - sys->bden) + inde2].real;
                        V_new[inde * (sys->N_edge - sys->bden) + inde2].imag = V_re[inde * (sys->N_edge - sys->bden) + inde2].imag;
                    }
                }
                for (myint inde = i_re; inde < i_re + i_nre; inde++){
                    for (myint inde2 = 0; inde2 < sys->N_edge - sys->bden; inde2++){
                        V_new[inde * (sys->N_edge - sys->bden) + inde2].real = V_nre[(inde - i_re) * (sys->N_edge - sys->bden) + inde2].real;
                        V_new[inde * (sys->N_edge - sys->bden) + inde2].imag = V_nre[(inde - i_re) * (sys->N_edge - sys->bden) + inde2].imag;
                    }
                }
                

                m_new = (lapack_complex_double*)calloc((i_nre + i_re) * (i_nre + i_re), sizeof(lapack_complex_double));
                ipiv = (lapack_int*)malloc((i_nre + i_re) * sizeof(lapack_int));
                //y_nre1 = (lapack_complex_double*)calloc(i_nre + i_re, sizeof(lapack_complex_double));
                status = matrix_multi('T', V_new, (sys->N_edge - sys->bden), i_nre + i_re, V_new, (sys->N_edge - sys->bden), i_nre + i_re, m_new);
                
                info = LAPACKE_zgesv(LAPACK_COL_MAJOR, i_nre + i_re, 1, m_new, i_nre + i_re, ipiv, y_new, i_nre + i_re);// , y_nre1, i_nre + i_re, &iter);
                x_re = 0;
                x_nre = 0;
                for (myint inde = 0; inde < i_re; inde++){
                    x_re += pow(y_new[inde].real, 2) + pow(y_new[inde].imag, 2);
                }
                for (myint inde = i_re; inde < i_nre + i_re; inde++){
                    x_nre += pow(y_new[inde].real, 2) + pow(y_new[inde].imag, 2);
                }
                cout << "The weight on step " << ind << " is " << x_nre / x_re << endl;
                //cout << "U0 size is " << U0_i << endl;
                cout << "Weight calculation time is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
                cout << endl;
                free(V_new); V_new = NULL;
                free(y_new); y_new = NULL;
                
                //free(y_nre1); y_nre1 = NULL;
                free(ipiv); ipiv = NULL;
                
                free(m_new); m_new = NULL;
                
                //free(V_nre); V_nre = NULL;
                free(Cr); Cr = NULL;
                
                free(D_epsr); D_epsr = NULL;
                
                free(D_sigr); D_sigr = NULL;
                free(Ar); Ar = NULL;
                free(Br); Br = NULL;
                free(alphar); alphar = NULL;
                free(alphai); alphai = NULL;
                free(beta); beta = NULL;
                if  ((x_nre / x_re) < eps1){
                    free(vr1); vr1 = NULL;
                    free(vl1); vl1 = NULL;
                    break;
                }
                else{
                    dp.clear();
                    for (myint inde = 0; inde < dp1.size(); inde++){
                        dp.push_back(dp1[inde]);
                    }
                    free(vr); vr = NULL;
                    vr = (double*)malloc(4 * U0_i * U0_i * sizeof(double));
                    for (myint inde = 0; inde < 2 * U0_i; inde++){
                        for (myint inde2 = 0; inde2 < 2 * U0_i; inde2++){
                            vr[inde * 2 * U0_i + inde2] = vr1[inde * 2 * U0_i + inde2];
                        }
                    }
                    free(vr1); vr1 = NULL;
                    free(vl1); vl1 = NULL;
                    //free(V_re); V_re = NULL;
                    
                }
                
            }
            
        }


        
        for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++){
            xr[inde] = xr[1 * (sys->N_edge - sys->bden) + inde];
            xr[1 * (sys->N_edge - sys->bden) + inde] = xr[2 * (sys->N_edge - sys->bden) + inde];
            xr[2 * (sys->N_edge - sys->bden) + inde] = 0;
        } 
    }
    free(xr); xr = NULL;
    temp.clear();
    temp1.clear();
    temp2.clear();
    U0.clear();
    /* deduct the u0 space */
    /*outfile.open("u0.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
        for (myint inde1 = 0; inde1 < 2; inde1++){
            outfile << u0[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " ";
        }
        outfile << endl;
    }
    outfile.close();*/
    
    sys->Vh = (lapack_complex_double*)malloc((sys->N_edge - sys->bden) * i_re * sizeof(lapack_complex_double));
    lapack_complex_double *V_re2 = (lapack_complex_double*)malloc((sys->N_edge - sys->bden) * i_re * sizeof(lapack_complex_double));
    for (myint inde = 0; inde < sys->N_edge - sys->bden; inde++){    // A*V_re
        for (myint inde2 = 0; inde2 < i_re; inde2++){
            sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].real = V_re[inde2 * (sys->N_edge - sys->bden) + inde].real;
            sys->Vh[inde2 * (sys->N_edge - sys->bden) + inde].imag = V_re[inde2 * (sys->N_edge - sys->bden) + inde].imag;
            //V_re2[inde2 * (sys->N_edge - sys->bden) + inde].real = V_re[inde2 * (sys->N_edge - sys->bden) + inde].real * (-sys->freqStart * sys->freqUnit * sys->freqStart * sys->freqUnit * 4 * pow(M_PI, 2)) * sys->eps[inde + sys->N_edge_s]
            //    - 2 * M_PI * sys->freqStart * sys->freqUnit * V_re[inde2 * (sys->N_edge - sys->bden) + inde].imag * sys->sig[inde + sys->N_edge_s];
            //V_re2[inde2 * (sys->N_edge - sys->bden) + inde].imag = V_re[inde2 * (sys->N_edge - sys->bden) + inde].real * (sys->freqStart * sys->freqUnit * 2 * M_PI) * sys->sig[inde + sys->N_edge_s]
            //    - V_re[inde2 * (sys->N_edge - sys->bden) + inde].imag * (sys->freqStart * sys->freqUnit * sys->freqStart * sys->freqUnit * 4 * pow(M_PI, 2)) * sys->eps[inde + sys->N_edge_s];
        }
    }

    //y_re = (lapack_complex_double*)calloc(2 * i_re, sizeof(lapack_complex_double));    // u0a'*A*V_re
    //status = matrix_multi('T', u0a, (sys->N_edge - 2 * sys->N_edge_s), 2, V_re2, (sys->N_edge - 2 * sys->N_edge_s), i_re, y_re);

    //tmp3 = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s) * 2, sizeof(lapack_complex_double));
    //for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){    // A*u0
    //    for (myint inde2 = 0; inde2 < 2; inde2++){
    //        tmp3[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real = u0[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real * (-sys->freqStart * sys->freqUnit * sys->freqStart * sys->freqUnit * 4 * pow(M_PI, 2)) * sys->eps[inde + sys->N_edge_s]
    //            - 2 * M_PI * sys->freqStart * sys->freqUnit * u0[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag * sys->sig[inde + sys->N_edge_s];
    //        tmp3[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag = u0[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real * (sys->freqStart * sys->freqUnit * 2 * M_PI) * sys->sig[inde + sys->N_edge_s]
    //            - u0[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag * (sys->freqStart * sys->freqUnit * sys->freqStart * sys->freqUnit * 4 * pow(M_PI, 2)) * sys->eps[inde + sys->N_edge_s];
    //    }
    //}

    //tmp4 = (lapack_complex_double*)calloc(2 * 2, sizeof(lapack_complex_double));    // u0a'*A*u0
    //status = matrix_multi('T', u0a, (sys->N_edge - 2 * sys->N_edge_s), 2, tmp3, (sys->N_edge - 2 * sys->N_edge_s), 2, tmp4);    // u0a'*A*u0
    //ipiv = (lapack_int*)malloc(2 * sizeof(lapack_int));
    //y_new = (lapack_complex_double*)calloc(2 * i_re, sizeof(lapack_complex_double));
    ///*outfile.open("ma.txt", std::ofstream::out | std::ofstream::trunc);
    //for (myint inde = 0; inde < 2; inde++){
    //    for (myint inde1 = 0; inde1 < 2; inde1++){
    //        outfile << tmp4[inde1 * 2 + inde].real << " " << tmp4[inde1 * 2 + inde].imag << " ";
    //    }
    //    outfile << endl;
    //}*/
    //info = LAPACKE_zcgesv(LAPACK_COL_MAJOR, 2, i_re, tmp4, 2, ipiv, y_re, 2, y_new, 2, &iter);    // (u0a'*A*u0)\(u0a'*A*V_re)
    //m_new = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s) * i_re, sizeof(lapack_complex_double));
    //status = matrix_multi('N', u0, (sys->N_edge - 2 * sys->N_edge_s), 2, y_new, 2, i_re, m_new);    // u0*((u0a'*A*u0)\(u0a'*A*V_re))
    //for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
    //    for (myint inde2 = 0; inde2 < i_re; inde2++){
    //        sys->Vh[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real = sys->Vh[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real -m_new[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real;
    //        sys->Vh[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag = sys->Vh[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag -m_new[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag;
    //    }
    //}
    //ofstream out;
    
    /*out.open("Vh.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
        for (myint inde2 = 0; inde2 < i_re; inde2++){
            out << sys->Vh[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " " << sys->Vh[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag << " ";
        }
        out << endl;
    }
    out.close();*/

    /*out.open("u0.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++){
        for (myint inde2 = 0; inde2 < 2; inde2++){
            out << u0[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real << " " << u0[inde2 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag << " ";
        }
        out << endl;
    }
    out.close();*/


    sys->leng_Vh = i_re;

    free(y_re); y_re = NULL;
    free(V_re2); V_re2 = NULL;
    free(tmp3); tmp3 = NULL;
    free(tmp4); tmp4 = NULL;
    free(ipiv); ipiv = NULL;
    free(y_new); y_new = NULL;
    free(m_new); m_new = NULL;
    

    return 0;

}

int matrix_multi(char operation, lapack_complex_double *a, myint arow, myint acol, lapack_complex_double *b, myint brow, myint bcol, lapack_complex_double *tmp3){
    /* operation = 'T' is first matrix conjugate transpose, operation = 'N' is first matrix non-conjugate-transpose*/
    if (operation == 'T'){
        for (myint ind = 0; ind < acol; ind++){
            for (myint ind1 = 0; ind1 < bcol; ind1++){
                tmp3[ind1 * acol + ind].real = 0;
                tmp3[ind1 * acol + ind].imag = 0;
                for (myint ind2 = 0; ind2 < arow; ind2++){
                    tmp3[ind1 * acol + ind].real = tmp3[ind1 * acol + ind].real + a[ind * arow + ind2].real * b[ind1 * brow + ind2].real + a[ind * arow + ind2].imag * b[ind1 * brow + ind2].imag;
                    tmp3[ind1 * acol + ind].imag = tmp3[ind1 * acol + ind].imag - a[ind * arow + ind2].imag * b[ind1 * brow + ind2].real + a[ind * arow + ind2].real * b[ind1 * brow + ind2].imag;
                }
            }
        }
    }
    else if (operation == 'N'){
        for (myint ind = 0; ind < arow; ind++){
            for (myint ind1 = 0; ind1 < bcol; ind1++){
                tmp3[ind1 * arow + ind].real = 0;
                tmp3[ind1 * arow + ind].imag = 0;
                for (myint ind2 = 0; ind2 < acol; ind2++){
                    tmp3[ind1 * arow + ind].real = tmp3[ind1 * arow + ind].real + a[ind2 * arow + ind].real * b[ind1 * brow + ind2].real - a[ind2 * arow + ind].imag * b[ind1 * brow + ind2].imag;
                    tmp3[ind1 * arow + ind].imag = tmp3[ind1 * arow + ind].imag + a[ind2 * arow + ind].imag * b[ind1 * brow + ind2].real + a[ind2 * arow + ind].real * b[ind1 * brow + ind2].imag;
                }
            }
        }
    }
    return 0;
}
