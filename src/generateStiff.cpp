/* Generate the stiffness matrix */
#include "fdtd.hpp"
using namespace std::complex_literals;


int generateStiff(fdtdMesh *sys){
    myint Senum, leng_Se;    // Se's size is (N_patch - N_patch_s) * (N_edge - N_edge_s) only lower boundary is considered as PEC
    myint Shnum, leng_Sh;    // Sh's size is (N_edge - N_edge_s) * (N_patch - N_patch_s) only lower boundary is considered as PEC
    myint *SeRowId, *SeColId;
    double *Seval;
    myint *ShRowId, *ShColId;
    double *Shval;
    myint indx, indy, indz;


    SeRowId = (myint*)malloc((sys->N_patch) * 4 * sizeof(myint));
    SeColId = (myint*)malloc((sys->N_patch) * 4 * sizeof(myint));
    Seval = (double*)malloc((sys->N_patch) * 4 * sizeof(double));
    ShRowId = (myint*)malloc((sys->N_patch) * 4 * sizeof(myint));
    ShColId = (myint*)malloc((sys->N_patch) * 4 * sizeof(myint));
    Shval = (double*)malloc((sys->N_patch) * 4 * sizeof(double));
    Senum = 0;
    Shnum = 0;
    leng_Se = 0;
    leng_Sh = 0;


    /* Generate Se */
    /* the lowest cell layer doesn't contain the lowest plane */

    /* the middle layers */
    for (indz = 1; indz <= sys->N_cell_z; indz++){
        for (indx = 1; indx <= sys->N_cell_x; indx++){
            for (indy = 1; indy <= sys->N_cell_y; indy++){
                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = -1 / (sys->xn[indx] - sys->xn[indx - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = 1 / (sys->xn[indx] - sys->xn[indx - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = 1 / (sys->yn[indy] - sys->yn[indy - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = -1 / (sys->yn[indy] - sys->yn[indy - 1]);
                Senum++;
                leng_Se++;


                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = 1 / (sys->zn[indz] - sys->zn[indz - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = -1 / (sys->yn[indy] - sys->yn[indy - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1 + 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = 1 / (sys->yn[indy] - sys->yn[indy - 1]);
                Senum++;

                SeRowId[Senum] = indz * (sys->N_edge_s + sys->N_edge_v) + (indx - 1) * sys->N_cell_y + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = -1 / (sys->zn[indz] - sys->zn[indz - 1]);
                Senum++;
                leng_Se++;


                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = -1 / (sys->zn[indz] - sys->zn[indz - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = 1 / (sys->xn[indx] - sys->xn[indx - 1]);
                Senum++;

                SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = -1 / (sys->xn[indx] - sys->xn[indx - 1]);
                Senum++;

                SeRowId[Senum] = indz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                SeColId[Senum] = leng_Se;
                Seval[Senum] = 1 / (sys->zn[indz] - sys->zn[indz - 1]);
                Senum++;
                leng_Se++;
            }
        }
    }

    /* the toppest layer doesn't contain the upper plane */
    indz = sys->N_cell_z + 1;
    for (indx = 1; indx <= sys->N_cell_x; indx++){
        for (indy = 1; indy <= sys->N_cell_y; indy++){
            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->xn[indx] - sys->xn[indx - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->xn[indx] - sys->xn[indx - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->yn[indy] - sys->yn[indy - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->yn[indy] - sys->yn[indy - 1]);
            Senum++;
            leng_Se++;
        }
    }

    /* the rightmost yz plane */
    indx = sys->N_cell_x + 1;
    for (indz = 1; indz <= sys->N_cell_z; indz++){
        for (indy = 1; indy <= sys->N_cell_y; indy++){
            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->zn[indz] - sys->zn[indz - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->yn[indy] - sys->yn[indy - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->yn[indy] - sys->yn[indy - 1]);
            Senum++;

            SeRowId[Senum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->zn[indz] - sys->zn[indz - 1]);
            Senum++;
            leng_Se++;
        }
    }

    /* the farthest xz plane */
    indy = sys->N_cell_y + 1;
    for (indz = 1; indz <= sys->N_cell_z; indz++){
        for (indx = 1; indx <= sys->N_cell_x; indx++){
            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->zn[indz] - sys->zn[indz - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->xn[indx] - sys->xn[indx - 1]);
            Senum++;

            SeRowId[Senum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->xn[indx] - sys->xn[indx - 1]);
            Senum++;

            SeRowId[Senum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->zn[indz] - sys->zn[indz - 1]);
            Senum++;
            leng_Se++;
        }
    }
    /* End of the generation of Se */

    /* Generate Sh */
    double *dxa = (double*)malloc((sys->N_cell_x + 1) * sizeof(double));
    double *dya = (double*)malloc((sys->N_cell_y + 1) * sizeof(double));
    double *dza = (double*)malloc((sys->N_cell_z + 1) * sizeof(double));
    dxa[0] = sys->xn[1] - sys->xn[0];
    dya[0] = sys->yn[1] - sys->yn[0];
    dza[0] = sys->zn[1] - sys->zn[0];
    dxa[sys->N_cell_x] = sys->xn[sys->N_cell_x] - sys->xn[sys->N_cell_x - 1];
    dya[sys->N_cell_y] = sys->yn[sys->N_cell_y] - sys->yn[sys->N_cell_y - 1];
    dza[sys->N_cell_z] = sys->zn[sys->N_cell_z] - sys->zn[sys->N_cell_z - 1];
    for (int i = 1; i < sys->N_cell_x; i++){
        dxa[i] = (sys->xn[i + 1] - sys->xn[i - 1]) / 2;
    }
    for (int i = 1; i < sys->N_cell_y; i++){
        dya[i] = (sys->yn[i + 1] - sys->yn[i - 1]) / 2;
    }
    for (int i = 1; i < sys->N_cell_z; i++){
        dza[i] = (sys->zn[i + 1] - sys->zn[i - 1]) / 2;
    }


    /* the middle layers */
    for (indz = 1; indz <= sys->N_cell_z; indz++){
        for (indx = 1; indx <= sys->N_cell_x; indx++){
            for (indy = 1; indy <= sys->N_cell_y; indy++){
                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -1 / (dxa[indx - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = 1 / (dxa[indx]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = 1 / (dya[indy - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -1 / (dya[indy]);
                Shnum++;
                leng_Sh++;


                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = 1 / (dza[indz - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -1 / (dya[indy - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1 + 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = 1 / (dya[indy]);
                Shnum++;

                ShRowId[Shnum] = indz * (sys->N_edge_s + sys->N_edge_v) + (indx - 1) * sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -1 / (dza[indz]);
                Shnum++;
                leng_Sh++;


                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -1 / (dza[indz - 1]);
                Shnum++;


                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = 1 / (dxa[indx - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -1 / (dxa[indx]);
                Shnum++;

                ShRowId[Shnum] = indz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = 1 / (dza[indz]);
                Shnum++;
                leng_Sh++;
            }
        }
    }

    /* the toppest layer doesn't contain the upper plane */
    indz = sys->N_cell_z + 1;
    for (indx = 1; indx <= sys->N_cell_x; indx++){
        for (indy = 1; indy <= sys->N_cell_y; indy++){
            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -1 / (dxa[indx - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = 1 / (dxa[indx]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = 1 / (dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -1 / (dya[indy]);
            Shnum++;
            leng_Sh++;

        }
    }

    /* the rightmost yz plane */
    indx = sys->N_cell_x + 1;
    for (indz = 1; indz <= sys->N_cell_z; indz++){
        for (indy = 1; indy <= sys->N_cell_y; indy++){
            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = 1 / (dza[indz - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -1 / (dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = 1 / (dya[indy]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -1 / (dza[indz]);
            Shnum++;
            leng_Sh++;
        }
    }

    /* the farthest xz plane */
    indy = sys->N_cell_y + 1;
    for (indz = 1; indz <= sys->N_cell_z; indz++){
        for (indx = 1; indx <= sys->N_cell_x; indx++){
            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -1 / (dza[indz - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = 1 / (dxa[indx - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -1 / (dxa[indx]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = 1 / (dza[indz]);
            Shnum++;
            leng_Sh++;
        }
    }
    /* End of the generation of Sh */

    
    /* Multiply Sh and Se */
    myint *ShColIdo, *SeColIdo, *ShRowIdn, *SeRowIdn;
    double *Shvaln, *Sevaln;
    myint ind;
    int status;
    SeColIdo = (myint*)malloc(Senum * sizeof(myint));
    SeRowIdn = (myint*)malloc(Senum * sizeof(myint));
    Sevaln = (double*)malloc(Senum * sizeof(double));
    for (ind = 0; ind < Senum; ind++){
        SeColIdo[ind] = SeColId[ind];
        SeRowIdn[ind] = SeRowId[ind];
        Sevaln[ind] = Seval[ind];
    }
    free(SeRowId); SeRowId = NULL;
    free(Seval); Seval = NULL;
    free(SeColId); SeColId = (myint*)malloc((leng_Se + 1) * sizeof(myint));
    status = COO2CSR_malloc(SeColIdo, SeRowIdn, Sevaln, Senum, leng_Se, SeColId);
    

    ShColIdo = (myint*)malloc(Shnum * sizeof(myint));
    ShRowIdn = (myint*)malloc(Shnum * sizeof(myint));
    Shvaln = (double*)malloc(Shnum * sizeof(double));
    for (ind = 0; ind < Shnum; ind++){
        ShColIdo[ind] = ShColId[ind];
        ShRowIdn[ind] = ShRowId[ind];
        Shvaln[ind] = Shval[ind];
    }
    free(ShRowId); ShRowId = NULL;
    free(Shval); Shval = NULL;
    free(ShColId); ShColId = (myint*)malloc((leng_Sh + 1) * sizeof(myint));
    status = COO2CSR_malloc(ShColIdo, ShRowIdn, Shvaln, Shnum, leng_Sh, ShColId);

    /* generate the matrix (-w^2*D_eps+iw*D_sig+S), alreayd decide the boundary condition as lower boundary PEC */
    status = mklMatrixMulti_nt(sys, sys->leng_S, ShRowIdn, ShColId, Shvaln, sys->N_edge, leng_Se, SeRowIdn, SeColId, Sevaln);    // matrix multiplication first matrix keep the same second matrix transpose
    

    free(dxa); dxa = NULL;
    free(dya); dya = NULL;
    free(dza); dza = NULL;
    free(SeRowId); SeRowId = NULL;
    free(SeColId); SeColId = NULL;
    free(SeColIdo); SeColIdo = NULL;
    free(Seval); Seval = NULL;
    free(ShRowId); ShRowId = NULL;
    free(ShColId); ShColId = NULL;
    free(ShColIdo); ShColIdo = NULL;
    free(Shval); Shval = NULL;

    /*ofstream out;
    out.open("S.txt", std::ofstream::out | std::ofstream::trunc);
    for (int index = 0; index < sys->leng_S; index++){
        out << sys->SRowId[index] << " " << sys->SColId[index] << " " << sys->Sval[index] << endl;
    }
    out.close();*/

    return 0;
}

int reference(fdtdMesh *sys, double freq, int sourcePort, myint *RowId, myint *ColId, double *val){
    int bdn;
#ifdef UPPER_BOUNDARY_PEC
    bdn = 2;
#else
    bdn = 1;
#endif

    myint size = sys->N_edge - bdn * sys->N_edge_s;
    myint *RowId1 = (myint*)malloc((size + 1) * sizeof(myint));
    int count = 0;
    int indi = 0;
    int k = 0;
    complex<double> *valc;
    valc = (complex<double>*)calloc(sys->leng_S, sizeof(complex<double>));
    complex<double> *J;
    J = (complex<double>*)calloc((sys->N_edge - bdn * sys->N_edge_s), sizeof(complex<double>));
    int indz, indy, temp;
    for (indi = 0; indi < sys->portEdge[sourcePort].size(); indi++){
        /* start hard code, add the function of multiple current injections */
        //indz = (sys->portEdge[sourcePort][indi] + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v);
        //indy = (sys->portEdge[sourcePort][indi] % (sys->N_edge_s + sys->N_edge_v) - sys->N_cell_y * (sys->N_cell_x + 1)) % (sys->N_cell_y + 1);
        /* end hard code */
        J[sys->portEdge[sourcePort][indi] - sys->N_edge_s] = 0. - (1i) * sys->portCoor[sourcePort].portDirection * freq * 2. * M_PI;
        /* start hard code */
        //indy++;
        //temp = sys->portEdge[sourcePort][indi] - sys->N_edge_s;
        //while (sys->yn[indy] <= -sys->portCoor[sourcePort].y1*1.01){
        //    temp++;
        //    //cout << sys->yn[indy] << endl;
        //    J[temp] = 0. - (1i) * sys->portCoor[sourcePort].portDirection * freq * 2. * M_PI;
        //    indy++;
        //}
        /* end hard code */
    }

    /* Used in plasma2D for upper and lower excitation */
    /*myint current_edge = sys->portEdge[sourcePort][indi - 1] + (sys->N_edge_s + sys->N_edge_v);
    
    while (current_edge < sys->N_edge - sys->N_edge_s){
        if (sys->markEdge[current_edge] == 0){
            J[current_edge - sys->N_edge_s] = 0. + (1i) * sys->portCoor[sourcePort].portDirection * freq * 2. * M_PI;
        }
        current_edge = current_edge + (sys->N_edge_s + sys->N_edge_v);
    }*/
	/* end of Used in plasma2D for upper and lower excitation */
    
    RowId1[k] = 0;
    k++;
    myint start;
    myint nnz = sys->leng_S;
    //cout << "Start to generate CSR form for S!\n";
    indi = 0;
    while (indi < nnz){
        start = RowId[indi];
        while (indi < nnz && RowId[indi] == start) {
            valc[indi] += val[indi]; // val[indi] is real
            if (RowId[indi] == ColId[indi]){
                if (sys->markEdge[RowId[indi] + sys->N_edge_s] != 0) {
                    complex<double> addedPart(-(2. * M_PI * freq) * sys->stackEpsn[(RowId[indi] + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0, SIGMA * sys->sig[RowId[indi] + sys->N_edge_s]);
                    valc[indi] += (2. * M_PI * freq) * addedPart;
                }
                else {
                    complex<double> addedPart(-(2. * M_PI * freq) * sys->stackEpsn[(RowId[indi] + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0, 0);
                    valc[indi] += (2. * M_PI * freq) * addedPart;
                }
            }
            count++;
            indi++;
        }
        RowId1[k] = (count);
        k++;
    }

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
    //iparm[59] = 2;    // out of core version to solve very large problem
    //iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */

    //cout << "Begin to solve (-w^2*D_eps+iwD_sig+S)x=-iwJ\n";
    complex<double> *xr;
    xr = (complex<double>*)calloc((sys->N_edge - bdn * sys->N_edge_s), sizeof(complex<double>));

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, valc, RowId1, ColId, &perm, &nrhs, iparm, &msglvl, J, xr, &error);
    if (error != 0){
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }

    //cout << "Solving (-w^2*D_eps+iwD_sig+S)x=-iwJ is complete!\n";

    // Just for debug purpose

    int inz, inx, iny;
    double leng;
    for (indi = 0; indi < sys->numPorts; indi++) {
        sys->x.push_back((0, 0));
        for (int j = 0; j < sys->portEdge[indi].size(); j++) {
            if (sys->portEdge[indi][j] % (sys->N_edge_s + sys->N_edge_v) >= sys->N_edge_s) {    // this edge is along z axis
                inz = sys->portEdge[indi][j] / (sys->N_edge_s + sys->N_edge_v);
                leng = sys->zn[inz + 1] - sys->zn[inz];
            }
            else if (sys->portEdge[indi][j] % (sys->N_edge_s + sys->N_edge_v) >= (sys->N_cell_y) * (sys->N_cell_x + 1)) {    // this edge is along x axis
                inx = ((sys->portEdge[indi][j] % (sys->N_edge_s + sys->N_edge_v)) - (sys->N_cell_y) * (sys->N_cell_x + 1)) / (sys->N_cell_y + 1);
                leng = sys->xn[inx + 1] - sys->xn[inx];
            }
            else {    // this edge is along y axis
                iny = (sys->portEdge[indi][j] % (sys->N_edge_s + sys->N_edge_v)) % sys->N_cell_y;
                leng = sys->yn[iny + 1] - sys->yn[iny];
            }

            /*leng = pow((sys->nodepos[sys->edgelink[sys->portEdge[indi][j] * 2] * 3] - sys->nodepos[sys->edgelink[sys->portEdge[indi][j] * 2 + 1] * 3]), 2);
            leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[indi][j] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->portEdge[indi][j] * 2 + 1] * 3 + 1]), 2);
            leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[indi][j] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->portEdge[indi][j] * 2 + 1] * 3 + 2]), 2);
            leng = sqrt(leng);*/
            
            complex<double> addedPart((xr[sys->portEdge[indi][j] - sys->N_edge_s].real() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection))), (xr[sys->portEdge[indi][j] - sys->N_edge_s].imag() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection))));
            
            /* start hard code */
            //complex<double> addedPart((xr[sys->portEdge[indi][j] - sys->N_edge_s].real() * leng / ((1.5e-4 * (sys->zn[indz + 1] - sys->zn[indz - 1]) / 2) * (-sys->portCoor[sourcePort].portDirection))), (xr[sys->portEdge[indi][j] - sys->N_edge_s].imag() * leng / ((1.5e-4 * (sys->zn[indz + 1] - sys->zn[indz - 1]) / 2) * (-sys->portCoor[sourcePort].portDirection))));
            /* end hard code */

            sys->x[sys->x.size() - 1] += addedPart;

        }
    }


    free(xr); xr = NULL;
    free(RowId1); RowId1 = NULL;
    free(valc); valc = NULL;
    return 0;
}


int reference1(fdtdMesh *sys, double freq, int sourcePort, myint *RowId, myint *ColId, double *val, complex<double> *xr){   // to check the accuracy of the result from the inverse model
    int bdn;
#ifdef UPPER_BOUNDARY_PEC
    bdn = 2;
#else
    bdn = 1;
#endif
    myint size = sys->N_edge - bdn * sys->N_edge_s;
    myint *RowId1 = (myint*)malloc((size + 1) * sizeof(myint));
    int count = 0;
    int indi = 0;
    int k = 0;
    complex<double> *valc;
    valc = (complex<double>*)calloc(sys->leng_S, sizeof(complex<double>));
    complex<double> *J;
    J = (complex<double>*)calloc((sys->N_edge - bdn * sys->N_edge_s), sizeof(complex<double>));
    for (indi = 0; indi < sys->portEdge[sourcePort].size(); indi++){
        J[sys->portEdge[sourcePort][indi] - sys->N_edge_s] = 0. - (1i) * 1 * freq * 2. * M_PI;    // portDirection is -1
        
    }

    RowId1[k] = 0;
    k++;
    myint start;
    myint nnz = sys->leng_S;
    //cout << "Start to generate CSR form for S!\n";
    indi = 0;
    while (indi < nnz){
        start = RowId[indi];
        while (indi < nnz && RowId[indi] == start) {
            valc[indi] += val[indi]; // val[indi] is real
            if (RowId[indi] == ColId[indi]){
                if (sys->markEdge[RowId[indi] + sys->N_edge_s] != 0) {
                    complex<double> addedPart(-(2. * M_PI * freq) * sys->stackEpsn[(RowId[indi] + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0, SIGMA);
                    valc[indi] += (2. * M_PI * freq) * addedPart;
                }
                else {
                    complex<double> addedPart(-(2. * M_PI * freq) * sys->stackEpsn[(RowId[indi] + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0, 0);
                    valc[indi] += (2. * M_PI * freq) * addedPart;
                }
            }
            count++;
            indi++;
        }
        RowId1[k] = (count);
        k++;
    }

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

    //cout << "Begin to solve (-w^2*D_eps+iwD_sig+S)x=-iwJ\n";

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, valc, RowId1, ColId, &perm, &nrhs, iparm, &msglvl, J, xr, &error);
    if (error != 0){
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }

    //cout << "Solving (-w^2*D_eps+iwD_sig+S)x=-iwJ is complete!\n";


    free(RowId1); RowId1 = NULL;
    free(valc); valc = NULL;
    return 0;
}

//int V0_reference(fdtdMesh *sys, int sourcePort, double freq) {
//    int bdn;
//#ifdef UPPER_BOUNDARY_PEC
//    bdn = 2;
//#else
//    bdn = 1;
//#endif
//    V0aTDV0 M;
//    myint k = 0;
//    int inz, inx, iny;
//    myint indi;
//    double *dxa, *dya, *dza;
//    dxa = (double*)calloc(sys->nx, sizeof(double));
//    dya = (double*)calloc(sys->ny, sizeof(double));
//    dza = (double*)calloc(sys->nz, sizeof(double));
//    
//    dxa[0] = sys->xn[1] - sys->xn[0];
//    dxa[sys->nx - 1] = sys->xn[sys->nx - 1] - sys->xn[sys->nx - 2];
//    for (indi = 1; indi < sys->nx - 1; indi++) {
//        dxa[indi] = (sys->xn[indi + 1] - sys->xn[indi - 1]) / 2;
//    }
//    dya[0] = sys->yn[1] - sys->yn[0];
//    dya[sys->ny - 1] = sys->yn[sys->ny - 1] - sys->yn[sys->ny - 2];
//    for (indi = 1; indi < sys->ny - 1; indi++) {
//        dya[indi] = (sys->yn[indi + 1] - sys->yn[indi - 1]) / 2;
//    }
//    dza[0] = sys->zn[1] - sys->zn[0];
//    dza[sys->nz - 1] = sys->zn[sys->nz - 1] - sys->zn[sys->nz - 2];
//    for (indi = 1; indi < sys->nz - 1; indi++) {
//        dza[indi] = (sys->zn[indi + 1] - sys->zn[indi - 1]) / 2;
//    }
//
//    myint count = 0;    // count the number of unknowns in M
//    for (indi = sys->N_node_s; indi < sys->N_node - (bdn - 1) * sys->N_node_s; indi++){
//        inz = indi / (sys->N_node_s);
//        inx = (indi % sys->N_node_s) / (sys->N_cell_y + 1);
//        iny = (indi % sys->N_node_s) % (sys->N_cell_y + 1);
//        
//        if (inz == 1 || inz == sys->nz - (bdn)) {
//            count += 1;
//        }
//        else {
//            count += 2;
//        }
//
//        if (inx == 0 || inx == sys->nx - 1) {
//            count += 1;
//        }
//        else {
//            count += 2;
//        }
//
//        if (iny == 0 || iny == sys->ny - 1) {
//            count += 1;
//        }
//        else {
//            count += 2;
//        }
//        count += 1;
//    }
//
//    M.RowId = (myint*)calloc(count, sizeof(myint));
//    M.RowId1 = (myint*)calloc(sys->N_node - bdn * sys->N_node_s + 1, sizeof(myint));
//    M.ColId = (myint*)calloc(count, sizeof(myint));
//    M.val = (complex<double>*)calloc(count, sizeof(complex<double>));
//	M.rhs = (complex<double>*)calloc(sys->N_node - bdn * sys->N_node_s, sizeof(complex<double>));
//
//    indi = 0;
//    M.RowId1[k] = 0;
//    k++;
//    myint start;
//    myint leng = 0;
//	for (indi = sys->N_node_s; indi < sys->N_node - (bdn - 1) * sys->N_node_s; indi++) {    // generate V0a'*D*V0
//		inz = indi / (sys->N_node_s);
//		inx = (indi % sys->N_node_s) / (sys->N_cell_y + 1);
//		iny = (indi % sys->N_node_s) % (sys->N_cell_y + 1);
//
//		if (inz == 1) {
//
//		}
//		else {
//			M.RowId[leng] = indi - sys->N_node_s;
//			M.ColId[leng] = indi - sys->N_node_s * 2;
//			if (sys->markEdge[(inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny] != 0) {
//				M.val[leng] = -1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] = -1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//			leng++;
//			M.rhs[indi - sys->N_egde_s].imag += -1 * freq * 2 * M_PI * (-1 / dza[inz]) * sys->J[(inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny];
//		}
//
//
//
//		if (inx == 0) {
//
//		}
//		else {
//			M.RowId[leng] = indi - sys->N_node_s;
//			M.ColId[leng] = indi - sys->N_node_s - sys->N_cell_y - 1;
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx - 1) + iny] != 0) {
//				M.val[leng] = -1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx - 1) + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] = -1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx - 1) + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//			leng++;
//		}
//
//		if (iny > 0) {
//			M.RowId[leng] = indi - sys->N_node_s;
//			M.ColId[leng] = indi - sys->N_node_s - 1;
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny - 1] != 0) {
//				M.val[leng] = -1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / dyz[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny - 1 + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] = -1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / dya[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny - 1 + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//			leng++;
//		}
//
//		// Row: this node to Col: this node
//		M.RowId[leng] = indi - sys->N_node_s;
//		M.ColId[leng] = indi - sys->N_node_s;
//		if (inz > 0) {
//			if (sys->markEdge[(inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny] != 0) {
//				M.val[leng] += 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += 1 / (sys->zn[inz] - sys->zn[inz - 1]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz - 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//		}
//		if (inz < sys->nz - 1) {
//			if (sys->markEdge[(inz) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny] != 0) {
//				M.val[leng] += 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz + 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += 1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz + 1) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//		}
//		if (inx > 0) {
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx - 1) + iny] != 0) {
//				M.val[leng] += 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx - 1) + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += 1 / (sys->xn[inx] - sys->xn[inx - 1]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx - 1) + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//		}
//		if (inx < sys->nx - 1) {
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx)+iny] != 0) {
//				M.val[leng] += 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx)+iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += 1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx)+iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//		}
//		if (iny > 0) {
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny - 1] != 0) {
//				M.val[leng] += 1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / dyz[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny - 1 + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += 1 / (sys->yn[iny] - sys->yn[iny - 1]) * 1 / dya[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny - 1 + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//		}
//		if (iny < sys->ny - 1) {
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny] != 0) {
//				M.val[leng] += 1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / dyz[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += 1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / dya[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//		}
//
//		leng++;
//
//
//		// Row: this node, Col: y+1 node
//		if (iny < sys->ny - 1) {
//			M.RowId[leng] = indi - sys->N_node_s;
//			M.ColId[leng] = indi - sys->N_node_s + 1;
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny] != 0) {
//				M.val[leng] += -1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / dya[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += -1 / (sys->yn[iny + 1] - sys->yn[iny]) * 1 / dya[iny] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + inx * sys->N_cell_y + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//			leng++;
//		}
//
//		// Row: this node, Col: x+1 node
//		if (inx < sys->nx - 1) {
//			M.RowId[leng] = indi - sys->N_node_s;
//			M.ColId[leng] = indi - sys->N_node_s + sys->N_cell_y + 1;
//			if (sys->markEdge[inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx)+iny] != 0) {
//				M.val[leng] += -1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx)+iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += -1 / (sys->xn[inx + 1] - sys->xn[inx]) * 1 / dxa[inx] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[(inz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (sys->N_cell_y + 1) * (inx)+iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//			leng++;
//		}
//
//		// Row: this node, Col: z+1 node
//		if (inz < sys->nz - bdn) {
//			M.RowId[leng] = indi - sys->N_node_s;
//			M.ColId[leng] = indi - sys->N_node_s + sys->N_node_s;
//			if (sys->markEdge[(inz) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny] != 0) {
//				M.val[leng] += -1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * freq * 2 * M_PI * SIGMA);
//			}
//			else {
//				M.val[leng] += -1 / (sys->zn[inz + 1] - sys->zn[inz]) * 1 / dza[inz] * (-pow(freq * 2 * M_PI, 2) * sys->stackEpsn[((inz) * (sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (sys->N_cell_y + 1) * inx + iny + sys->N_edge_v) % (sys->N_edge_s + sys->N_edge_v)] * EPSILON0
//					+ 1i * 0);
//			}
//			leng++;
//		}
//
//		k++;
//		M.RowId1[k] = leng;
//	}
//
//    return 0;
//}

int plotTime(fdtdMesh *sys, int sourcePort, double *u0d, double *u0c){
    int bdn;
#ifdef UPPER_BOUNDARY_PEC
    bdn = 2;
#else
    bdn = 1;
#endif

    clock_t t1 = clock();
    int t0n = 3;
    double dt = DT;
    double tau = 1.e-11;
    double t0 = t0n * tau;
    myint nt = 2 * t0 / dt;
    double *rsc = (double*)calloc((sys->N_edge - bdn * sys->N_edge_s), sizeof(double));
    double *xr = (double*)calloc((sys->N_edge - bdn * sys->N_edge_s) * 3, sizeof(double));
    myint start;
    myint index;
    double leng;
    double *y0np1 = (double*)calloc(sys->N_edge - bdn * sys->N_edge_s, sizeof(double));
    double *y0n = (double*)calloc(sys->N_edge - bdn * sys->N_edge_s, sizeof(double));
    double *y0nm1 = (double*)calloc(sys->N_edge - bdn * sys->N_edge_s, sizeof(double));
    double *temp = (double*)calloc(sys->N_edge - bdn * sys->N_edge_s, sizeof(double));
    double *temp1 = (double*)calloc(sys->N_edge - bdn * sys->N_edge_s, sizeof(double));
    lapack_complex_double *yh1 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *yh2 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *yh3 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *y_h = (lapack_complex_double*)calloc(sys->N_edge - bdn * sys->N_edge_s, sizeof(lapack_complex_double));
    lapack_complex_double *temp2 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *temp3;
    lapack_complex_double *temp4;
    lapack_complex_double *temp5;
    lapack_complex_double *temp6 = (lapack_complex_double*)calloc((sys->N_edge - bdn * sys->N_edge_s) * sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *temp7;
    lapack_complex_double *m_h;
    int status;
    lapack_int *ipiv;
    lapack_int info;
    lapack_complex_double *y = (lapack_complex_double*)calloc(nt, sizeof(lapack_complex_double));
    double *y0 = (double*)calloc(nt, sizeof(double));
    double *yr = (double*)calloc(nt, sizeof(double));
    double nn;
    int inx, iny, inz;
    for (int ind = 1; ind <= nt; ind++){    // time marching to find the repeated eigenvectors
        for (int inde = 0; inde < sys->portEdge[sourcePort].size(); inde++){
            rsc[sys->portEdge[sourcePort][inde] - sys->N_edge_s] = 2000 * exp(-pow((((dt * ind) - t0) / tau), 2)) + 2000 * (dt * ind - t0) * exp(-pow(((dt * ind - t0) / tau), 2)) * (-2 * (dt * ind - t0) / pow(tau, 2));
        }
        
        //for (myint inde = 0; inde < sys->leng_Vh; inde++){
        //    temp2[inde].real = yh1[inde].real - 2 * yh2[inde].real;
        //    temp2[inde].imag = yh1[inde].imag - 2 * yh2[inde].imag;
        //}
        //temp3 = (lapack_complex_double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(lapack_complex_double));
        //temp4 = (lapack_complex_double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(lapack_complex_double));
        //temp5 = (lapack_complex_double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(lapack_complex_double));
        //status = matrix_multi('N', sys->Vh, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, temp2, sys->leng_Vh, 1, temp3);    // V_re1*(yh(:,1)-2*yh(:,2))
        //status = matrix_multi('N', sys->Vh, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, yh1, sys->leng_Vh, 1, temp4);    // V_re1*(yh(:,1))
        //status = matrix_multi('N', sys->Vh, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, yh2, sys->leng_Vh, 1, temp5);    // V_re1*yh(:,2)
        cout << "Step " << ind << endl;
        for (myint inde = 0; inde < sys->N_edge - bdn * sys->N_edge_s; inde++)
        {
            //y0np1[inde] = 1000 * pow(tau, 2) * (exp(-pow(t0, 2) / pow(tau, 2)) - exp(-pow((dt * (ind + 1) - t0) / tau, 2))) * u0d[inde] + 2000 * ((ind + 1) * dt - t0) * exp(-pow(((ind + 1) * dt - t0) / tau, 2)) * u0c[inde];
            y0n[inde] = 1000 * pow(tau, 2) * (exp(-pow(t0, 2) / pow(tau, 2)) - exp(-pow((dt * (ind)-t0) / tau, 2))) * u0d[inde] + 2000 * ((ind)* dt - t0) * exp(-pow(((ind)* dt - t0) / tau, 2)) * u0c[inde];
            //y0nm1[inde] = 1000 * pow(tau, 2) * (exp(-pow(t0, 2) / pow(tau, 2)) - exp(-pow((dt * (ind - 1) - t0) / tau, 2))) * u0d[inde] + 2000 * ((ind - 1) * dt - t0) * exp(-pow(((ind - 1) * dt - t0) / tau, 2)) * u0c[inde];
        //    temp[inde] = (y0np1[inde] + y0nm1[inde] - 2 * y0n[inde]) * sys->eps[inde + sys->N_edge_s] * 2;    // 2 * D_eps*(y0np1+y0nm1-2*y0n)
        //    temp1[inde] = (y0np1[inde] - y0nm1[inde]) * sys->sig[inde + sys->N_edge_s] * dt;    // dt * D_sig*(y0np1-y0nm1)
        //    for (myint inde1 = 0; inde1 < sys->leng_Vh; inde1++){
        //        temp6[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real = sys->Vh[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde].real * (dt * sys->sig[inde + sys->N_edge_s] + 2 * sys->eps[inde + sys->N_edge_s]);
        //        temp6[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag = sys->Vh[inde1 * (sys->N_edge - 2 * sys->N_edge_s) + inde].imag * (dt * sys->sig[inde + sys->N_edge_s] + 2 * sys->eps[inde + sys->N_edge_s]);
        //    }
            
        }
        
        
        //m_h = (lapack_complex_double*)calloc(sys->leng_Vh * sys->leng_Vh, sizeof(lapack_complex_double));
        //status = matrix_multi('T', sys->Vh, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, temp6, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, m_h);
        //index = 0;
        //while (index < sys->leng_S){
        //    start = sys->SRowId[index];
        //    while (index < sys->leng_S && sys->SRowId[index] == start){
        //        y_h[sys->SRowId[index]].real += sys->Sval[index] * temp5[sys->SColId[index]].real * (-2) * pow(dt, 2);
        //        y_h[sys->SRowId[index]].imag += sys->Sval[index] * temp5[sys->SColId[index]].imag * (-2) * pow(dt, 2);
        //        index++;
        //    }
        //    y_h[start].real += temp4[start].real * sys->sig[start + sys->N_edge_s] * dt + temp3[start].real * sys->eps[start + sys->N_edge_s] * 2 - temp1[start] - temp[start] - rsc[start] * 2 * pow(dt, 2);
        //    y_h[start].imag += temp4[start].imag * sys->sig[start + sys->N_edge_s] * dt + temp3[start].imag * sys->eps[start + sys->N_edge_s] * 2;
        //}
        //status = matrix_multi('T', sys->Vh, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, y_h, (sys->N_edge - 2 * sys->N_edge_s), 1, yh3);
        //ipiv = (lapack_int*)malloc(sys->leng_Vh * sizeof(lapack_int));
        //info = LAPACKE_zgesv(LAPACK_COL_MAJOR, sys->leng_Vh, 1, m_h, sys->leng_Vh, ipiv, yh3, sys->leng_Vh);
        //temp7 = (lapack_complex_double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(lapack_complex_double));
        //status = matrix_multi('N', sys->Vh, (sys->N_edge - 2 * sys->N_edge_s), sys->leng_Vh, yh2, sys->leng_Vh, 1, temp7);    // yhh = V_re1*yh(:,2)

        //for (myint inde = 0; inde < sys->leng_Vh; inde++){
        //    yh1[inde].real = yh2[inde].real;
        //    yh1[inde].imag = yh2[inde].imag;
        //    yh2[inde].real = yh3[inde].real;
        //    yh2[inde].imag = yh3[inde].imag;
        //    yh3[inde].real = 0;
        //    yh3[inde].imag = 0;
        //}
        
        /* central difference */
        index = 0;
        nn = 0;
        while (index < sys->leng_S){
            start = sys->SRowId[index];
            while (index < sys->leng_S && sys->SRowId[index] == start){
                xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + sys->SRowId[index]] += sys->Sval[index] * xr[1 * (sys->N_edge - bdn * sys->N_edge_s) + sys->SColId[index]] * (-2) * pow(dt, 2);
                index++;
            }
           
            if (sys->markEdge[start + sys->N_edge_s] != 0) {
                xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + start] += -rsc[start] * 2 * pow(dt, 2) + dt * SIGMA * xr[start] - 2 * sys->stackEpsn[(start + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * xr[start] + 4 * sys->stackEpsn[(start + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * xr[1 * (sys->N_edge - bdn * sys->N_edge_s) + start];
                xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + start] /= (2 * sys->stackEpsn[(start + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 + dt * SIGMA);
                nn += pow(xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + start], 2);
                
            }
            else {
                xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + start] += -rsc[start] * 2 * pow(dt, 2) - 2 * sys->stackEpsn[(start + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * xr[start] + 4 * sys->stackEpsn[(start + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0 * xr[1 * (sys->N_edge - bdn * sys->N_edge_s) + start];
                xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + start] /= (2 * sys->stackEpsn[(start + sys->N_edge_s + sys->N_edge_v) / (sys->N_edge_s + sys->N_edge_v)] * EPSILON0);
                nn += pow(xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + start], 2);
                
            }
            
        }
        cout << "The norm of xr is " << sqrt(nn) << endl;
        for (myint inde = 0; inde < (sys->N_edge - bdn * sys->N_edge_s); inde++){
            xr[inde] = xr[1 * (sys->N_edge - bdn * sys->N_edge_s) + inde];
            xr[1 * (sys->N_edge - bdn * sys->N_edge_s) + inde] = xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + inde];
            xr[2 * (sys->N_edge - bdn * sys->N_edge_s) + inde] = 0;
        }
        
        for (myint j = 0; j < sys->portEdge[sourcePort].size(); j++){
            
            if (sys->portEdge[sourcePort][j] % (sys->N_edge_s + sys->N_edge_v) >= sys->N_edge_s) {    // this edge is along z axis
                inz = sys->portEdge[sourcePort][j] / (sys->N_edge_s + sys->N_edge_v);
                leng = sys->zn[inz + 1] - sys->zn[inz];
            }
            else if (sys->portEdge[sourcePort][j] % (sys->N_edge_s + sys->N_edge_v) >= (sys->N_cell_y) * (sys->N_cell_x + 1)) {    // this edge is along x axis
                inx = ((sys->portEdge[sourcePort][j] % (sys->N_edge_s + sys->N_edge_v)) - (sys->N_cell_y) * (sys->N_cell_x + 1)) / (sys->N_cell_y + 1);
                leng = sys->xn[inx + 1] - sys->xn[inx];
            }
            else {    // this edge is along y axis
                iny = (sys->portEdge[sourcePort][j] % (sys->N_edge_s + sys->N_edge_v)) % sys->N_cell_y;
                leng = sys->yn[iny + 1] - sys->yn[iny];
            }
            
            //y[ind - 1].real += (temp7[sys->portEdge[sourcePort][j] - sys->N_edge_s].real * leng) + y0n[sys->portEdge[sourcePort][j] - sys->N_edge_s] * leng;
            //y[ind - 1].imag += (temp7[sys->portEdge[sourcePort][j] - sys->N_edge_s].imag * leng);
            y0[ind - 1] += y0n[sys->portEdge[sourcePort][j] - sys->N_edge_s] * leng;
            yr[ind - 1] += xr[sys->portEdge[sourcePort][j] - sys->N_edge_s] * leng;
        }
        
        //free(temp3);temp3 =NULL;
        //free(temp4); temp4 = NULL;
        //free(temp5); temp5 = NULL;
        //free(temp7); temp7 = NULL;
        //free(ipiv); ipiv = NULL;
        //free(m_h); m_h = NULL;
        
    }
    
    free(temp); temp = NULL;
    free(temp1); temp1 = NULL;
    free(temp2); temp2 = NULL;
    free(temp6); temp6 = NULL;
    free(y0n); y0n = NULL;
    free(y0np1); y0np1 = NULL;
    free(y0nm1); y0nm1 = NULL;
    free(yh1); yh1 = NULL;
    free(yh2); yh2 = NULL;
    free(yh3); yh3 = NULL;
    free(xr); xr = NULL;
    free(rsc); rsc = NULL;
    ofstream out;
    /*out.open("y1.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint ind = 0; ind < nt; ind++){
        out << y[ind].real << " " << y[ind].imag << endl;
    }
    out.close();*/

    out.open("y01.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint ind = 0; ind < nt; ind++){
        out << y0[ind] << endl;
    }
    out.close();

    out.open("yr1.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint ind = 0; ind < nt; ind++){
        out << yr[ind] << endl;
    }
    out.close();

    return 0;
}

int mklMatrixMulti_nt(fdtdMesh *sys, myint &leng_A, myint *aRowId, myint *aColId, double *aval, myint arow, myint acol, myint *bRowId, myint *bColId, double *bval){
    // ArowId, AcolId, and Aval should be in the COO format
    sparse_status_t s0;
    sparse_matrix_t a, a_csr;
    sparse_index_base_t indexing1 = SPARSE_INDEX_BASE_ZERO;
    myint row, col;
    myint *cols, *cole, *rowi;
    double *val;
    MKL_INT *AcolId;
    double *Aval;
    myint k;
    s0 = mkl_sparse_d_create_csr(&a, SPARSE_INDEX_BASE_ZERO, acol, arow, &aColId[0], &aColId[1], aRowId, aval);
    
    sparse_matrix_t b, b_csr;
    s0 = mkl_sparse_d_create_csr(&b, SPARSE_INDEX_BASE_ZERO, acol, arow, &bColId[0], &bColId[1], bRowId, bval);
    
    sparse_matrix_t A;
    s0 = mkl_sparse_spmm(SPARSE_OPERATION_TRANSPOSE, a, b, &A);
    sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;
    myint ARows, ACols;
    MKL_INT *ArowStart, *ArowEnd;
        
    s0 = mkl_sparse_d_export_csr(A, &indexing, &ARows, &ACols, &ArowStart, &ArowEnd, &AcolId, &Aval);
    leng_A = ArowEnd[ARows - 1];    // how many non-zeros are in S

    sys->SRowId = (myint*)malloc(leng_A * sizeof(myint));
    sys->SColId = (myint*)malloc(leng_A * sizeof(myint));
    sys->Sval = (double*)calloc(leng_A, sizeof(double));
    myint count, num, j;
    j = 0;
    vector<pair<myint, double>> col_val;
    for (myint i = 0; i < leng_A; i++){
        col_val.push_back(make_pair(AcolId[i], Aval[i]));
    }
    for (myint i = 0; i < ARows; i++){
        if (i < sys->N_edge_s){
            continue;
        }
#ifdef UPPER_BOUNDARY_PEC
        if (i >= sys->N_edge - sys->N_edge_s){
            continue;
        }
#else

#endif
        num = ArowEnd[i] - ArowStart[i];
        count = 0;
        vector<pair<myint, double>> v(col_val.begin() + ArowStart[i], col_val.begin() + ArowEnd[i]);
        sort(v.begin(), v.end());
        while (count < num){
            if (v[count].first < sys->N_edge_s){
                count++;
                continue;
            }
#ifdef UPPER_BOUNDARY_PEC
            if (v[count].first >= sys->N_edge - sys->N_edge_s){
                count++;
                continue;
            }
#else

#endif
            
            sys->SRowId[j] = i - sys->N_edge_s;
            sys->SColId[j] = v[count].first - sys->N_edge_s;
            sys->Sval[j] = v[count].second / MU;
            j++;
            count++;
        }
        v.clear();
    }
    leng_A = j;

    mkl_sparse_destroy(a);
    mkl_sparse_destroy(b);
    mkl_sparse_destroy(A);
    return 0;
}

int matrix_multi_cd(char operation, lapack_complex_double *a, myint arow, myint acol, double *b, myint brow, myint bcol, lapack_complex_double *tmp3){    //complex multiply double
    /* operation = 'T' is first matrix conjugate transpose, operation = 'N' is first matrix non-conjugate-transpose*/
    if (operation == 'T'){
        for (myint ind = 0; ind < acol; ind++){
            for (myint ind1 = 0; ind1 < bcol; ind1++){
                for (myint ind2 = 0; ind2 < arow; ind2++){
                    tmp3[ind1 * acol + ind].real = tmp3[ind1 * acol + ind].real + a[ind * arow + ind2].real * b[ind1 * brow + ind2];
                    tmp3[ind1 * acol + ind].imag = tmp3[ind1 * acol + ind].imag - a[ind * arow + ind2].imag * b[ind1 * brow + ind2];
                }
            }
        }
    }
    else if (operation == 'N'){
        for (myint ind = 0; ind < arow; ind++){
            for (myint ind1 = 0; ind1 < bcol; ind1++){
                for (myint ind2 = 0; ind2 < acol; ind2++){
                    tmp3[ind1 * arow + ind].real = tmp3[ind1 * arow + ind].real + a[ind2 * arow + ind].real * b[ind1 * brow + ind2];
                    tmp3[ind1 * arow + ind].imag = tmp3[ind1 * arow + ind].imag + a[ind2 * arow + ind].imag * b[ind1 * brow + ind2];
                }
            }
        }
    }
    return 0;
}
