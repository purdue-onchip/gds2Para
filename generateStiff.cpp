/* Generate the stiffness matrix */
#include "fdtd.h"
using namespace std::complex_literals;


int generateStiff(fdtdMesh *sys){
    myint Senum, leng_Se;    // Se's size is (N_patch - 2 * N_patch_s) * (N_edge - 2 * N_edge_s)
    myint Shnum, leng_Sh;    // Sh's size is (N_edge - 2 * N_edge_s) * (N_patch - 2 * N_patch_s)
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
    for (indz = 1; indz < sys->N_cell_z; indz++){
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
    indz = sys->N_cell_z;
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

            SeRowId[Senum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = -1 / (sys->xn[indx] - sys->xn[indx - 1]);
            Senum++;

            SeRowId[Senum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->xn[indx] - sys->xn[indx - 1]);
            Senum++;

            SeRowId[Senum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            SeColId[Senum] = leng_Se;
            Seval[Senum] = 1 / (sys->yn[indy] - sys->yn[indy - 1]);
            Senum++;

            SeRowId[Senum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
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
    for (indz = 1; indz < sys->N_cell_z; indz++){
        for (indx = 1; indx <= sys->N_cell_x; indx++){
            for (indy = 1; indy <= sys->N_cell_y; indy++){
                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -dza[indz - 1] / (dxa[indx - 1] * dza[indz - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = dza[indz - 1] / (dxa[indx] * dza[indz - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = dza[indz - 1] / (dza[indz - 1] * dya[indy - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -dza[indz - 1] / (dza[indz - 1] * dya[indy]);
                Shnum++;
                leng_Sh++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = dxa[indx - 1] / (dxa[indx - 1] * dza[indz - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -dxa[indx - 1] / (dxa[indx - 1] * dya[indy - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1 + 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = dxa[indx - 1] / (dxa[indx - 1] * dya[indy]);
                Shnum++;

                ShRowId[Shnum] = indz * (sys->N_edge_s + sys->N_edge_v) + (indx - 1) * sys->N_cell_y + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -dxa[indx - 1] / (dxa[indx - 1] * dza[indz]);
                Shnum++;
                leng_Sh++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -dya[indy - 1] / (dza[indz - 1] * dya[indy - 1]);
                Shnum++;


                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = dya[indy - 1] / (dxa[indx - 1] * dya[indy - 1]);
                Shnum++;

                ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = -dya[indy - 1] / (dya[indy - 1] * dxa[indx]);
                Shnum++;

                ShRowId[Shnum] = indz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
                ShColId[Shnum] = leng_Sh;
                Shval[Shnum] = dya[indy - 1] / (dya[indy - 1] * dza[indz]);
                Shnum++;
                leng_Sh++;
            }
        }
    }

    /* the toppest layer doesn't contain the upper plane */
    indz = sys->N_cell_z;
    for (indx = 1; indx <= sys->N_cell_x; indx++){
        for (indy = 1; indy <= sys->N_cell_y; indy++){
            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dza[indz - 1] / (dza[indz - 1] * dxa[indx - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dza[indz - 1] / (dza[indz - 1] * dxa[indx]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dza[indz - 1] / (dza[indz - 1] * dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dza[indz - 1] / (dza[indz - 1] * dya[indy]);
            Shnum++;
            leng_Sh++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dxa[indx - 1] / (dxa[indx - 1] * dza[indz - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dxa[indx - 1] / (dxa[indx - 1] * dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1 + 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dxa[indx - 1] / (dxa[indx - 1] * dya[indy]);
            Shnum++;

            ShRowId[Shnum] = indz * (sys->N_edge_s + sys->N_edge_v) + (indx - 1) * sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dxa[indx - 1] / (dxa[indx - 1] * dza[indz]);
            Shnum++;
            leng_Sh++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dya[indy - 1] / (dza[indz - 1] * dya[indy - 1]);
            Shnum++;


            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dya[indy - 1] / (dxa[indx - 1] * dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dya[indy - 1] / (dya[indy - 1] * dxa[indx]);
            Shnum++;

            ShRowId[Shnum] = indz * (sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dya[indy - 1] / (dya[indy - 1] * dza[indz]);
            Shnum++;
            leng_Sh++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dza[indz - 1] / (dza[indz - 1] * dxa[indx - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + indx*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dza[indz - 1] / (dza[indz - 1] * dxa[indx]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dza[indz - 1] / (dza[indz - 1] * dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y * (sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dza[indz - 1] / (dza[indz - 1] * dya[indy]);
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
            Shval[Shnum] = dxa[indx - 1] / (dxa[indx - 1] * dza[indz - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dxa[indx - 1] / (dxa[indx - 1] * dya[indy - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy + 1 - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dxa[indx - 1] / (dxa[indx - 1] * dya[indy]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + (indx - 1)*sys->N_cell_y + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dxa[indx - 1] / (dxa[indx - 1] * dza[indz]);
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
            Shval[Shnum] = -dya[indy - 1] / (dya[indy - 1] * dza[indz - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dya[indy - 1] / (dya[indy - 1] * dxa[indx - 1]);
            Shnum++;

            ShRowId[Shnum] = (indz - 1)*(sys->N_edge_s + sys->N_edge_v) + sys->N_edge_s + indx*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = -dya[indy - 1] / (dya[indy - 1] * dxa[indx]);
            Shnum++;

            ShRowId[Shnum] = (indz)*(sys->N_edge_s + sys->N_edge_v) + sys->N_cell_y*(sys->N_cell_x + 1) + (indx - 1)*(sys->N_cell_y + 1) + indy - 1;
            ShColId[Shnum] = leng_Sh;
            Shval[Shnum] = dya[indy - 1] / (dya[indy - 1] * dza[indz]);
            Shnum++;
            leng_Sh++;
        }
    }
    /* End of the generation of Sh */
    /*for (int i = 0; i < Shnum; i++){
        cout << ShRowId[i] << " " << ShColId[i] << " " << Shval[i] << endl;
    }*/
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


    /* generate the matrix (-w^2*D_eps+iw*D_sig+S) */
    status = mklMatrixMulti_nt(sys, sys->leng_S, ShRowIdn, ShColId, Shvaln, sys->N_edge, leng_Se, SeRowIdn, SeColId, Sevaln);    // matrix multiplication first matrix keep the same second matrix transpose
  
    //cout << "Length of S is " << sys->leng_S << endl;
    //cout << "S generation is done!\n";

    /* use pardiso to solve (-w^2*D_eps+iw*D_sig+S)\(-1i*w*J) */
    //int sourcePort = 0;
    //double leng;
    //complex<double> *J;
    //complex<double> *Yr = (complex<double>*)malloc(sys->numPorts * sys->numPorts * sizeof(complex<double>));
    //while (sourcePort < sys->numPorts){
    //    J = (complex<double>*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(complex<double>));
    //    complex<double> *xr = (complex<double>*)calloc((sys->N_edge - 2 * sys->N_edge_s), sizeof(complex<double>));
    //    for (int i = 0; i < sys->portEdge[sourcePort].size(); i++){
    //        J[sys->portEdge[sourcePort][i] - sys->N_edge_s] = -1i * (2 * PI*sys->freqStart * sys->freqUnit) * sys->portCoor[sourcePort].portDirection;
    //    }
    //    
    //    status = pardisoSolve_r(sys, J, sys->SRowId, sys->SColId, sys->Sval, sys->leng_S, sys->N_edge - 2 * sys->N_edge_s, xr);
    //    

    //    for (int i = 0; i < sys->numPorts; i++){
    //        Yr[i + sys->numPorts * sourcePort] = 0 + 1i * 0;
    //        for (int j = 0; j < sys->portEdge[i].size(); j++){
    //            leng = pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3]), 2);
    //            leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3 + 1]), 2);
    //            leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->portEdge[i][j] * 2 + 1] * 3 + 2]), 2);
    //            leng = sqrt(leng);
    //            Yr[i + sys->numPorts * sourcePort] = Yr[i + sys->numPorts * sourcePort].real() + xr[sys->portEdge[i][j] - sys->N_edge_s].real() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection)) + 1i * (xr[sys->portEdge[i][j] - sys->N_edge_s].imag() * leng / (sys->portArea[sourcePort] * (-sys->portCoor[sourcePort].portDirection)) + Yr[i + sys->numPorts * sourcePort].imag());

    //        }
    //        cout << Yr[i + sys->numPorts * sourcePort] << endl;
    //    }

    //    sourcePort++;
    //    free(J); J = NULL;
    //    free(xr); xr = NULL;
    //}

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

    /*for (int i = 0; i < sys->numPorts; i++){
        for (int j = 0; j < sys->numPorts; j++){
            cout << Yr[i * sys->numPorts + j].real() << " " << Yr[i * sys->numPorts + j].imag() << endl;
        }
    }*/

    return 0;
}

int reference(fdtdMesh *sys, lapack_complex_double *x, myint *RowId, myint *ColId, double *val){
    myint size = sys->N_edge - 2 * sys->N_edge_s;
    myint *RowId1 = (myint*)malloc((size + 1) * sizeof(myint));
    int count = 0;
    int indi = 0;
    int k = 0;
    complex<double> *valc;
    valc = (complex<double>*)calloc(sys->leng_S, sizeof(complex<double>));
    complex<double> *J;
    J = (complex<double>*)malloc((sys->N_edge - 2 * sys->N_edge_s) * sizeof(complex<double>));
    for (indi = sys->N_edge_s; indi < sys->N_edge - sys->N_edge_s; indi++){
        J[indi - sys->N_edge_s] = 0. + (1i) * -sys->J[indi] * sys->freqEnd * sys->freqUnit * 2 * PI;
        //J[i - sys->N_edge_s].i = -sys->J[i] * sys->freqEnd * sys->freqUnit * 2 * PI;
        //J[i - sys->N_edge_s].re = 0;
    }
    
    RowId1[k] = 0;
    k++;
    myint start;
    myint nnz = sys->leng_S;
    cout << "Start to generate CSR form for S!\n";
    indi = 0;
    while (indi < nnz){
        start = RowId[indi];
        while (indi < nnz && RowId[indi] == start) {
            valc[indi] += val[indi]; // val[indi] is real
            if (RowId[indi] == ColId[indi]){
                // valc[indi] needs omega * (-omega * epsilon  + j * sigma) added to it
                complex<double> addedPart(-(2 * PI * sys->freqEnd * sys->freqUnit) * sys->eps[RowId[indi] + sys->N_edge_s], sys->sig[RowId[indi] + sys->N_edge_s]);
                valc[indi] += (2 * PI * sys->freqEnd * sys->freqUnit) * addedPart;
                //valc[indi] += -std::complex::pow((2 * PI*sys->freqEnd * sys->freqUnit), 2) * sys->eps[RowId[indi] + sys->N_edge_s] + (1i) * (2 * PI*sys->freqEnd * sys->freqUnit) * sys->sig[RowId[indi] + sys->N_edge_s];
                //valc[i].re += - pow((2 * PI*sys->freqEnd * sys->freqUnit), 2) * sys->eps[RowId[i] + sys->N_edge_s];
                //valc[i].i += (2 * PI*sys->freqEnd * sys->freqUnit) * sys->sig[RowId[i] + sys->N_edge_s];
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
    
    cout << "Begin to solve (-w^2*D_eps+iwD_sig+S)x=-iwJ\n";
    complex<double> *xr;
    xr = (complex<double>*)calloc((sys->N_edge - 2 * sys->N_edge_s), sizeof(complex<double>));

    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, valc, RowId1, ColId, &perm, &nrhs, iparm, &msglvl, J, xr, &error);
    if (error != 0){
        printf("\nERROR during numerical factorization: %d", error);
        exit(2);
    }

    cout << "Solving (-w^2*D_eps+iwD_sig+S)x=-iwJ is complete!\n";
    double nn = 0.;
    double nnn = 0.;
    for (indi = 0; indi < sys->N_edge - 2 * sys->N_edge_s; indi++){
        nn += pow((xr[indi].real() - x[indi].real), 2) + pow((xr[indi].imag() - x[indi].imag), 2);
        nnn += pow(xr[indi].real(), 2) + pow(xr[indi].imag(), 2);
        //nn += pow((xr[i].re - x[i].real), 2) + pow((xr[i].i - x[i].imag), 2);
        //nnn += pow(xr[i].re, 2) + pow(xr[i].i, 2);
    }
    ofstream out;
    /* out.open("x.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_edge - sys->N_edge_s * 2; i++){
        out << x[i].real << " " << x[i].imag << endl;
    }
    out.close();*/

    /*out.open("xr.txt", std::ofstream::out | std::ofstream::trunc);
    for (i = 0; i < sys->N_edge - sys->N_edge_s * 2; i++){
        out << xr[i].re << " " << xr[i].i << endl;
    }
    out.close();*/

    nn = sqrt(nn);
    nnn = sqrt(nnn);
    /*cout << "Relative residual of xr is " << numeratorr / denominator << endl;*/
    cout << "Relative residual  is " << nn / nnn << endl;

    free(RowId1); RowId1 = NULL;
    free(valc); valc = NULL;

}

int plotTime(fdtdMesh *sys, int sourcePort, double *u0d, double *u0c){

    clock_t t1 = clock();
    int t0n = 3;
    double dt = DT;
    double tau = 1.e-11;
    double t0 = t0n * tau;
    myint nt = 2 * t0 / dt;
    double *rsc = (double*)calloc((sys->N_edge - 2 * sys->N_edge_s), sizeof(double));
    double *xr = (double*)calloc((sys->N_edge - 2 * sys->N_edge_s) * 3, sizeof(double));
    myint start;
    myint index;
    double leng;
    double *y0np1 = (double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(double));
    double *y0n = (double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(double));
    double *y0nm1 = (double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(double));
    double *temp = (double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(double));
    double *temp1 = (double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(double));
    lapack_complex_double *yh1 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *yh2 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *yh3 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *y_h = (lapack_complex_double*)calloc(sys->N_edge - 2 * sys->N_edge_s, sizeof(lapack_complex_double));
    lapack_complex_double *temp2 = (lapack_complex_double*)calloc(sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *temp3;
    lapack_complex_double *temp4;
    lapack_complex_double *temp5;
    lapack_complex_double *temp6 = (lapack_complex_double*)calloc((sys->N_edge - 2 * sys->N_edge_s) * sys->leng_Vh, sizeof(lapack_complex_double));
    lapack_complex_double *temp7;
    lapack_complex_double *m_h;
    int status;
    lapack_int *ipiv;
    lapack_int info;
    lapack_complex_double *y = (lapack_complex_double*)calloc(nt, sizeof(lapack_complex_double));
    double *y0 = (double*)calloc(nt, sizeof(double));
    double *yr = (double*)calloc(nt, sizeof(double));
    double nn;
    for (int ind = 1; ind <= 10; ind++){    // time marching to find the repeated eigenvectors
        /*for (int inde = 0; inde < sys->portEdge[sourcePort].size(); inde++){
            rsc[sys->portEdge[sourcePort][inde] - sys->N_edge_s] = 2000 * exp(-pow((((dt * ind) - t0) / tau), 2)) + 2000 * (dt * ind - t0) * exp(-pow(((dt * ind - t0) / tau), 2)) * (-2 * (dt * ind - t0) / pow(tau, 2));
        }*/
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
        //
        for (myint inde = 0; inde < sys->N_edge - 2 * sys->N_edge_s; inde++)
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
        //
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
        
        // central difference
        //index = 0;
        ////nn = 0;
        //while (index < sys->leng_S){
        //    start = sys->SRowId[index];
        //    while (index < sys->leng_S && sys->SRowId[index] == start){
        //        xr[2 * (sys->N_edge - 2 * sys->N_edge_s) + sys->SRowId[index]] += sys->Sval[index] * xr[1 * (sys->N_edge - 2 * sys->N_edge_s) + sys->SColId[index]] * (-2) * pow(dt, 2);
        //        index++;
        //    }
        //    xr[2 * (sys->N_edge - 2 * sys->N_edge_s) + start] += -rsc[start] * 2 * pow(dt, 2) + dt * sys->sig[start + sys->N_edge_s] * xr[start] - 2 * sys->eps[start + sys->N_edge_s] * xr[start] + 4 * sys->eps[start + sys->N_edge_s] * xr[1 * (sys->N_edge - 2 * sys->N_edge_s) + start];
        //    xr[2 * (sys->N_edge - 2 * sys->N_edge_s) + start] /= (2 * sys->eps[start + sys->N_edge_s] + dt * sys->sig[start + sys->N_edge_s]);
        //    //nn += pow(xr[2 * (sys->N_edge - 2 * sys->N_edge_s) + start], 2);
        //}
        ////cout << "The norm of xr is " << sqrt(nn) << endl;
        //for (myint inde = 0; inde < (sys->N_edge - 2 * sys->N_edge_s); inde++){
        //    xr[inde] = xr[1 * (sys->N_edge - 2 * sys->N_edge_s) + inde];
        //    xr[1 * (sys->N_edge - 2 * sys->N_edge_s) + inde] = xr[2 * (sys->N_edge - 2 * sys->N_edge_s) + inde];
        //    xr[2 * (sys->N_edge - 2 * sys->N_edge_s) + inde] = 0;
        //}
        
        /*for (myint j = 0; j < sys->portEdge[sourcePort].size(); j++){
            
            leng = pow((sys->nodepos[sys->edgelink[sys->portEdge[sourcePort][j] * 2] * 3] - sys->nodepos[sys->edgelink[sys->portEdge[sourcePort][j] * 2 + 1] * 3]), 2);
            leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[sourcePort][j] * 2] * 3 + 1] - sys->nodepos[sys->edgelink[sys->portEdge[sourcePort][j] * 2 + 1] * 3 + 1]), 2);
            leng = leng + pow((sys->nodepos[sys->edgelink[sys->portEdge[sourcePort][j] * 2] * 3 + 2] - sys->nodepos[sys->edgelink[sys->portEdge[sourcePort][j] * 2 + 1] * 3 + 2]), 2);
            leng = sqrt(leng);
            y[ind - 1].real += (temp7[sys->portEdge[sourcePort][j] - sys->N_edge_s].real * leng) + y0n[sys->portEdge[sourcePort][j] - sys->N_edge_s] * leng;
            y[ind - 1].imag += (temp7[sys->portEdge[sourcePort][j] - sys->N_edge_s].imag * leng);
            y0[ind - 1] += y0n[sys->portEdge[sourcePort][j] - sys->N_edge_s] * leng;
            yr[ind - 1] += xr[sys->portEdge[sourcePort][j] - sys->N_edge_s] * leng;
        }
        
        free(temp3);temp3 =NULL;
        free(temp4); temp4 = NULL;
        free(temp5); temp5 = NULL;
        free(temp7); temp7 = NULL;
        free(ipiv); ipiv = NULL;
        free(m_h); m_h = NULL;*/
        
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
    cout << "Time for 10 steps our proposed time marching is " << (clock() - t1) * 1.0 / CLOCKS_PER_SEC << endl;
    /*ofstream out;
    out.open("y1.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint ind = 0; ind < nt; ind++){
        out << y[ind].real << " " << y[ind].imag << endl;
    }
    out.close();

    out.open("y01.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint ind = 0; ind < nt; ind++){
        out << y0[ind] << endl;
    }
    out.close();

    out.open("yr1.txt", std::ofstream::out | std::ofstream::trunc);
    for (myint ind = 0; ind < nt; ind++){
        out << yr[ind] << endl;
    }
    out.close();*/

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
        if (i < sys->N_edge_s || i >= sys->N_edge - sys->N_edge_s){
            continue;
        }
        num = ArowEnd[i] - ArowStart[i];
        count = 0;
        vector<pair<myint, double>> v(col_val.begin() + ArowStart[i], col_val.begin() + ArowEnd[i]);
        sort(v.begin(), v.end());
        while (count < num){
            if (v[count].first < sys->N_edge_s || v[count].first >= sys->N_edge - sys->N_edge_s){
                count++;
                continue;
            }
            sys->SRowId[j] = i - sys->N_edge_s;
            sys->SColId[j] = v[count].first - sys->N_edge_s;
            sys->Sval[j] = v[count].second / MU;
            //if (sys->SRowId[j] == sys->SColId[j]){
            //    sys->Sval[j] = sys->Sval[j].real();// -pow((2 * PI*sys->freqStart * sys->freqUnit), 2) * sys->eps[i + sys->N_edge_s] + 1i * ((2 * PI*sys->freqStart * sys->freqUnit) * sys->sig[i + sys->N_edge_s] + sys->Sval[j].imag());
            //}
            j++;
            count++;
        }
        v.clear();
    }
    leng_A = j;
    

    mkl_sparse_destroy(a);
    mkl_sparse_destroy(b);
    mkl_sparse_destroy(A);
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