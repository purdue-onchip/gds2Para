/* Generate the stiffness matrix */
#include "fdtd.h"



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
    cout << "Length of S is " << sys->leng_S << endl;
    ofstream outfile;
    outfile.open("S.txt", std::ofstream::out | std::ofstream::trunc);
    for (ind = 0; ind < sys->leng_S; ind++){
        outfile << sys->SRowId[ind] + 1 << " " << sys->SColId[ind] + 1 << " " << sys->Sval[ind] << endl;
    }
    outfile.close();
    cout << "S generation is done!\n";
    /*  */


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
    //    
    //    /*ofstream outfile;
    //    outfile.open("xr.txt", std::ofstream::out | std::ofstream::trunc);
    //    for (int i = 0; i < sys->N_edge - 2 * sys->N_edge_s; i++){
    //        outfile << xr[i].re << " " << xr[i].i << endl;
    //    }
    //    outfile.close();*/

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

//int pardisoSolve_r(fdtdMesh *sys, complex<double> *rhs, int *RowId, int *ColId, complex<double> *val, int nnz, int size, complex<double> *solution){
//    
//    int *RowId1 = (int*)malloc((size + 1) * sizeof(int));
//    int count = 0;
//    int i = 0;
//    int k = 0;
//    int start;
//    RowId1[k] = 0;
//    k++;
//    while (i < nnz){
//        start = RowId[i];
//        while (i < nnz && RowId[i] == start) {
//            count++;
//            i++;
//        }
//        RowId1[k] = (count);
//        k++;
//    }
//
//
//    
//    int mtype = 13;    /* Real complex unsymmetric matrix */
//    int nrhs = 1;    /* Number of right hand sides */
//    void *pt[64];
//
//    /* Pardiso control parameters */
//    int iparm[64];
//    int maxfct, mnum, phase, error, msglvl, solver;
//    double dparm[64];
//    int v0csin;
//    int perm;
//
//    /* Number of processors */
//    int num_procs;
//
//    /* Auxiliary variables */
//    char *var;
//
//    msglvl = 1;    /* print statistical information */
//    solver = 0;    /* use sparse direct solver */
//    error = 0;
//    maxfct = 1;
//    mnum = 1;
//    phase = 13;
//
//    pardisoinit(pt, &mtype, iparm);
//    iparm[38] = 1;
//    iparm[34] = 1;    // 0-based indexing
//    //iparm[10] = 0;        /* Use nonsymmetric permutation and scaling MPS */
//    iparm[23] = 1;        /// !PARDISO uses new two - level factorization algorithm
//    iparm[12] = 1;    // a maximum weighted matching algorithm 
//
//    pardiso(pt, &maxfct, &mnum, &mtype, &phase, &size, val, RowId1, ColId, &perm, &nrhs, iparm, &msglvl, rhs, solution, &error);
//    if (error != 0){
//        printf("\nERROR during numerical factorization: %d", error);
//        exit(2);
//    }
//
//    i = 0;   // row number
//    int j = 0;
//    int nz = RowId1[size];
//    double numeratorr = 0, numerator = 0;
//    complex<double> temp;
//    
//    /*while (i < size){
//        temp = 0 + 1i * 0;
//        while (RowId[j] == i && j < nz){
//            temp += val[j] * solution[ColId[j]];
//            j++;
//        }
//        temp = temp.real() - rhs[i].real() + 1i * (temp.imag() - rhs[i].imag());
//        numeratorr += pow(temp.real(), 2) + pow(temp.imag(), 2); 
//        i++;
//    }
//    numeratorr = sqrt(numeratorr);
//
//    i = 0;
//    while (i < size){
//        temp = 0 + 1i * 0;
//        while (RowId[j] == i && j < nz){
//            temp += val[j] * sys->y[ColId[j] + sys->N_edge_s];
//            j++;
//        }
//        temp = temp.real() - rhs[i].real() + 1i * (temp.imag() - rhs[i].imag());
//        numerator += pow(temp.real(), 2) + pow(temp.imag(), 2);
//        i++;
//    }
//    numerator = sqrt(numerator);*/
//    while (i < size){
//        temp = sys->y[i + sys->N_edge_s].real() - solution[i].real() + 1i * (sys->y[i + sys->N_edge_s].imag() - solution[i].imag());
//        numerator += pow(temp.real(), 2) + pow(temp.imag(), 2);
//        i++;
//    }
//    numerator = sqrt(numerator);
//    
//    i = 0;
//    double denominator = 0;
//    while (i < size){
//        denominator += pow(rhs[i].real(), 2) + pow(rhs[i].imag(), 2);
//        
//        i++;
//    }
//    cout << endl;
//    cout << denominator << endl;
//    denominator = sqrt(denominator);
//
//    /*cout << "Relative residual of xr is " << numeratorr / denominator << endl;*/
//    cout << "Relative residual  is " << numerator / denominator << endl;
//    cout << "numerator is " << numerator << " and denominator is " << denominator << endl;
//
//    free(RowId1); RowId1 = NULL;
//
//}

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
    leng_A = ArowEnd[ARows - 1];

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