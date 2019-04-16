#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


int hypreSolve(fdtdMesh *sys, HYPRE_IJMatrix A, HYPRE_ParCSRMatrix parcsr_A, myint leng_A, double *bin, myint leng_v0, double *solution);