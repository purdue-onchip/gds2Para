#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


int hypreSolve(fdtdMesh *sys, HYPRE_IJMatrix A, HYPRE_ParCSRMatrix parcsr_A, int leng_A, double *bin, int leng_v0, double *solution);