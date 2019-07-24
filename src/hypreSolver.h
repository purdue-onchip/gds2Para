#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"


//int hypreSolve(fdtdMesh *sys, HYPRE_IJMatrix A, HYPRE_ParCSRMatrix parcsr_A, myint leng_A, double *bin, myint leng_v0, double *solution);
int hypreSolve(fdtdMesh *sys, myint *ARowId, myint *AColId, double *Aval, myint leng_A, double *bin, myint leng_v0, double *solution);
int hypre_FlexGMRESModifyPCAMGExample(void *precond_data, HYPRE_Int iterations, double rel_residual_norm);