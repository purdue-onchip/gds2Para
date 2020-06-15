#include "_hypre_utilities.h"
#include "HYPRE_krylov.h"
#include "HYPRE.h"
#include "HYPRE_parcsr_ls.h"

// HYPRE control macros
#define HYPRE_METHOD (3) // 1 = AMG, 2 = PCG with AMG Preconditioner, 3 = Flexible GMRES with AMG Preconditioner
#define HYPRE_CONV_TOL (1.e-5) // Convergence relative tolerance for HYPRE
#define HYPRE_PC_TOL (1000. * HYPRE_CONV_TOL) // Preconditioner tolerance (relative residual norm) for HYPRE
#define HYPRE_PC_MOD_SWEEPS (10) // Modified number of preconditioner sweeps if preconditioner tolerance not met on first preconditioner sweep for HYPRE
#define HYPRE_MAX_ITER (100) // Maximum iterations for HYPRE

//int hypreSolve(fdtdMesh *sys, HYPRE_IJMatrix A, HYPRE_ParCSRMatrix parcsr_A, myint leng_A, double *bin, myint leng_v0, double *solution);
int hypreSolve(myint *ARowId, myint *AColId, double *Aval, myint leng_A, double *bin, myint leng_v0, double *solution, int info_show, int solve);
int hypre_FlexGMRESModifyPCAMG(void *precond_data, HYPRE_Int iterations, double rel_residual_norm);