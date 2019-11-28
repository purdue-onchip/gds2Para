#include <cmath>
#include <ctime>
#include "fdtd.hpp"
#include "hypreSolver.h"
#include "vis.c"


//int hypreSolve(fdtdMesh *sys, HYPRE_IJMatrix A, HYPRE_ParCSRMatrix parcsr_A, myint leng_A, double *bin, myint leng_v0, double *solution) {
int hypreSolve(fdtdMesh *sys, myint *ARowId, myint *AColId, double *Aval, myint leng_A, double *bin, myint leng_v0, double *solution, int info_show, int solve_id) {
	/* ARowId : matrix rowId
	   AColId : matrix colId
	   Aval : matrix value
	   leng_A : nnz in A
	   bin : right hand side
	   leng_v0 : matrix size, row number or column number
	   solution : the result 
       info_show : whether show the information of iteration number and relative residual, 0 not show, 1 show
       solve : which algorithm is implemented */
    /* Initialize MPI */
    int num_procs, myid;

    // Get the number of processors
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    /* Initialize HYPRE for multiple processes */
    HYPRE_Int indi = 0;
    int N = leng_v0;
    HYPRE_Int local_size = N / num_procs; // Integer division
    HYPRE_Int extra = N - local_size * num_procs; // Remaining length of matrix row not evenly given to processes

    HYPRE_Int ilower = local_size * myid; // Index of matrix row where this process starts
    ilower += hypre_min(myid, extra);
    HYPRE_Int iupper = local_size * (myid + 1); // Index of matrix row where this process ends
    iupper += hypre_min(myid + 1, extra);
    iupper -= 1;

    local_size = iupper - ilower + 1; // Corrected number of rows in process
    //cout << "iupper is " << iupper << " and ilower is " << ilower << endl;
    /* Create the square matrix. [Square matrix => indicate the
    row partition size twice (since N_rows = N_cols)] */
    HYPRE_IJMatrix A;
    HYPRE_ParCSRMatrix parcsr_A;
    HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &A);

    /* Choose a parallel CSR format storage (see the User's Manual) */
    HYPRE_IJMatrixSetObjectType(A, HYPRE_PARCSR);

    /* Initialize before setting coefficients matrix values */
    HYPRE_IJMatrixInitialize(A);

    vector<double> values;
    vector<HYPRE_Int> cols;
    HYPRE_Int index = 0;
    for (indi = ilower; indi <= iupper; indi++) {
        HYPRE_Int nnz = 0;   // Number of non-zeros on row indi
        while (index < leng_A && ARowId[index] == indi) {
            cols.push_back(AColId[index]);
            values.push_back(Aval[index]);
            nnz++;
            index++;
        }

        /* Set the values for row indi */
        HYPRE_IJMatrixSetValues(A, 1, &nnz, &indi, &cols[0], &values[0]);
        cols.clear();
        values.clear();
    }

    /* Assemble after setting the coefficients matrix */
    HYPRE_IJMatrixAssemble(A);

    /* Get the parcsr matrix object to use */
    HYPRE_IJMatrixGetObject(A, (void**)&parcsr_A);

    /* Create the rhs vector and solution vector */
    HYPRE_IJVector b;
    HYPRE_ParVector par_b;
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &b);
    HYPRE_IJVectorSetObjectType(b, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(b);

    HYPRE_IJVector x;
    HYPRE_ParVector par_x;
    HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper, &x);
    HYPRE_IJVectorSetObjectType(x, HYPRE_PARCSR);
    HYPRE_IJVectorInitialize(x);

    /* Set the RHS values to vector bin and the initial solution to zero vector */
    double *rhs_values = (double*)calloc(local_size, sizeof(double));
    double *x_values = (double*)calloc(local_size, sizeof(double));
    HYPRE_Int *rows = (HYPRE_Int*)calloc(local_size, sizeof(HYPRE_Int));

    for (indi = 0; indi < local_size; indi++) {
        rhs_values[indi] = bin[indi];
        x_values[indi] = 0.0;
        rows[indi] = ilower + indi;
    }
    HYPRE_IJVectorSetValues(b, local_size, rows, rhs_values);
    HYPRE_IJVectorSetValues(x, local_size, rows, x_values);

    /* Release heap memory for the two vectors */
    free(x_values);
    free(rhs_values);
    free(rows);

    /* Assemble and get par vector objects after setting the RHS vector and solution vector*/
    HYPRE_IJVectorAssemble(b);
    HYPRE_IJVectorGetObject(b, (void **)&par_b);
    HYPRE_IJVectorAssemble(x);
    HYPRE_IJVectorGetObject(x, (void **)&par_x);

    /* Call the appropriate preconditioner and solver */
//#if HYPRE_METHOD == 1
    /* AMG */
    if (solve_id == 1) {
        int num_iterations;
        double final_res_norm;

        /* Create solver */
        HYPRE_Solver solver;
        HYPRE_BoomerAMGCreate(&solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_BoomerAMGSetPrintLevel(solver, 0);           /* no printout */
        HYPRE_BoomerAMGSetOldDefault(solver);              /* Falgout coarsening with modified classical interpolaiton */
        HYPRE_BoomerAMGSetRelaxType(solver, 4);            /* G-S/Jacobi hybrid relaxation */
        HYPRE_BoomerAMGSetRelaxOrder(solver, 1);           /* uses C/F relaxation */
        HYPRE_BoomerAMGSetNumSweeps(solver, 1);            /* Sweeps on each level */
        HYPRE_BoomerAMGSetMaxLevels(solver, 20);           /* maximum number of levels */
        HYPRE_BoomerAMGSetTol(solver, HYPRE_CONV_TOL);     /* conv. tolerance */
        HYPRE_BoomerAMGSetMaxIter(solver, HYPRE_MAX_ITER); /* maximum number of iterations */

        /* Now setup and solve! */
        HYPRE_BoomerAMGSetup(solver, parcsr_A, par_b, par_x);
        HYPRE_BoomerAMGSolve(solver, parcsr_A, par_b, par_x);

        /* Run info - needed logging turned on */
        HYPRE_BoomerAMGGetNumIterations(solver, &num_iterations);
        HYPRE_BoomerAMGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (myid == 0 && info_show == 1)
        {
            cout << endl << " Iterations = " << num_iterations << endl;
            cout << " Final Relative Residual Norm = " << final_res_norm << endl;
        }

        /* Output the final solution */
        HYPRE_Complex v;
        for (indi = ilower; indi <= iupper; indi++) {
            HYPRE_IJVectorGetValues(x, 1, &indi, &v);
            solution[indi] = v;
        }

        /* Destroy solver */
        HYPRE_BoomerAMGDestroy(solver);
    }
//#elif HYPRE_METHOD == 2
    /* PCG with AMG preconditioner */
    else if (solve_id == 2) {
        int num_iterations;
        double final_res_norm;

        /* Create solver and preconditioner */
        HYPRE_Solver solver, precond;
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);
        HYPRE_BoomerAMGCreate(&precond);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, HYPRE_MAX_ITER); /* max iterations */
        HYPRE_PCGSetTol(solver, HYPRE_CONV_TOL);     /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1);              /* use the 2-norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 0);           /* print solve info */
        HYPRE_PCGSetLogging(solver, 1);              /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print AMG solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 6);
        HYPRE_BoomerAMGSetOldDefault(precond);
        HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero for preconditioner */
        HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */

        /* Set the AMG preconditioner for the PCG solver */
        HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
        HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);

        /* Run info - needed logging turned on */
        HYPRE_PCGGetNumIterations(solver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (myid == 0 && info_show == 1)
        {
            cout << endl << " Iterations = " << num_iterations << endl;
            cout << " Final Relative Residual Norm = " << final_res_norm << endl;
        }

        /* Output the final solution */
        HYPRE_Complex v;
        for (indi = ilower; indi <= iupper; indi++) {
            HYPRE_IJVectorGetValues(x, 1, &indi, &v);
            solution[indi] = v;
        }

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
//#elif HYPRE_METHOD == 3
    /* Flexible GMRES with AMG Preconditioner */
    else if (solve_id == 3) {
        HYPRE_Int    num_iterations;
        double final_res_norm;
        HYPRE_Int    restart = 10;
        bool modify = true;

        /* Create solver */
        HYPRE_Solver solver, precond;
        HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);
        HYPRE_BoomerAMGCreate(&precond);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_FlexGMRESSetKDim(solver, restart);
        HYPRE_FlexGMRESSetMaxIter(solver, HYPRE_MAX_ITER); /* max iterations */
        HYPRE_FlexGMRESSetTol(solver, HYPRE_CONV_TOL);     /* conv. tolerance */
        HYPRE_FlexGMRESSetPrintLevel(solver, 0);           /* print solve info */
        HYPRE_FlexGMRESSetLogging(solver, 1);              /* needed to get run info later */

        /* Now set up the AMG preconditioner and specify any parameters */
        HYPRE_BoomerAMGSetPrintLevel(precond, 0); /* print AMG solution info */
        HYPRE_BoomerAMGSetCoarsenType(precond, 6);
        HYPRE_BoomerAMGSetOldDefault(precond);
        HYPRE_BoomerAMGSetRelaxType(precond, 6); /* Sym G.S./Jacobi hybrid */
        HYPRE_BoomerAMGSetNumSweeps(precond, 1);
        HYPRE_BoomerAMGSetTol(precond, 0.0); /* conv. tolerance zero for preconditioner */
        HYPRE_BoomerAMGSetMaxIter(precond, 1); /* do only one iteration! */
        //HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);
        //HYPRE_EuclidSetMem(precond, 1);

        /* Set the AMG preconditioner for the Flexible GMRES solver */
        HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSolve,
            (HYPRE_PtrToSolverFcn)HYPRE_BoomerAMGSetup, precond);
        //HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_EuclidSolve, (HYPRE_PtrToSolverFcn)HYPRE_EuclidSetup, precond);

        /* Modify AMG parameters at runtime if needed and allowed */
        if (modify) {
            /* Optional call: if not called, hypre_FlexGMRESModifyPCDefault is used, which does nothing.
            Otherwise, the custom defined hypre_FlexGMRESModifyPCAMG() is called. */
            HYPRE_FlexGMRESSetModifyPC(solver, (HYPRE_PtrToModifyPCFcn)hypre_FlexGMRESModifyPCAMG);
        }

        /* Now setup and solve! */
        //cout << "Setup and solve in GMRES block" << endl;
        HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_A, par_b, par_x);
        //cout << "GMRES has been set up, awaiting solve" << endl;
        HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_A, par_b, par_x);

        /* Run info - needed logging turned on */
        HYPRE_FlexGMRESGetNumIterations(solver, &num_iterations);
        HYPRE_FlexGMRESGetFinalRelativeResidualNorm(solver, &final_res_norm);

        if (myid == 0 && info_show == 1)
        {
            cout << endl << " Iterations = " << num_iterations << endl;
            cout << " Final Relative Residual Norm = " << final_res_norm << endl;
        }

        /* Output the final solution */
        HYPRE_Complex v;
        for (indi = ilower; indi <= iupper; indi++) {
            HYPRE_IJVectorGetValues(x, 1, &indi, &v);
            solution[indi] = v;
        }

        /* Destroy solver and preconditioner */
        HYPRE_ParCSRFlexGMRESDestroy(solver);
        HYPRE_BoomerAMGDestroy(precond);
    }
    /* PCG with Parasails Preconditioner */
    else if (solve_id == 4) {
        HYPRE_Int num_iterations;
        double final_res_norm;

        int      sai_max_levels = 1;
        double   sai_threshold = 0.1;
        double   sai_filter = 0.05;
        int      sai_sym = 1;

        /* Create solver */
        HYPRE_Solver solver, precond;
        HYPRE_ParCSRPCGCreate(MPI_COMM_WORLD, &solver);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_PCGSetMaxIter(solver, 100); /* max iterations */
        HYPRE_PCGSetTol(solver, 1e-5); /* conv. tolerance */
        HYPRE_PCGSetTwoNorm(solver, 1); /* use the two norm as the stopping criteria */
        HYPRE_PCGSetPrintLevel(solver, 2); /* print solve info */
        HYPRE_PCGSetLogging(solver, 1); /* needed to get run info later */

        /* Now set up the ParaSails preconditioner and specify any parameters */
        HYPRE_ParaSailsCreate(MPI_COMM_WORLD, &precond);

        /* Set some parameters (See Reference Manual for more parameters) */
        HYPRE_ParaSailsSetParams(precond, sai_threshold, sai_max_levels);
        HYPRE_ParaSailsSetFilter(precond, sai_filter);
        HYPRE_ParaSailsSetSym(precond, sai_sym);
        HYPRE_ParaSailsSetLogging(precond, 3);

        /* Set the PCG preconditioner */
        HYPRE_PCGSetPrecond(solver, (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSolve,
            (HYPRE_PtrToSolverFcn)HYPRE_ParaSailsSetup, precond);

        /* Now setup and solve! */
        HYPRE_ParCSRPCGSetup(solver, parcsr_A, par_b, par_x);
        HYPRE_ParCSRPCGSolve(solver, parcsr_A, par_b, par_x);


        /* Run info - needed logging turned on */
        HYPRE_PCGGetNumIterations(solver, &num_iterations);
        HYPRE_PCGGetFinalRelativeResidualNorm(solver, &final_res_norm);
        if (myid == 0)
        {
            printf("\n");
            printf("Iterations = %lld\n", num_iterations);
            printf("Final Relative Residual Norm = %e\n", final_res_norm);
            printf("\n");
        }

        /* Destory solver and preconditioner */
        HYPRE_ParCSRPCGDestroy(solver);
        HYPRE_ParaSailsDestroy(precond);
    }

//#endif

    /* Clean up */
    HYPRE_IJMatrixDestroy(A);
    HYPRE_IJVectorDestroy(b);
    HYPRE_IJVectorDestroy(x);
    return(0);
}

int hypre_FlexGMRESModifyPCAMG(void *precond_data, HYPRE_Int iterations, double rel_residual_norm) {
    if (rel_residual_norm > HYPRE_PC_TOL) {
        HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data, HYPRE_PC_MOD_SWEEPS);
    }
    else {
        HYPRE_BoomerAMGSetNumSweeps((HYPRE_Solver)precond_data, 1);
    }
    return 0;
}
