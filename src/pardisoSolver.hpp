#ifndef GDS2PARA_PARDISO_SOLVER_H_
#define GDS2PARA_PARDISO_SOLVER_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

// Cal all the computed freq points and store in a vector
vector<double> CalAllFreqPointsHz(const fdtdMesh &sys) {
	vector<double> vFreqHz;

	vFreqHz.push_back(sys.freqStart * sys.freqUnit);          // First frequency in sweep

	for (int id = 1; id < sys.nfreq; id++) {                  // When nfreq > 1
		if (sys.freqScale == 1) {
			vFreqHz.push_back((sys.freqStart + id * (sys.freqEnd - sys.freqStart) / (sys.nfreq - 1)) * sys.freqUnit);
		}   // Linear interpolation of frequency sweep
		else {
			vFreqHz.push_back(sys.freqStart * sys.freqUnit * pow(sys.freqEnd / sys.freqStart, (id * 1.0 / (sys.nfreq - 1))));
		}   // Logarithmic interpolation of frequency sweep
	}

	return vFreqHz;
}

// Assign source current density for a source port
void AssignSourceCurrentForSourcePort(fdtdMesh *psys, int sourcePort) {
	if (psys->J != nullptr) {
		free(psys->J);
	}	// delete previous J assignment

	psys->J = (double*)calloc(psys->N_edge, sizeof(double));
	for (int sourcePortSide = 0; sourcePortSide < psys->portCoor[sourcePort].multiplicity; sourcePortSide++) {
		for (int indEdge = 0; indEdge < psys->portCoor[sourcePort].portEdge[sourcePortSide].size(); indEdge++) {
			/* Set current density for all edges within sides in port to prepare solver */
			/*cout << " port #" << sourcePort + 1 << ", side #" << sourcePortSide + 1 
				<< ", edge #" << sys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge] 
				<< ": J = " << sys->portCoor[sourcePort].portDirection[sourcePortSide] << " A/m^2" << endl;*/
			psys->J[psys->portCoor[sourcePort].portEdge[sourcePortSide][indEdge]] = psys->portCoor[sourcePort].portDirection[sourcePortSide];
		}
	}
}

// Solve E field and Z-parameters in Pardiso, solve layer by layer. (under developing)
void Solve_E_Zpara_InPardiso_layered(fdtdMesh *psys) {
    int bdl = 0, bdu = 0;
#ifdef UPPER_BOUNDARY_PEC
    bdu = 1;
#endif
#ifdef LOWER_BOUNDARY_PEC
    bdl = 1;
#endif

	// All computed freq points
	vector<double> vFreqHz = CalAllFreqPointsHz(*psys);

    // S matrix in CSR format
    const myint *pcsrS_colId = psys->SColId;

	// Right-Hand-Side term -iwJ (A * m^-2 / s) and electric field eField (V/m), SI unit
	complex<double> *pRhsJ_SI = new complex<double>[psys->N_edge - 2 * psys->N_edge_s]();    // -iw{j}
	complex<double> *peField_SI = new complex<double>[psys->N_edge - 2 * psys->N_edge_s]();  // {e}
    
    // Use complex double constructor to assign initial output matrix for single-frequency solve
    psys->x.assign(psys->numPorts * psys->numPorts * psys->nfreq, complex<double>(0., 0.));

    for (int indFreq = 0; indFreq < vFreqHz.size(); indFreq++) {
        //for (int sourcePort = 0; sourcePort < psys->numPorts; sourcePort++) {
        //    AssignSourceCurrentForSourcePort(psys, sourcePort);
        //    Solve_E_InPardiso(psys, peField_SI, pRhsJ_SI, vFreqHz[indFreq], psys->SRowId, psys->SColId, psys->Sval);
        //    psys->Construct_Z_V0_Vh(peField_SI, indFreq, sourcePort, bdl, bdu);
        //}   // solve for each port
        reference(psys, indFreq, psys->SRowId, psys->SColId, psys->Sval);
        
	}   // solve for each frequency

    // Print Z-parameters
    psys->print_z_V0_Vh();

	delete[] pRhsJ_SI, peField_SI;
}

// Solve E field and Z-parameters in Pardiso, solve the whole structure as reference
void Solve_E_Zpara_InPardiso_reference(fdtdMesh *psys) {
    int bdl = 0, bdu = 0;
#ifdef UPPER_BOUNDARY_PEC
    bdu = 1;
#endif
#ifdef LOWER_BOUNDARY_PEC
    bdl = 1;
#endif

    // All computed freq points
    vector<double> vFreqHz = CalAllFreqPointsHz(*psys);

    // Use complex double constructor to assign initial output matrix for single-frequency solve
    psys->x.assign(psys->numPorts * psys->numPorts * psys->nfreq, complex<double>(0., 0.));

    // At each frequency, solve all ports together
    for (int indFreq = 0; indFreq < vFreqHz.size(); indFreq++) {
        reference(psys, indFreq, psys->SRowId, psys->SColId, psys->Sval);
    }   

    // Print Z-parameters
    psys->print_z_V0_Vh();

}

#endif