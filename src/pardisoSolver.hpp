#ifndef GDS2PARA_PARDISO_SOLVER_H_
#define GDS2PARA_PARDISO_SOLVER_H_

#include <iostream> 
#include <fstream> 
#include <string>

#include "fdtd.hpp"

using namespace std;

vector<myint> Map_eInd_GrowZ2Y(const myint Nx, const myint Ny, const myint Nz) {
    /* This function changes the global {e} index from layer growth along Z to that along Y.

    The index follows y - x - z ordering, and {e} is stacked layer by layer as
    {e_surface, e_volumn}.T at each layer. For the two cases here :
    case 1 ~ grow along Z: {e} = {e_s, | e_v}.T = {ey,ex, | ez}.T, in which vector
    {ey} first runs through y for each x, then runs through x as by
    y - x - z ordering. Same for {ex} ~y->x and {ez} ~y->x.
    case 2 ~ grow along Y: {e} = {ex,ez, | ey}.T at every layer.
    {ex}, {ez} and {ey} ~x->z, frist all x for each z then all z

    Input: num of bricks Nx*Ny*Nz
    Output: vector "eInd_map_z2y" of size N_e mapping global indices of {e}, eInd_map_z2y[oldInd] = newInd

    Yee's grid is used, with E at edge center and H at face center. Outmost
    boundaries are all E edges. Removal of {e} due to PEC BC has not been considered. */

    // Num of surface or volumetric unknown e at each layer
    myint n_surfEy_growZ        = Ny*(Nx + 1);
    myint n_surfEyEx            = n_surfEy_growZ + Nx*(Ny + 1);
    myint n_volEz               = (Nx + 1)*(Ny + 1);
    myint n_layerE_growZ        = n_surfEyEx + n_volEz;

    myint n_surfEx_growY        = Nx*(Nz + 1);
    myint n_surfExEz            = n_surfEx_growY + Nz*(Nx + 1);
    myint n_volEy               = (Nx + 1)*(Nz + 1);
    myint n_layerE_growY        = n_surfExEz + n_volEy;

    myint N_tot_E               = n_surfEyEx*(Nz + 1) + n_volEz*Nz;


    // y - x - z ordering
    // Grow along Z : {e} = { ey,ex, | ez }.T, y->x frist all y for each x then all x
    // Grow along Y : {e} = { ex,ez, | ey }.T, x->z frist all x for each z then all z

    vector<myint> eInd_map_z2y(N_tot_E, 0);

    // Map the edges along y direction, Ey
    for (myint iz = 0; iz < Nz + 1; iz++) {
        for (myint iy = 0; iy < Ny; iy++) {
            for (myint ix = 0; ix < Nx + 1; ix++) {
                eInd_map_z2y[iz*n_layerE_growZ + ix*(Ny)+iy] =
                    iy*n_layerE_growY + n_surfExEz + iz*(Nx + 1) + ix;
            }
        }
    }

    // Map the edges along x direction, Ex
    for (myint iz = 0; iz < Nz + 1; iz++) {
        for (myint iy = 0; iy < Ny + 1; iy++) {
            for (myint ix = 0; ix < Nx; ix++) {
                eInd_map_z2y[iz*n_layerE_growZ + n_surfEy_growZ + ix*(Ny + 1) + iy] =
                    iy*n_layerE_growY + iz*(Nx)+ix;
            }
        }
    }

    // Map the edges along z direction, Ez
    for (myint iz = 0; iz < Nz; iz++) {
        for (myint iy = 0; iy < Ny + 1; iy++) {
            for (myint ix = 0; ix < Nx + 1; ix++) {
                eInd_map_z2y[iz*n_layerE_growZ + n_surfEyEx + ix*(Ny + 1) + iy] =
                    iy*n_layerE_growY + n_surfEx_growY + iz*(Nx + 1) + ix;
            }
        }
    }

    return eInd_map_z2y;
}

// Reverse the map by swapping the index and value of the mapping array
vector<myint> Reverse_Map_eInd(const vector<myint> &eInd_map_z2y) {
    
    myint N_tot_E = eInd_map_z2y.size();
    vector<myint> eInd_map_y2z(N_tot_E, 0);

    for (myint id = 0; id < N_tot_E; id++) {
        eInd_map_y2z[eInd_map_z2y[id]] = id;
    }

    return eInd_map_y2z;
}

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