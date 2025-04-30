import sys

import numpy as np

import metal

from energy import Energy
from species import Species
from config import Config

class Field:

    def __init__(self):
        
        self.model = None
        self.struc = None
        self.cutoff = 0.0
        self.species = None
        #self.atoms = None

        #self.device = "gpu"

        self.first_setup = True

        self.fld = metal.Metal("potentials")

    def readPotential(self, fieldFileName, spec: Species, out_stream):
        #read the species and potentials
        self.fld.readPotential(fieldFileName)

        out_stream.write(f"finished reading external fieled from {fieldFileName} \n")
        out_stream.flush()    
        
    def get_calculator(self):
        pass

    def finalise(self):
        self.fld.finalise()
        
    def setup(self, basin: Config):

        if self.first_setup == False:
            return
        
        cell_prop = basin.cell_properties(basin.vectors)
        rcpvec = np.linalg.inv(basin.vectors)
        rcp_prop = basin.cell_properties(rcpvec)
        min_dimension = np.min(cell_prop[0:3])
        
        self.fld.setup(cell_prop, rcpvec, rcp_prop, min_dimension, basin.get_volume(), basin.natoms)
        #void setup(const Eigen::VectorXd &, const Eigen::VectorXd &, double, double, int)
        self.first_setup = False

    def calculate_forces(self, pos_r: np.ndarray, force: np.ndarray, stress: np.ndarray, lat_vector: np.ndarray, rcp_vector: np.ndarray, 
                         charge: np.ndarray, atm_label, symbls, frozen, natoms, total_energy: Energy):
        
        #the call to c++ is
        #matrix = calculateForces(const Eigen::MatrixXd &pos, Eigen::MatrixXd &stress, const Eigen::MatrixXd &latVector, 
        #const Eigen::MatrixXd &rcpVector, const Eigen::VectorXi &atmLabel, const Eigen::VectorXi &frozen,
        #                    int numAtoms)
        new_data = self.fld.calculateForces(pos_r, stress, lat_vector, rcp_vector, atm_label, frozen, natoms)
        
        total_energy.manyEnergy = new_data[natoms][0]
        total_energy.vdwEnergy = new_data[natoms][1]
        #print("python energies ", total_energy.rcpEnergy, total_energy.realEnergy, total_energy.vdwEnergy)

        for i in range(natoms):
            force[i][:] = new_data[i][:]

        force[frozen == 1] = 0.0
        #print("force", force)
        #exit(-1)

    
