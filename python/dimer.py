import os
import numpy as np
from scipy.linalg import eigh_tridiagonal

from config import Config
from energy import Energy
from search_control import SearchCntrl
from field import Field
from status import Status
from job_control import JobControl

class Dimer(object):

    def __init__(self):
        pass

    def _initialise(self, frozen, natoms, num_frozen, indeces, mask):

        for i in range(natoms):
            if frozen[i] == 1:
                indeces.append(i)
                mask.append(False)

    """
    def _initialise(self, frozen, natoms, indeces, mask):
    indeces.extend(i for i in range(natoms) if frozen[i] == 1)
    mask.extend([False] * len(indeces))
    """
    def run(self, saddle: Config, r_vector, displacement_vector, saddle_energy: Energy, fld: Field, job: JobControl, saddle_status: Status, search_peprpindicular: bool, out_stream):

        saddle_status.set_status_failed()

        ase_search = False
        flag = (os.environ.get("USE_ASE_SEARCH", ase_search) == "True")
        
        if flag == True:  # use ase dimer method
            from ase import Atoms
            from ase.constraints import FixAtoms
            from ase.mep import DimerControl, MinModeAtoms, MinModeTranslate

            basin1 = Atoms(positions=saddle.pos, cell=saddle.vectors, symbols=saddle.symbol, pbc=saddle.pbc)
            #fld.get_calculator(basin1)
            basin1.set_calculator(fld.get_calculator())

            try:
                indeces = []
                mask = []
                self._initialise(saddle.frozen, saddle.natoms, saddle.num_frozen, indeces, mask)

                #fixt atom positions given list of fixed atoms
                c = FixAtoms(indeces)
                #print("fixed atoms ", indices)
                basin1.set_constraint(c)

                # Set up the dimer
                d_control = DimerControl(initial_eigenmode_method='displacement',
                                         displacement_method='vector', logfile=None, mask=mask)
                d_atoms = MinModeAtoms(basin1, d_control)

                # Displace the atoms
                #print("disp vec ",displacement_vector)
                d_atoms.displace(displacement_vector=displacement_vector.reshape((-1,3)))

                # Converge to a saddle point
                with MinModeTranslate(d_atoms, trajectory='dimer_method.traj',
                                      logfile=None) as dim_rlx:
                    new_flag = dim_rlx.run(fmax=job.search_parameters.force_tol, steps=job.search_parameters.max_par_iter)
            
            except Exception as e:
                out_stream.write("\n exception occured in ASE dimer method : ")
                out_stream.write(f"{e}\n")
                exit(-1)
                
            new_basin = d_atoms.atoms
            saddle_energy.vdwEnergy = d_atoms.get_potential_energy()
            saddle.pos = new_basin.get_positions()
            saddle.vectors = new_basin.cell[:]

            if new_flag:
                saddle_status.set_status_success()
                out_stream.write(f"\n saddle energy {saddle_energy.vdwEnergy} \n")
            else:
                out_stream.write("\n dimer method failed")
    
        else:
            raise(NotImplementedError)

    
