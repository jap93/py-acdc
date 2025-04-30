import sys

import numpy as np

from nequip.ase.nequip_calculator import NequIPCalculator
from ase import Atoms

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

        self.device = "cuda"

        self.atoms = Atoms()
        self.NequIPCalc = None #NequIPCalculator()

        self.first_setup = True

    def finalise(self):
        pass

    def readPotential(self, fieldFileName, spec: Species, out_stream):
        #read the species and potentials
        try:
            with open(fieldFileName, 'r') as in_stream:

                while True:
                    line = in_stream.readline()
                    if not line:
                        break
            
                    #out_stream.write(f" input line: {line.strip()}\n")
                    if line[0] == '#':
                        continue

                    words = line.split()

                    if not words:
                        continue

                    if words[0].lower() == "close":
                        return
            
                    elif words[0] == "species":
                        num = int(words[1])
                        self.species = Species()
                        self.species.load_species(in_stream, num)

                    elif words[0].lower() == "device":
                        self.device = words[1].lower()

                    elif words[0].lower() == "nequip":
                        self.model = words[1]
        except FileNotFoundError:
            #the c++ needs to get a copy of the potential parameters  
            out_stream.write("\n*** could not find potentials file\n")
            sys.exit(1)

    def get_calculator(self):
        return self.NequIPCalc     
        
    def setup(self, basin: Config):

        if self.first_setup == False:
            return
        
        if self.model is None:
            print("a model name is required")
            exit(-1)
        
        self.NequIPCalc = NequIPCalculator.from_deployed_model(
            model_path=self.model,
            device=self.device,)

        
        #self.atoms.set_pbc((True, True, True))

        self.first_setup = False

    def calculate_forces(self, pos_r: np.ndarray, force: np.ndarray, stress: np.ndarray, lat_vector: np.ndarray, rcp_vector: np.ndarray, 
                         charge: np.ndarray, atm_label, symbls, frozen, natoms, total_energy: Energy):
        
        atoms = Atoms(positions=pos_r, cell=lat_vector, symbols=symbls, pbc=True)
        
        atoms.set_calculator(self.NequIPCalc)
        #self.atoms.set_chemical_symbols(symbls)

        total_energy.vdwEnergy = atoms.get_potential_energy()
        force = atoms.get_forces()
        stress = atoms.get_stress()
        #print("total energy", atoms.get_potential_energy())
