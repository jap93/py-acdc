import sys

import numpy as np

from janus_core.helpers.mlip_calculators import choose_calculator
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

        #self.atoms = Atoms()
        self.janCalc = None #NequIPCalculator()

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

                    elif words[0].lower() == "model":
                        self.model = words[1]
        except FileNotFoundError:
            #the c++ needs to get a copy of the potential parameters  
            out_stream.write("\n*** could not find potentials file\n")
            sys.exit(1)

    def get_calculator(self):
        return self.janCalc
        
    def setup(self, basin: Config):

        if self.first_setup == False:
            return
        
        if self.model is None:
            print("a model name is required")
            exit(-1)
        
        try:
            self.janCalc = choose_calculator(architecture="mace_mp", model=self.model, precision="float64", device=self.device)

        except Exception as e:
            print(f"{e}\n")
            exit()

        self.first_setup = False

    def calculate_forces(self, pos_r: np.ndarray, force: np.ndarray, stress: np.ndarray, lat_vector: np.ndarray, rcp_vector: np.ndarray, 
                         charge: np.ndarray, atm_label, symbls, frozen, natoms, total_energy: Energy):
        
        atoms = Atoms(positions=pos_r, cell=lat_vector, symbols=symbls, pbc=True)
        
        atoms.set_calculator(self.janCalc)

        total_energy.vdwEnergy = atoms.get_potential_energy()
        np.copyto(force, atoms.get_forces())

        force[frozen == 1] = 0.0 # if an atom is frozen then set force to zero
        np.copyto(stress, atoms.get_stress())
        #print("total energy", atoms.get_potential_energy())

    
