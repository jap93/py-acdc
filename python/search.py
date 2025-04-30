import os
import sys

import numpy as np

from mpi4py import MPI

from config import Config
from field import Field
from species import Species
from job_control import JobControl
from KMCControl import KMCCntrl
from search_control import SearchCntrl
from energy import Energy
from status import Status
from ART import ART
from dimer import Dimer
from relax import Relax


def get_random_number():
    return np.random.random()

class Saddle(object):
    
    def __init__(self, kmc_iteration, event_number, method):
        self.vector_size = 0
        self.r_vector = None
        self.delta_pos = None

        self.kmc_iteration = kmc_iteration
        self.event_number = event_number
        self.method = method
        self.searcher = None

        #get the current directory and create a new directory for data to be placed here

        #create working directory and open a log file

        if method == "art":
            self.searcher = ART()
        elif method == "dimer":
            self.searcher = Dimer()
        else:
            print("kinetic monte carlo method ", method, " not recognised")
            exit(-1)
        
    # Assuming these are functions or classes defined elsewhere in the code
    # from custom_module import min_saddle_point, constructActiveRegion, makeGaussDisplacement, makeRandomDisplacement, checkVectorSize

    def find_saddle_point(self, spec, job: JobControl, fld: Field, basin1: Config, saddle_energy: Energy, continue_search, new_path, out_stream):
    
        saddle_status = Status()
        basin2_energy = Energy()
    
        #num_lanczos_vectors = job.search_parameters.numLanczosVectors
        saddle_status.set_status_failed()

        centre = np.zeros(3)

        # Make a copy of the initial basin for the saddle
        saddle = Config()
        basin1.copy_config(saddle)
        
        # Construct active region
        try:
            self.construct_active_region(saddle, spec, job, centre, out_stream)

            # Randomly displace unfrozen ions
            if job.kmc_parameters.useGauss:
                self.make_gauss_displacement(saddle, centre, job.search_parameters.initial_displacement, job.kmc_parameters.gaussWidth, out_stream)
            else:
                self.make_random_displacement(saddle, job.search_parameters.initial_displacement, out_stream)
        except Exception as e:
            out_stream.write("\n exception occured building active search region")
            out_stream.write(f"{e}\n")
            out_stream.flush()
            MPI.COMM_WORLD.Abort(4)
            exit(-1)

        if self.method == "art":
            self.check_vector_size(out_stream=out_stream)

        if job.search_parameters.go_perpendicular:
            out_stream.write("\n\n\n SEARCH using perpendicular forces ")
            self.searcher.run(saddle, self.r_vector, self.delta_pos, saddle_energy, fld, job, saddle_status, True, out_stream)
            
            #if saddle_status.get_status() == Status.FAILED:
            #    return saddle_energy, saddle_status

        # Second call to min_saddle_point
        out_stream.write("\n\n\n TRUE SEARCH for saddle point ")
        self.searcher.run(saddle, self.r_vector, self.delta_pos, saddle_energy, fld, job, saddle_status, False, out_stream)
        
        out_stream.write(f"\n status on return {saddle_status.get_status()} \n")
        if saddle_status.get_status() == Status.FAILED or saddle_status.get_status() == Status.TERMINATED:
            out_stream.write(f"\n failed to find saddle point {saddle_status.get_status()} \n")
            return basin2_energy, centre, saddle_status
        
        basin2 = Config()
        saddle.copy_config(basin2)

        if job.search_parameters.displace_basin2 == True:
            dpos = (saddle.pos - basin1.pos) * job.search_parameters.basin2_nudge
            basin2.pos = basin2.pos + dpos

        if job.search_parameters.unfreeze_basin2:
            basin2.frozen[:] = basin1.frozen[:] #use basin1 to preserve frozen types

        rel = Relax()
        basin2_energy, rel_status = rel.relax_structure(basin2, fld, basin2.natoms, job.min_parameters, out_stream)

        #make sure relaxation is ok
        if rel_status.get_status == Status.FAILED:
            out_stream.write("\n minimisation of basin2 failed \n")
            saddle_status.set_status_failed()
            return basin2_energy, centre, saddle_status

        #check whether it relaxes into new basin
        found_new_basin = self.check_for_new_basin(basin1, basin2, job.kmc_parameters.kmcBasinRadius)
        #print("sercher complete: new basin = ", found_new_basin)

        if found_new_basin == False:
            out_stream.write(f"\n failed to find new basin with distance greater than {job.kmc_parameters.kmcBasinRadius}")
            saddle_status.set_status_failed()

        # Check if the transition point was found successfully
        if saddle_status.get_status() == Status.SUCCESS:
            basin1.report_basis_difference(basin2, spec, job.kmc_parameters.kmcBasinRadius, out_stream)
        
            #write out the basin1, basin2 and saddle the correct directory
            basin2.config_write(new_path, "basin2.xyz", basin2_energy.get_total_energy())   
            saddle.config_write(new_path, "saddle.xyz", saddle_energy.get_total_energy())   
        
        return basin2_energy, centre, saddle_status
    
    def check_vector_size(self, out_stream):
        if self.vector_size == 0:
            out_stream.write("the vector size is zero and something has gone wrong!")
            exit(-1)

    def check_for_new_basin(self, basin1: Config, saddle: Config, radius: float)-> bool:
        new_basin = False

        min_radius = radius * radius
        latvector = basin1.vectors
        rcpvec = np.linalg.inv(basin1.vectors)
        
        for i in range(basin1.get_number_of_atoms()):
            rx = basin1.pos[i][0] - saddle.pos[i][0]
            ry = basin1.pos[i][1] - saddle.pos[i][1]
            rz = basin1.pos[i][2] - saddle.pos[i][2]

            xx = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
            yy = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
            zz = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]
                
            xx -= np.rint(xx)
            yy -= np.rint(yy)
            zz -= np.rint(zz)

            rx = xx * latvector[0][0] + yy * latvector[0][1] + zz * latvector[0][2]
            ry = xx * latvector[1][0] + yy * latvector[1][1] + zz * latvector[1][2]
            rz = xx * latvector[2][0] + yy * latvector[2][1] + zz * latvector[2][2]

            rsq = rx ** 2 + ry ** 2 + rz ** 2

            if rsq > min_radius:
                new_basin = True
                break

        return new_basin
        

    def construct_active_region(self, cfg : Config, spec: Species, job: JobControl, centre, out_stream):
        j = 0
        
        # Unfreeze all atoms, then freeze specific types
        #cfg.unfreeze_atoms()
        #cfg.freeze_atom_types(job.frozen_types)
        
        
        #if len(job.kmc_parameters.maskHigh) > 0:
        #    for i in range(len(job.kmc_parameters['maskHigh'])):
        #        out_stream.write(f"\n masking region between {job.kmc_parameters['maskLow'][i]} and {job.kmc_parameters['maskHigh'][i]}")
        #        cfg.freeze_region(job.kmc_parameters['maskLow'][i], job.kmc_parameters['maskHigh'][i])
        
        # Calculate vector size
        self.vector_size = (cfg.get_number_of_atoms() - cfg.num_frozen) * 3
        
        if job.kmc_parameters.regionStyle == 1:
            # All atoms except frozen ones
            if job.kmc_parameters.useGauss:
                found = False
                while not found:
                    j = int(get_random_number() * cfg.get_number_of_atoms())
                    if cfg.frozen[j] == 0:
                        found = True
                
                centre[:] = cfg.pos[j][:]
                out_stream.write(f"\n Gaussian weighting is used at the centre {centre[0]} {centre[1]} {centre[2]}")
            else:
                centre[:] = 0.0
                
            out_stream.write(f"\n no of atoms: {cfg.get_number_of_atoms()}, number frozen: {cfg.num_frozen}\n")
        
        elif job.kmc_parameters.regionStyle == 2:
            
            found = False
            radius = job.kmc_parameters.regionRadius ** 2
            while not found:
                j = int(get_random_number() * cfg.get_number_of_atoms())
                if cfg.frozen[j] == 0:
                    #print("selected ",j,cfg.frozen[j])
                    found = True

            centre[:] = cfg.pos[j][:]
            out_stream.write(f"\n atom {j} selected at random with centre {centre[0]} {centre[1]} {centre[2]}")
            latvector = cfg.vectors
            rcpvec = np.linalg.inv(cfg.vectors)
            for i in range(cfg.get_number_of_atoms()):
                if i == j or cfg.frozen[i] == 1:
                    continue

                rx = cfg.pos[i][0] - centre[0]
                ry = cfg.pos[i][1] - centre[1]
                rz = cfg.pos[i][2] - centre[2]

                xx = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
                yy = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
                zz = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]
                
                xx -= np.rint(xx)
                yy -= np.rint(yy)
                zz -= np.rint(zz)

                rx = xx * latvector[0][0] + yy * latvector[0][1] + zz * latvector[0][2]
                ry = xx * latvector[1][0] + yy * latvector[1][1] + zz * latvector[1][2]
                rz = xx * latvector[2][0] + yy * latvector[2][1] + zz * latvector[2][2]

                rsq = rx ** 2 + ry ** 2 + rz ** 2
                if rsq > radius:
                    cfg.frozen[i] = 1
                    cfg.num_frozen += 1
                    

            out_stream.write(f"\n no of atoms: {cfg.get_number_of_atoms()}, number frozen: {cfg.num_frozen}\n")
            self.vector_size = (cfg.get_number_of_atoms() - cfg.num_frozen) * 3

        elif job.kmc_parameters.regionStyle == 5:
            
            # creates list for selection of atoms that do not have the ideal coordination number i.e. a defect
            defects = self.create_coordination_list(cfg, job)
        
            #select defect atom
            if len(defects) != 0:
                j = int(get_random_number() * len(defects))
                out_stream.write(f"\n atom {j} selected from coordination defect list {len(defects)}")
            else:
                while not found:
                    j = int(get_random_number() * cfg.get_number_of_atoms())
                    if cfg.frozen[j] == 0:
                        #print("selected ",j,cfg.frozen[j])
                        found = True
                out_stream.write(f"\n atom {j} selected at random")
            #create active region and freeze the rest
            radius = job.kmc_parameters.regionRadius ** 2
            
            centre[:] = cfg.pos[j][:]
            out_stream.write(f"\n atom centre {centre[0]} {centre[1]} {centre[2]}")
            latvector = cfg.vectors
            rcpvec = np.linalg.inv(cfg.vectors)
            for i in range(cfg.get_number_of_atoms()):
                if i == j or cfg.frozen[i] == 1:
                    continue

                rx = cfg.pos[i][0] - centre[0]
                ry = cfg.pos[i][1] - centre[1]
                rz = cfg.pos[i][2] - centre[2]

                xx = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
                yy = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
                zz = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]
                
                xx -= np.rint(xx)
                yy -= np.rint(yy)
                zz -= np.rint(zz)

                rx = xx * latvector[0][0] + yy * latvector[0][1] + zz * latvector[0][2]
                ry = xx * latvector[1][0] + yy * latvector[1][1] + zz * latvector[1][2]
                rz = xx * latvector[2][0] + yy * latvector[2][1] + zz * latvector[2][2]

                rsq = rx ** 2 + ry ** 2 + rz ** 2
                if rsq > radius:
                    cfg.frozen[i] = 1
                    cfg.num_frozen += 1

            out_stream.write(f"\n no of atoms: {cfg.get_number_of_atoms()}, number frozen: {cfg.num_frozen}\n")
            self.vector_size = (cfg.get_number_of_atoms() - cfg.num_frozen) * 3

        elif job.kmc_parameters.regionStyle == 6:
            
            #reads template of atom positions and creates list for selection
            defects = self.create_defect_list(cfg, job.kmc_parameters.min_defect_distance)
            #select defect atom
            if len(defects) != 0:
                j = int(get_random_number() * len(defects))
                out_stream.write(f"\n atom {j} selected from defect list {len(defects)}")
            else:
                while not found:
                    j = int(get_random_number() * cfg.get_number_of_atoms())
                    if cfg.frozen[j] == 0:
                        #print("selected ",j,cfg.frozen[j])
                        found = True
                out_stream.write(f"\n atom {j} selected at arndom")
            #create active region and freeze the rest
            radius = job.kmc_parameters.regionRadius ** 2
            
            centre[:] = cfg.pos[j][:]
            out_stream.write(f"\n atom centre {centre[0]} {centre[1]} {centre[2]}")
            latvector = cfg.vectors
            rcpvec = np.linalg.inv(cfg.vectors)
            for i in range(cfg.get_number_of_atoms()):
                if i == j or cfg.frozen[i] == 1:
                    continue

                rx = cfg.pos[i][0] - centre[0]
                ry = cfg.pos[i][1] - centre[1]
                rz = cfg.pos[i][2] - centre[2]

                xx = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
                yy = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
                zz = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]
                
                xx -= np.rint(xx)
                yy -= np.rint(yy)
                zz -= np.rint(zz)

                rx = xx * latvector[0][0] + yy * latvector[0][1] + zz * latvector[0][2]
                ry = xx * latvector[1][0] + yy * latvector[1][1] + zz * latvector[1][2]
                rz = xx * latvector[2][0] + yy * latvector[2][1] + zz * latvector[2][2]

                rsq = rx ** 2 + ry ** 2 + rz ** 2
                if rsq > radius:
                    cfg.frozen[i] = 1
                    cfg.num_frozen += 1

            out_stream.write(f"\n no of atoms: {cfg.get_number_of_atoms()}, number frozen: {cfg.num_frozen}\n")
            self.vector_size = (cfg.get_number_of_atoms() - cfg.num_frozen) * 3

    def make_random_displacement(self, cfg: Config, initial_displacement, out_stream):
        
        if self.vector_size <= 0:
            out_stream.write("\n vectorSize is zero - something has gone wrong with atom selection")
            raise ValueError("vectorSize is zero - something has gone wrong with atom selection")

        self.delta_pos = np.zeros((cfg.get_number_of_atoms(),3))
        self.r_vector = np.zeros(self.vector_size)

        
        for i in range(cfg.get_number_of_atoms()):
            if cfg.frozen[i] == 0:
                self.delta_pos[i,0] = get_random_number() - 0.5
                self.delta_pos[i,1] = get_random_number() - 0.5
                self.delta_pos[i,2] = get_random_number() - 0.5
                

        # Prevent drift
        self.delta_pos = self.centreOfMass(self.delta_pos)

        # Normalize the displacement
        v_norm = np.linalg.norm(self.delta_pos)
        inv_norm = 1.0 / v_norm

        j = 0
        for i in range(cfg.get_number_of_atoms()):
            if cfg.frozen[i] == 0:
                self.r_vector[j] = inv_norm * self.delta_pos[i,0]
                self.r_vector[j+1] = inv_norm * self.delta_pos[i,1]
                self.r_vector[j+2] = inv_norm * self.delta_pos[i,2]
                
                j += 3
        self.delta_pos *= initial_displacement

        #self.r_vector = self.r_vector.flatten()
        
        # Displace atoms - the ASE dimer method requires this so best do it later
        #or there will be a double displacement
        if self.method == "art":
            cfg.displace_atoms(self.delta_pos)

    def make_gauss_displacement(self, cfg, centre, initial_displacement, gauss_width, out_stream):
        
        if self.vector_size <= 0:
            out_stream.write("\n vectorSize is zero - something has gone wrong with atom selection")
            raise ValueError("vectorSize is zero - something has gone wrong with atom selection")

        self.delta_pos = np.zeros((cfg.get_number_of_atoms(),3))
        self.r_vector = np.zeros(self.vector_size)

        for i in range(cfg.get_number_of_atoms()):
            if cfg.frozen[i] == 0:
                dist = cfg.get_distance(centre, i)
                weight = np.exp(-dist / (2.0 * gauss_width ** 2))
                self.delta_pos[i,0] = (get_random_number() - 0.5) * weight
                self.delta_pos[i,1] = (get_random_number() - 0.5) * weight
                self.delta_pos[i,2] = (get_random_number() - 0.5) * weight
                
        # Prevent drift
        self.delta_pos = self.centreOfMass(self.delta_pos)

        # Normalize the displacement
        v_norm = np.linalg.norm(self.delta_pos)
        inv_norm = 1.0 / v_norm

        j = 0
        for i in range(cfg.get_number_of_atoms()):
            if cfg.frozen[i] == 0:
                self.r_vector[j] = inv_norm * self.delta_pos[i,0]
                self.r_vector[j+1] = inv_norm * self.delta_pos[i,1]
                self.r_vector[j+2] = inv_norm * self.delta_pos[i,2]
                
                j += 3
        self.delta_pos *= initial_displacement

        #self.r_vector = self.r_vector.flatten()
        
        # Displace atoms - the ASE dimer method requires this so best do it later
        #or there will be a double displacement
        if self.method == "art":
            cfg.displace_atoms(self.delta_pos)

    def centreOfMass(self, vec):

        natoms = int(len(vec) / 3)

        sumX = 0.0
        sumY = 0.0
        sumZ = 0.0

        j = 0  #take account that atom positions are in separate vectors

        #determine the centre of mass
        for i in range(natoms):
            sumX += vec[j]
            sumY += vec[j+1]
            sumZ += vec[j+2]
            j += 3
    
    
        sumX /= natoms
        sumY /= natoms
        sumZ /= natoms

        #subtract the total
        j = 0
        for i in range(natoms):
            vec[j] -= sumX
            vec[j+1] -= sumY
            vec[j+2] -= sumZ

            j += 3
    
        return vec
    
    def create_coordination_list(self, cfg: Config, job: JobControl):

        defects = []
        latvector = cfg.vectors
        rcpvec = np.linalg.inv(cfg.vectors)

        
        for p in range(len(job.kmc_parameters.artTypes1)):
                        
            typ1 = job.kmc_parameters.artTypes1[p]
            typ2 = job.kmc_parameters.artTypes2[p]
            radius = job.kmc_parameters.coordDistance[p] * job.kmc_parameters.coordDistance[p]
            ideal_cn = float(job.kmc_parameters.idealCoordNumber[p])

            for i in range(cfg.natoms):

                if typ1 != cfg.symbol[i]:
                        continue
                
                if cfg.frozen[i] == 1:
                    continue
                
                num_nbrs = 0.0

                for j in range(cfg.natoms):  
                   
                    if typ2 != cfg.symbol[j]:
                        continue

                    if i == j:
                        continue

                    rx = cfg.pos[i][0] - cfg.pos[j][0]
                    ry = cfg.pos[i][1] - cfg.pos[j][1]
                    rz = cfg.pos[i][2] - cfg.pos[j][2]

                    xx = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
                    yy = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
                    zz = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]
                
                    xx -= np.rint(xx)
                    yy -= np.rint(yy)
                    zz -= np.rint(zz)

                    rx = xx * latvector[0][0] + yy * latvector[0][1] + zz * latvector[0][2]
                    ry = xx * latvector[1][0] + yy * latvector[1][1] + zz * latvector[1][2]
                    rz = xx * latvector[2][0] + yy * latvector[2][1] + zz * latvector[2][2]

                    rsq = rx * rx + ry * ry + rz * rz
                    if rsq < radius:
                        num_nbrs += 1.0

                #print("coordination number ", p, i, typ1, typ2, num_nbrs, ideal_cn)
                if np.abs(num_nbrs - ideal_cn) > 1.0e-6:
                    defects.append(i)
                
        #print("finished")
        #print("defect list ", defects)
        return defects
    
    def create_defect_list(self, cfg, min_defect_distance)-> list:

        defects = []
        radius = min_defect_distance * min_defect_distance
        latvector = cfg.vectors
        rcpvec = np.linalg.inv(cfg.vectors)
        template_pos = np.zeros((cfg.natoms,3))
        # first read in template positions
        try:
            instream = open("template", "r")
        except Exception as e:
            print(e)
            print("can not find the template file ")
            sys.exit(1)

        for i in range(cfg.natoms):
            
            line = instream.readline()

            #just in case we want to search on types as well in future
            #words = line.split()
            #name = words[0]
            
            line = instream.readline()
            words = line.split()
            x = float(words[0])
            y = float(words[1])
            z = float(words[2])
            
            template_pos[i][0] = x
            template_pos[i][1] = y
            template_pos[i][2] = z

            line = instream.readline()
            line = instream.readline()

        for i in range(cfg.natoms):
            #for j in range(cfg.natoms):     
            rx = cfg.pos[i][0] - template_pos[i][0]
            ry = cfg.pos[i][1] - template_pos[i][1]
            rz = cfg.pos[i][2] - template_pos[i][2]

            xx = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
            yy = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
            zz = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]
                
            xx -= np.rint(xx)
            yy -= np.rint(yy)
            zz -= np.rint(zz)

            rx = xx * latvector[0][0] + yy * latvector[0][1] + zz * latvector[0][2]
            ry = xx * latvector[1][0] + yy * latvector[1][1] + zz * latvector[1][2]
            rz = xx * latvector[2][0] + yy * latvector[2][1] + zz * latvector[2][2]

            rsq = rx * rx + ry * ry + rz * rz
            if rsq > radius and cfg.frozen[i] == 0:
                defects.append(i)
                
        return defects
