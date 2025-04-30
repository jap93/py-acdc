import os

import numpy as np
from scipy import linalg

from species import Species

class Config (object):
    
    def __init__(self):
        
        self.title = 'untitled config'
        self.levcfg = 0
        self.pbc = True
        self.natoms = 0
        self.num_frozen = 0
        self.symbol = []
        self.vectors = np.zeros((3,3))
        self.charge = None
        self.mass = None
        self.pos = None
        self.frozen = None
        self.label = None
        self.total_energy = 0.0

    def reset_config(self):

        self.title = 'untitled config'
        self.levcfg = 0
        self.pbc = 3
        self.natoms = 0
        self.vectors = np.zeros((3,3))
        self.total_energy = 0.0
        del self.symbol [:]
        self.pos = None
        self.label = None
        self.frozen = None
        self.charge = None
        self.mass = None

    def setup_configuration(self, spec: Species):

        self.charge = np.zeros(self.natoms, dtype=np.float64)
        self.mass = np.zeros(self.natoms, dtype=np.float64)
        self.label = np.zeros(self.natoms, dtype=np.int32)

        for i in range(self.natoms):

            for k in range(spec.get_num_species()):
                element = spec.get_species(k)
                if self.symbol[i] == element.name:
                    self.label[i] = k 
                    self.charge[i] = element.charge
                    self.mass[i] = element.mass


    def freeze_atom_types(self, typ):

        self.num_frozen = 0
        for i in range(self.natoms):
            #if self.symbol[i] == typ:
            if typ in self.symbol[i]:
                self.frozen[i] = 1
                self.num_frozen += 1

    def unfreeze_atoms(self):
        self.num_frozen = 0
        self.frozen.fill(0)

    #returns the number of atoms in a list to controling program
    def get_number_of_atoms(self):
        return self.natoms
        
    def get_atom(self, i):
        return self.pos[i]

    def reset_number_of_atoms(self):
        self.natoms = len(self.pos)
    
        #calculates volume assuming an orthorhombic cell
    def get_volume(self):
        axb = np.cross(self.vectors[0][:], self.vectors[1][:])
        self.volume = np.dot(axb, self.vectors[2][:])

        return self.volume
        
                  
    #create a copy of the config
    def copy_config(self, c):
        
        c.title = self.title
        c.levcfg = self.levcfg
        c.pbc = self.pbc
        c.natoms = self.natoms
        c.num_frozen = self.num_frozen
        c.symbol = self.symbol

        c.pos = np.zeros((self.natoms,3), dtype=np.float64)
        c.charge = np.zeros(self.natoms, dtype=np.float64)
        c.mass = np.zeros(self.natoms, dtype=np.float64)
        c.label = np.zeros(self.natoms, dtype=np.int32)
        c.frozen = np.zeros(self.natoms, dtype=np.int32)
        
        np.copyto(c.vectors, self.vectors)

        c.total_energy = self.total_energy
    
        np.copyto(c.pos, self.pos)
        np.copyto(c.charge, self.charge)
        np.copyto(c.mass, self.mass)

        np.copyto(c.label, self.label)

        np.copyto(c.frozen, self.frozen)
            
    
    #dtermines how far atoms have moved - it assumes both configurations
    #have the same order of atoms and both cells are of the same size
    def report_basis_difference(self, new_cfg, spec, radius, out_stream):
        min_radius = radius * radius
        latvector = self.vectors
        rcpvec = np.linalg.inv(self.vectors)
        
        for i in range(self.get_number_of_atoms()):
            rx = self.pos[i][0] - new_cfg.pos[i][0]
            ry = self.pos[i][1] - new_cfg.pos[i][1]
            rz = self.pos[i][2] - new_cfg.pos[i][2]

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

            if rsq > min_radius:
                out_stream.write(f"\n\n displacement of atom {i} {self.pos[i][0]} {self.pos[i][1]} {self.pos[i][2]}")
                out_stream.write(f"\n displacement of atom {i} {new_cfg.pos[i][0]} {new_cfg.pos[i][1]} {new_cfg.pos[i][2]}")
                out_stream.write(f"\n displacement of atom {i} {rx} {ry} {rz} {np.sqrt(rsq)}")

    def find_active_region(self, centre, radius):
        min_radius = radius * radius
        latvector = self.vectors
        rcpvec = np.linalg.inv(self.vectors)

        num_list = []
        px = []
        py = []
        pz = []
        
        for i in range(self.get_number_of_atoms()):
            rx = self.pos[i][0] - centre[0]
            ry = self.pos[i][1] - centre[1]
            rz = self.pos[i][2] - centre[2]

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

            if rsq <= min_radius:
                num_list.append(i)
                px.append(self.pos[i][0])
                py.append(self.pos[i][1])
                pz.append(self.pos[i][2])
                
        return num_list, px, py, pz

    def get_distance(self, centre, i):
        
        latvector = self.vectors
        rcpvec = np.linalg.inv(self.vectors)

        rx = self.pos[i][0] - centre[0]
        ry = self.pos[i][1] - centre[1]
        rz = self.pos[i][2] - centre[2]

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
                
        return np.sqrt(rsq)
    
    def difference(self, cfg):
        
        #check the configs have the same number of atoms
        if len(self.pos) != len(cfg.pos):
            print (' configurations have different no of atoms')
            return
        
        difference = np.zeros(len(self.pos), 3)
        
        for i in range(len(self.pos)):
            
            p1 = self.pos[i][:]
            p2 = cfg.pos[i][:]
            #convert atoms to fractional coordinates
            
            #calculate distance
            print("needs to be done")
            exit(1)
            #difference[i][:] = get_distance(p1, p2)

        return difference
        
    def displace_atoms(self, delta):
        for i in range(self.natoms):
            if self.frozen[i] == 0: #i.e. it is free to move
                self.pos[i,0] += delta[i,0]
                self.pos[i,1] += delta[i,1]
                self.pos[i,2] += delta[i,2]
        
    def displace_atoms_mapped(self, delta):
        j = 0
        for i in range(self.natoms):
            if self.frozen[i] == 0: #i.e. it is free to move
                self.pos[i,0] += delta[j]
                self.pos[i,1] += delta[j+1]
                self.pos[i,2] += delta[j+2]
                j +=3

    #dtermines how far atoms have moved - it assumes both configurations
    #have the same order of atoms and both cells are of the same size
    def distance_moved(self, cfg, distance):
        
        radius = distance * distance
        difference = np.zeros(3)
        
        #check the configs have the same number of atoms
        if len(self.pos) != len(cfg.data):
            print (' configurations have different no of atoms')
            
        for i in range(len(self.pos)):
            
            atm1 = self.pos[i]
            atm2 = cfg.data[i]
            #convert atoms to fractional coordinates
            
            #calculate distance
            difference = atm1.get_distance(atm2.get_position())
            
            dist = difference[0] * difference[0] + difference[1] * difference[1] + difference[2] * difference[2]
            
            if dist > radius:
                
                dist = np.sqrt(dist)
                print ('displacement ', atm1.symbol, i, atm1.posn, atm2.posn, dist)
                #print ('displacement ', atm1.symbol, i, difference, dist)

    #calculates running total for rdf's
    def calculate_interactions(self, nbins, r_min, r_max, tmp_gofr):
        
        deltar = (r_max - r_min) / float(nbins)
        natoms = self.natoms
        
        box_x = self.vectors[0][0]
        box_y = self.vectors[1][1]
        box_z = self.vectors[2][2]
        
        p = np.zeros((natoms,3))
        for i in range(natoms):
            atmi = self.pos[i]
            p[i][0] = atmi.posn[0]
            p[i][1] = atmi.posn[1]
            p[i][2] = atmi.posn[2]
        
        for i in range(natoms - 1):
            
            ix = p[i][0]
            iy = p[i][1]
            iz = p[i][2]
            
            for j in range(i + 1, natoms):
                
                #atmj = self.pos[j]
                
                xx = p[j][0] - ix
                yy = p[j][1] - iy
                zz = p[j][2] - iz
                
                #allow for periodic boundary conditions
                rx = xx - box_x * np.rint(xx / box_x)
                ry = yy - box_y * np.rint(yy / box_y)
                rz = zz - box_z * np.rint(zz / box_z)
                
                rsq = rx * rx + ry * ry + rz * rz
                r = np.sqrt(rsq)
                
                if r <= r_max and r >= r_min:
                    
                    #try:
                    ibin = int(r / deltar) - 1
                        # determine the correct bin
                    tmp_gofr[ibin] += 2        
                        
                    #except:
                    #    print "ibin out of range"
                    #    print r, deltar, ibin
                    #    sys.exit(1)
            
    

    def find_numatom_of_type(self, typ):
        
        n = 0
        
        for i in range(self.natoms):
            
            tmp = self.pos[i]
            
            if (tmp.symbol == typ):
                n += 1
        
        return n

    
    def centre_configuration(self):
        """
        the config is setup to run from -0.5 to 0.5
        """
        
        self.cartesian_to_fractional()
        
        for i in range(self.natoms):
            
            self.pos[i][0] -= np.rint(self.pos[i][0])
            self.pos[i][1] -= np.rint(self.pos[i][1])
            self.pos[i][2] -= np.rint(self.pos[i][2])
            
        self.fractional_to_cartesian()

    def start_at_origin(self):

        self.cartesian_to_fractional()

        for i in range(self.natoms):

            a = self.pos[i]

            if self.pos[i][0] < 0.0: 
                self.pos[i][0] += 1.0

            if self.pos[i][1] < 0.0: 
                self.pos[i][1] += 1.0

            if self.pos[i][2] < 0.0: 
                self.pos[i][2] += 1.0

        self.fractional_to_cartesian()
        
    def cartesian_to_fractional(self):
        
        rcpvec = linalg.inv(self.vectors)

        for i in range(self.natoms):
            
            rx = self.pos[i][0]
            ry = self.pos[i][1]
            rz = self.pos[i][2]

            self.pos[i][0] = rx * rcpvec[0][0] + ry * rcpvec[0][1] + rz * rcpvec[0][2]
            self.pos[i][1] = rx * rcpvec[1][0] + ry * rcpvec[1][1] + rz * rcpvec[1][2]
            self.pos[i][2] = rx * rcpvec[2][0] + ry * rcpvec[2][1] + rz * rcpvec[2][2]

    def fractional_to_cartesian(self):
        
        latvector = self.vectors
        
        for i in range(self.natoms):
            
            rx = self.pos[i][0]
            ry = self.pos[i][1]
            rz = self.pos[i][2]

            self.pos[i][0] = rx * latvector[0][0] + ry * latvector[0][1] + rz * latvector[0][2]
            self.pos[i][1] = rx * latvector[1][0] + ry * latvector[1][1] + rz * latvector[1][2]
            self.pos[i][2] = rx * latvector[2][0] + ry * latvector[2][1] + rz * latvector[2][2]      

    def read_config(self, instream):
        """
        reads in the simplified xyz file
        """
        restart_iteration = 0
        restart_time = 0.0
        restart_energy = 0.0

        line = instream.readline()
        words = line.split()
        self.natoms = int(words[0])

        self.pos = np.zeros((self.natoms, 3), dtype=np.float64)
        self.label = np.zeros(self.natoms, dtype=np.int32)
        self.frozen = np.zeros(self.natoms, dtype=np.int32)

        #read in the information line containing vectors etc
        info = instream.readline()
        info = info.lower()
        words = info.split()

        idx = 0
        for w in words:
            if "lattice=" in w:
                break
            idx +=1

        tmp = words[idx].replace('lattice=','')
        tmp = tmp.replace('"','')
        self.vectors[0,0] = float(tmp)
        self.vectors[0][1] = float(words[idx+1])
        self.vectors[0][2] = float(words[idx+2])
        self.vectors[1][0] = float(words[idx+3])
        self.vectors[1][1] = float(words[idx+4])
        self.vectors[1][2] = float(words[idx+5])
        self.vectors[2][0] = float(words[idx+6])
        self.vectors[2][1] = float(words[idx+7])
        tmp = words[idx+8].replace('"','')
        self.vectors[2][2] = float(tmp)

        for w in words:
            if "iteration" in w:
                restart_iteration = int(w.replace('iteration=',''))
                break

        for w in words:
            if "time" in w:
                restart_time = float(w.replace('time=',''))
                break

        for w in words:
            if "energy" in w:
                restart_energy = float(w.replace('energy=',''))
                break
    
        #next read in the atom positions
        for i in range(self.natoms):
            
            line = instream.readline()
            
            words = line.split()
            self.symbol.append(words[0])
            
            x = float(words[1])
            y = float(words[2])
            z = float(words[3])
            
            self.pos[i][0] = x
            self.pos[i][1] = y
            self.pos[i][2] = z

       
        return restart_iteration, restart_time, restart_energy
        
    def config_write(self, new_path, out_name, total_energy=None, restart_iteration = None, restart_time = None):
        """
        writes out a simplified form of the extended xyz format used by ASE
        """
        file_name = os.path.join(new_path, out_name)
        outstream = open(file_name, "w")

        outstream.write(f"{self.natoms} \n")
        
        #prepare the info line
        info = 'Lattice="'
        tmp = self.vectors.flatten()
        for i in range(9):
            info = info + str(tmp[i]) + ' '
        info = info + '"  '

        if total_energy != None:
            info = info + "energy=" + str(total_energy) + "  "

        if restart_iteration != None:
            info = info + "iteration=" + str(restart_iteration) + " "

        if restart_time != None:
            info = info + "time=" + str(restart_time) + " "

        info = info + 'pbc="T T T"  '  #not used within this prog but may be externally so included

        outstream.write(f"{info} \n")
    
        for i in range(self.natoms):
            outstream.write(f"{self.symbol[i]}    {self.pos[i][0]}  {self.pos[i][1]}  {self.pos[i][2]}\n")


        outstream.flush()
        outstream.close()
        
    def cell_properties(self, aaa_matrix):
        # aaa is a flat array of size 9 (3x3 matrix in row-major order)
        bbb = np.zeros(10)
    
        
        # Calculate lengths of cell vectors
        bbb[0] = np.linalg.norm(aaa_matrix[0, :])  # Length of first row vector
        bbb[1] = np.linalg.norm(aaa_matrix[1, :])  # Length of second row vector
        bbb[2] = np.linalg.norm(aaa_matrix[2, :])  # Length of third row vector

        # Calculate cosines of cell angles
        bbb[3] = np.dot(aaa_matrix[0, :], aaa_matrix[1, :]) / (bbb[0] * bbb[1])
        bbb[4] = np.dot(aaa_matrix[0, :], aaa_matrix[2, :]) / (bbb[0] * bbb[2])
        bbb[5] = np.dot(aaa_matrix[1, :], aaa_matrix[2, :]) / (bbb[1] * bbb[2])

        # Calculate vector products of cell vectors (cross products)
        axb = np.cross(aaa_matrix[0, :], aaa_matrix[1, :])
        bxc = np.cross(aaa_matrix[1, :], aaa_matrix[2, :])
        cxa = np.cross(aaa_matrix[2, :], aaa_matrix[0, :])

        # Calculate volume of the cell
        bbb[9] = abs(np.dot(aaa_matrix[0, :], bxc))

        # Calculate cell perpendicular widths
        bbb[6] = bbb[9] / np.linalg.norm(bxc)
        bbb[7] = bbb[9] / np.linalg.norm(cxa)
        bbb[8] = bbb[9] / np.linalg.norm(axb)
    
        return bbb
