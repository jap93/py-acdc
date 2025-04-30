def split(s):
    return s.split()

class SearchCntrl:
    def __init__(self):
        self.print = 1
        self.max_eigenvalue = 0.0
        self.go_perpendicular = False
        self.initial_displacement = 0.02
        self.displace_basin2 = False
        self.basin2_nudge = 0.05
        self.eig_val_tol = 1.0e-3
        self.delta_x = 0.01
        self.max_par_iter = 200
        self.max_perp_iter = 100
        self.min_method = ""
        
        self.alpha_decrease = 0.99
        self.step_increase = 1.1
        self.step_decrease = 0.5
        
        self.max_time_step = 0.005
        self.init_time_step = 0.0005
        self.max_time_step_perp = 0.005
        self.init_time_step_perp = 0.001
        self.min_time_step = 0.00001
        self.N_min = 5
        self.alpha_start = 0.1
        self.norm_tol = 1.0e-3
        self.force_tol = 1.0e-4
        self.par_damp = 0.5
        self.unfreeze_basin2 = False
        self.debug = False
        self.stop_pos_eigenval = False
        self.verlet_integration = False
        self.num_lanczos_vectors = 25
        self.N_max_uphill = 50
        self.N_delay_steps = 20
        self.initial_delay = True
        
        

    def write_search_details(self, out_stream):
        num_errors = 0  # not used at the moment

        out_stream.write("\n\n *****************************************************************************************************")
        out_stream.write("\n Search Details")
        out_stream.write("\n ***************************************************************************************************** \n")

        out_stream.write(f" the number of Lanczos vectors {self.num_lanczos_vectors}\n")

        if self.num_lanczos_vectors <= 0:
            out_stream.write("\n the number of Lanczos vectors <= 0\n")
            out_stream.flush()
            exit(1)

        out_stream.write(f" the maximum value of the eigen value {self.max_eigenvalue}\n")
        if self.go_perpendicular:
            out_stream.write(" the ART forces perpendicular will be relaxed first\n")

        out_stream.write(f" the magnitude of the initial displacement {self.initial_displacement} A\n")
        out_stream.write(f" the tolenence for the eigenvalue convergence {self.eig_val_tol} energy units\n")
        out_stream.write(f" the step size for numerical derivatives {self.delta_x} A\n")
        out_stream.write(f" the maximum number of minimisation steps {self.max_par_iter}\n")
        out_stream.write(f" the maximum number of perpendicular minimisation steps {self.max_perp_iter}\n")
        out_stream.write(f" the mimimisation method (only FIRE available) {self.min_method}\n")
        out_stream.write(f" the maximum step size in FIRE {self.max_time_step}\n")
        out_stream.write(f" the initial step size in FIRE {self.init_time_step}\n")
        out_stream.write(f" the value of alpha in FIRE {self.alpha_start}\n")
        out_stream.write(f" the tolerence in the FIRE norm {self.norm_tol}\n")
        out_stream.write(f" the tolerence in the FIRE max force {self.force_tol}\n")
        # out_stream.write(f" the tolerence in the FIRE energy convergence {self.energyTol}\n")
        out_stream.write(f" the parallel force damping factor for FIRE {self.par_damp}\n")
        

        out_stream.write(f" the maximum value of the eigen value {self.max_eigenvalue}\n")
        
        out_stream.write(f" the magnitude of the initial displacement {self.initial_displacement} A\n")
        out_stream.write(f" the tolenence for the eigenvalue convergence {self.eig_val_tol} energy units\n")
        
        return num_errors

    def read_search(self, inStream, out_stream):
        readMore = True

        while readMore:
            line = inStream.readline()
            if not line:
                break

            if line[0] == '#':
                continue

            line = line.lower()
            words = split(line)

            if len(words) == 0:
                continue

            keyWord = words[0]

            if keyWord == "}":
                readMore = False
                break
            elif keyWord == "numvectors":
                self.num_lanczos_vectors = int(words[1])
            elif keyWord == "maxeigenvalue":
                self.max_eigenvalue = float(words[1])
            elif keyWord == "relaxperpendicular":
                self.go_perpendicular = True
            elif keyWord == "initialdisplacement":
                self.initial_displacement = float(words[1])
            elif keyWord == "eigentolerence":
                self.eig_val_tol = float(words[1])
            elif keyWord == "maxpariter":
                self.max_par_iter = int(words[1])
            elif keyWord == "maxperpiter":
                self.max_perp_iter = int(words[1])
            elif keyWord == "minmethod":
                self.min_method = words[1]
            elif keyWord == "maxstep":
                self.max_time_step = float(words[1])
            elif keyWord == "timestep":
                self.init_time_step = float(words[1])
            elif keyWord == "maxstepperp":
                self.max_time_step_perp = float(words[1])
            elif keyWord == "timestepperp":
                self.init_time_step_perp = float(words[1])
            elif keyWord == "alpha":
                self.alpha_start = float(words[1])
            elif keyWord == "stopposeigen":
                self.stop_pos_eigenval = True
            elif keyWord == "damp":
                self.par_damp = float(words[1])
            elif keyWord == "normtol":
                self.norm_tol = float(words[1])
            elif keyWord == "forcetol":
                self.force_tol = float(words[1])
            elif keyWord == "lanczosdisplacement":
                self.deltaX = float(words[1])
            elif keyWord == "dispbasin2":
                self.displace_basin2 = True
                self.basin2_nudge = float(words[1])
            elif keyWord == "unfreezebasin2":
                self.unfreeze_basin2 = True
            elif keyWord == "debug":
                self.debug = True
            elif keyWord == "verlet":
                self.verlet_integration = True
            else:
                out_stream.write(f"\n Dimer: keyword not found {words[0]}")
                out_stream.flush()
                exit(1)