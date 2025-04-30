import sys
from mpi4py import MPI

class RelaxCntrl:
    def __init__(self):
        self.max_iter = 500
        self.print = 1
        self.norm_tol = 1.0e-5
        self.force_tol = 1.0e-4
        self.N_min = 5
        self.init_time_step = 0.0005
        self.alpha_decrease = 0.99
        self.step_increase = 1.1
        self.step_decrease = 0.5
        self.alpha_start = 0.1
        self.max_time_step = 0.001
        self.debug = False
        self.method = "fire"
        self.verlet_integration = False
        self.N_maxUphill = 50
        self.N_delaySteps = 20
        self.min_time_step = 0.00001
        self.initialDelay = True


    def split(self, s):
        return s.split()

    def readRelax(self, inStream, out_stream):
        readMore = True

        while readMore:
            line = inStream.readline()
            if not line:
                break

            if line[0] == '#':
                continue

            line = line.lower()
            words = self.split(line)
            if not words:
                continue

            keyWord = words[0]

            if keyWord == "}":
                readMore = False
                break
            elif keyWord == "forcetol":
                self.force_tol = float(words[1])
            elif keyWord == "debug":
                self.debug = True
            elif keyWord == "method":
                self.method = words[1]
            elif keyWord == "relaxsteps":
                self.max_iter = int(words[1])
            elif keyWord == "maxstep":
                self.max_time_step = float(words[1])
            elif keyWord == "timestep":
                self.init_time_step = float(words[1])
            elif keyWord == "alpha":
                self.alpha_start = float(words[1])
            elif keyWord == "verlet":
                self.verlet_integration = True
            else:
                out_stream.write("\n RELAX: keyword not found " + words[0] + "\n")
                out_stream.flush()
                MPI.COMM_WORLD.Abort()
                sys.exit(1)

    def __eq__(self, src):
        if self is src:
            return self

        self.max_iter = src.max_iter
        self.print = src.print
        self.norm_tol = src.norm_tol
        self.force_tol = src.force_tol
        self.N_min = src.N_min
        self.init_time_step = src.init_time_step
        self.alpha_decrease = src.alpha_decrease
        self.step_increase = src.step_increase
        self.step_decrease = src.step_decrease
        self.alpha_start = src.alpha_start
        self.max_time_step = src.max_time_step
        self.debug = src.debug
        self.method = src.method
        self.verlet_integration = src.verlet_integration
        self.N_maxUphill = src.N_maxUphill
        self.N_delaySteps = src.N_delaySteps
        self.min_time_step = src.min_time_step
        self.initialDelay = src.initialDelay

        return self