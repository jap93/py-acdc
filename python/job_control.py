import sys

from mpi4py import MPI

from KMCControl import KMCCntrl
from relax_control import RelaxCntrl
from search_control import SearchCntrl

class JobControl(object):

    def __init__(self):
        #the number of cores per calculation */
        #self.number_cores = 0,

        #the number of workwork groups
        self.num_workgroups = 1

        #seed for random number generator */
        self.seed = 0

        #time delay between queries of job status
        #self.time_delay = 60 #seconds

        #flag to indicate a restart of the program */
        self.restart = False

        #format for the starting geometry
        #self.inFormat = "dlpoly"
        
        # vector of names of atom types to be frozen in simulations */
        self.frozen_types = []

        # how the jobs are executed
        self.batchStyle = "archer2"

        #the number of excess kmc jobs to submit
        self.excessImages = 0

        #flag to indicate restrt of the Kinetic Monte Carlo part of the code
        self.restart = False

        #flag to request that the code attempts to recover from checkpoint files
        self.recover = False

        #flag to use batch system
        self.batch = False

        #flag indicates that the code will stop when it reaches a checkpont
        self.exit_at_checkpoint = False

        self.kmc_parameters = None
        self.min_parameters = None
        self.search_parameters = None

    def writeSimulationDetails(self, out_stream):

        num_errors = 0

        out_stream.write("\n\n *****************************************************************************************************")
        out_stream.write("\n Simulation Details")
        out_stream.write("\n *****************************************************************************************************")


        # key words for both internal and external energy evaluation */

        if len(self.frozen_types):
    
            out_stream.write("\n atom types that will be frozen by ACDC ")
            for i in range(len(self.frozen_types)):
                out_stream.write("\n type    {}".format(self.frozen_types[i]))
        
        out_stream.write("\n")  

    
        out_stream.write("\n the number of MPI workgroups {:d} \n".format(self.num_workgroups))

        #out_stream.write("\n time delay between job status query {:d} seconds \n".format(self.time_delay))

        if self.restart == True:
            out_stream.write("\n the simulation will be restarted ")

        #if self.preRelax == True:
        #    out_stream.write("\n the cell positions will be relaxed prior to KMC start \n")
        #else:
        #    out_stream.write("\n only the energy will be calculated prior to KMC start \n")

        # key words for external energy evaluation */
        
        #out_stream.write("\n the excess number of KMC events that will be launched {:8d} \n".format(self.excessImages))

        if self.restart == True:
            out_stream.write("\n the simulation will be restarted")
    
        #if self.kmcMethod != "art" and self.kmcMethod != "dimer":
        #    num_errors += 1
        #    out_stream.write("\n only ART+ Dimer is supported at the moment! \n")

        #if self.recover == True:
        #    out_stream.write("\n the calculation will commence from checkpoint file")

        #if self.exit_at_checkpoint == True:
        #    out_stream.write("\n the calculation will stop after checkpoint")
        
        num_errors += self.kmc_parameters.writeKMCDetails(out_stream)
        
        num_errors += self.search_parameters.write_search_details(out_stream)
        
        #if there have been any input errors / conflicts abort */
        if num_errors > 0:
            print("errors have been found in the input")
            out_stream.write("\n\n there are {:d} errors in job control input that must be corrected before continuing".format(num_errors))
            out_stream.flush()
            sys.exit(-1)
    
    

    def readJobControl(self, inStream, out_stream):
        readMore = True

        self.search_parameters = SearchCntrl()
        self.min_parameters = RelaxCntrl()
        self.kmc_parameters = KMCCntrl()

        #start reading the job control parameters
        while readMore == True:
    
            #keep on reading until end of file
            line = inStream.readline()
        
            if line[0] == '#':
                continue

            words = line.split()

            if len(words) == 0:
                continue

            keyWord = words[0].lower()
        

            if keyWord == "close":
                readMore = False
           
            #elif keyWord == "prerelax":
            #    self.preRelax = True

            elif keyWord ==  "externalfile":
                self.externalFile = words[1]
        
            elif keyWord ==  "rootname":
                self.pubRootName = words[1]
        
            elif keyWord ==  "restart":
                self.restart = True

            elif keyWord ==  "timedelay":
                self.time_delay = int(words[1])

            elif keyWord == "batch":
                self.batchStyle = words[1]
        
            elif keyWord ==  "freeze":
            
                num = int(words[1])
                for i in range(num):
                    line = inStream.readline()
                    words = line.split()
                    self.frozen_types.append(words[0])
            
            elif keyWord == "numgroups":
                self.numGrps = int(words[1])
        
            elif keyWord == "seed":
                self.seed = int(words[1])
        
            elif keyWord == "restart":
                self.restart = True

            elif keyWord == "kmc{":
                self.kmc_parameters.readKineticMonteCarlo(inStream, out_stream)

            elif keyWord == "relax{":
                self.min_parameters.readRelax(inStream, out_stream)

            elif keyWord == "search{":
                self.search_parameters.read_search(inStream, out_stream)
            

            #else:
            #    out_stream.write("\n keyword not found " + words[0])
            #    sys.exit(-1)



 
