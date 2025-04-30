import sys
import os

import numpy as np
import time as pytime

from dotenv import load_dotenv

from mpi4py import MPI

from config import Config as cfg
from job_control import JobControl as cntrl

from KineticMonteCarlo import KineticMonteCarlo
from species import Species
from field import Field

def main(args=None):


    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    #get the environment that is being used
    load_dotenv("search.env")
    #########################################################################################################################################
    # start of the kmc program
    #########################################################################################################################################
    restart_iteration = 0
    restart_time = 0.0

    outstream = None
    if rank == 0:
        try:
            outstream = open("kmc.out", 'w')

        except:
            print('cant open output file')
            sys.exit(1)

    #open the job control file and read in
    job = cntrl()
    path = os.getcwd()
    if rank == 0:
        outstream.write(f"\n the current working directory is : {path} \n")

    file_name = os.path.join(path, "control")
    try:
        instream = open(file_name, "r")
        job.readJobControl(instream, outstream)
    except Exception as e:
        print(e)
        outstream.write(f"can not find the control file or there is an error in the control file: {file_name} \n")
        outstream.flush()
        outstream.close()
        sys.exit(1)

    if rank == 0:
        job.writeSimulationDetails(outstream)

    #the positions of basin1 - note these need to be fully relaxed
    try:
        instream = open("basis.xyz", "r")

    except Exception as e:
        print("cant open configuration file")
        outstream.write(f"can not find the configuration file or there is an error reading the file: {file_name} \n")
        outstream.flush()
        outstream.close()
        sys.exit(1)

    basin1 = cfg()
    restart_iteration, restart_time, restart_energy = basin1.read_config(instream)

    instream.close()

    #if job.restart == True and np.abs(basin1.total_energy) < 1.0e-3:
    #    if rank == 0:
    #        outstream.write("\n the absolute total energy is less that 1.0e-3 and a valid total energy is required")
    #    sys.exit(1)

    #run the kinetic MC
    if job.batch == True:
        
        #start setting up and run
        #kmc = KineticMonteCarloBatch()
        #creat the batch system
        #if job.batchStyle == "archer2":
        #    bsys = ba()
        #elif job.batchStyle == "slurm":
        #    bsys = bsp()
        #else:
        #    outstream.write("\n unrecognised batch system has been specified")
        #    sys.exit(1)
        #kmc.run_kmc_batch(restart_iteration, restart_time, spec, bsys,job, basin1, outstream)
        outstream.write(f"this option has not been implemented yet \n")
        outstream.flush()
        outstream.close()
        sys.exit(1)

    else:
        fieldFileName = "potentials"
        spec = Species()
        
        fld = Field()
        kmc = KineticMonteCarlo()
        #read in species and field file
        fld.readPotential(fieldFileName, spec, outstream)
    
        try:
            with open(fieldFileName, 'r') as inStream:
                spec.read_species(inStream)
        except FileNotFoundError:
            outstream.write("\n*** could not find potentials file\n")
            sys.exit(1)
        
        if rank == 0:
            spec.print_species(outstream)

        #setup the configuration and field files
        basin1.setup_configuration(spec)
        fld.setup(basin1)

        
        #run the kmc calculation
        kmc.run_kmc(restart_energy, restart_iteration, restart_time, spec, fld, job, basin1, outstream)

    #tidy up memory and files as necessary
    fld.finalise()

    MPI.Finalize()


if __name__ == "__main__":
    main()

