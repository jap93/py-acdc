"""
this uses search method to find saddle point using batch submit method (i.e. qm codes)
"""
import os
import sys
import textwrap
import logging
import argparse
import time as pytime


#from mpi4py import MPI
#
from dotenv import load_dotenv
#


from search import Saddle
from config import Config as cfg
from job_control import JobControl as cntrl
from species import Species
from field import Field
from status import Status
from energy import Energy
from event import Event

def main(args=None):

    load_dotenv("search.env")

    log_stream = None
    try:
        log_stream = open("logfile", 'w')

    except:
        print('cant open output file')
        sys.exit(1)

    #read the job control file as this has the search parameters
    job = cntrl()
    path = os.getcwd()
    log_stream.write(f"\n the current working directory is : {path} \n")

    file_name = os.path.join(path, "control")
    try:
        instream = open(file_name, "r")
        job.readJobControl(instream, log_stream)
    except Exception as e:
        print(e)
        print("can not find the control file or there is an error in the control file: ", file_name)
        sys.exit(1)

    #the positions of basin1 - note these need to be fully relaxed
    try:
        instream = open("basis.xyz", "r")

    except FileNotFoundError:
        print("cant open configuration file")
        sys.exit(1)

    basin1 = cfg()
    restart_iteration, restart_time, restart_energy = basin1.read_config(instream)

    instream.close()

    #load the field file and species
    spec = Species()  
    fld = Field()
    
    #read in species and field file
    fieldFileName = "potentials"
    fld.readPotential(fieldFileName, spec, log_stream)
    
    try:
        with open(fieldFileName, 'r') as inStream:
            spec.read_species(inStream)
    except FileNotFoundError:
        log_stream.write("\n*** could not find potentials file\n")
        sys.exit(1)
        
    
    #setup the configuration and field files
    basin1.setup_configuration(spec)
    fld.setup(basin1)

    start_time = pytime.time()
    #do the calculation now
    searcher = Saddle(0, 0, job.kmc_parameters.kmc_method)
    basin1_energy = Energy()
    basin1_energy.vdwEnergy = restart_energy
    saddle_energy = Energy()
    continue_search = True
    new_path = "./"
    
    basin1.unfreeze_atoms()
    for typ in job.frozen_types:
        log_stream.write(f"\n freezeing atoms of type {typ}")
        basin1.freeze_atom_types(typ)
    
    basin2_energy, status = searcher.find_saddle_point(spec, job, fld, basin1, saddle_energy, continue_search, new_path, log_stream)  
    log_stream.write(f"\n\n finished saddle search {status.get_status()}\n")
    
    if status.get_status() == Status.FAILED:
        log_stream.write("\n search failed\n")
        log_stream.flush()
        exit(-1)

    act_energy = saddle_energy.get_total_energy() - basin1_energy.get_total_energy()

    if job.kmc_parameters.kmc_min_cap <= act_energy <= job.kmc_parameters.kmc_max_cap:
        log_stream.write(f"\n\n*** valid activation energy\n"
                                     f" basin 1 energy      : {basin1_energy.get_total_energy()}\n"
                                     f" basin 2 energy      : {basin2_energy.get_total_energy()}\n"
                                     f" saddle point energy : {saddle_energy.get_total_energy()}\n"
                                     f" activation energy   : {act_energy}\n")

        finish_time = pytime.time()
        elapsed_time = finish_time - start_time
        log_stream.write(f" time to calculate event : {elapsed_time:.3f} seconds\n")

        search_number = 0        
        ev = Event(basin1_energy.get_total_energy(), basin2_energy.get_total_energy(), saddle_energy.get_total_energy(), act_energy, search_number)
        ev.calculate_rate(job.kmc_parameters.kmcTemperature, job.kmc_parameters.kmcPreFactor)
        
        #save the events
        fp = open("events", "w")
        ev.write_event(fp)
        fp.flush()
        fp.close()
                
             
    else:
        log_stream.write(f"\n invalid activation energy {act_energy}\n")
        status.set_status_failed()
            
    #successful or otherwise write out the status for use later
    status.write("")
    log_stream.flush()
    log_stream.close()

    del searcher

if __name__ == "__main__":
    main()
