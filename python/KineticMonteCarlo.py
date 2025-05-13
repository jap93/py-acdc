from mpi4py import MPI
import numpy as np
import time
import random
import sys
import os
import traceback

from event import Event
from config import Config
from energy import Energy
from relax import Relax
from search import Saddle
from event import Event
from status import Status
from job_control import JobControl
from constants import SUCCESS, TERMINATED

EXIT_FAILURE = 1

class KineticMonteCarlo:
    def __init__(self):
        pass
    
    def run_kmc(self, restart_energy, restart_iteration, restart_time, spec, fld, job, basin1, out_stream):
        ievent = 0
        iteration = 0
        comm = MPI.COMM_WORLD
        world_rank = comm.Get_rank()
        num_procs = comm.Get_size()

        
        time_increase = 0.0
        kmc_time = 0.0
        ran_number = 0.0
        total_rate = 0.0
        temp_rate = 0.0

        converged = False
        prime = (world_rank == 0)

        event_list = []
        excluded_list = []
        #excluded_list = [0] * basin1.natoms if job.kmc_parameters.super_basin else []
        
        if job.restart:
            iteration = restart_iteration
            kmc_time = restart_time
            if world_rank == 0:
                out_stream.write("\n\n *****************************************************************************************************\n")
                out_stream.write(f"\n Restart undertaken at iteration : {iteration} for a kmc time of : {kmc_time}\n")
                out_stream.write(" *****************************************************************************************************\n")
                out_stream.flush()
        
        while not converged:
            iteration += 1

            if job.restart and job.kmc_parameters.recycle_events and world_rank == 0:
                event_list = self.recover_event_list(iteration, job.kmc_parameters.kmcTemperature, job.kmc_parameters.kmcPreFactor, out_stream)
                out_stream.write(f"\n then event list recovered from iteration {iteration} is {len(event_list)} \n")

            if world_rank == 0:
                out_stream.write("\n\n *****************************************************************************************************\n")
                out_stream.write(f" Kinetic Monte Carlo iteration : {iteration}\n")
                out_stream.write(" *****************************************************************************************************\n")
                out_stream.flush()

            basin1.unfreeze_atoms()
            for typ in job.frozen_types:
                basin1.freeze_atom_types(typ)

            rel = Relax()
            relaxed_energy, rel_status = rel.relax_structure(basin1, fld, basin1.natoms, job.min_parameters, out_stream)
            
            basin1_energy = relaxed_energy.get_total_energy()

            if rel_status.get_status() == Status.FAILED:
                out_stream.write("\n failed to relax atoms at the start of the KMC iteration\n")
                out_stream.flush()
                MPI.COMM_WORLD.Abort(EXIT_FAILURE)
                sys.exit(EXIT_FAILURE)

            if world_rank == 0:
                out_stream.write(f"\n total energy of initial relaxation {basin1_energy:.10e}\n")

            start_time = MPI.Wtime()

            if num_procs == 1:
                self.calculate_rate_constant_single_core(spec, fld, job, basin1, relaxed_energy, event_list, iteration, kmc_time, out_stream)
            else:
                self.calculate_rate_constant_experimental(spec, fld, job, basin1, relaxed_energy, event_list, excluded_list, iteration, kmc_time, out_stream)
                MPI.COMM_WORLD.Barrier()

            if prime:
                total_rate = sum(event.get_rate() for event in event_list)
                ran_number = np.random.random() * total_rate

                temp_rate = 0.0
                search_number = 0
                for i, event in enumerate(event_list):
                    temp_rate += event.get_rate()
                    if temp_rate > ran_number:
                        ievent = i
                        search_number = event.get_search_number()
                        activation_energy = event.get_activation_energy()
                        break
                
                out_stream.write(f"\n\n kmc iteration : {iteration} chosen event : {ievent}\n")

                time_increase = self.calculate_time(total_rate)
                kmc_time += time_increase

                for i, e in enumerate(event_list):
                    out_stream.write(f" number: {i} rate : {e.get_rate():.5e} e_act : {e.get_activation_energy():.5e} search: {e.get_search_number()} \n")

                out_stream.write(f" total rate {total_rate:.5e}\n")
                out_stream.write(f"\n kinetic monte time step : {time_increase:.5e} now : {kmc_time:.5e}\n")

                finish_time = MPI.Wtime()
                diff = finish_time - start_time
                out_stream.write(f"\n elapsed time for this iteration : {diff:.5f} seconds\n")
                start_time = finish_time

                basin1.reset_config()
                tmp = "iteration" + str(iteration)
                file_name = os.path.join(tmp,"selection")
                with open(file_name, "w") as selectionstream:
                    selectionstream.write(f"event {ievent} search {search_number} E_a {activation_energy}")

                tmp1 = "search" + str(search_number)
                file_name = os.path.join(tmp,tmp1,"basin2.xyz")
                out_stream.write(f"\n resetting basin1 from file  {file_name} \n")

                try:
                    with open(file_name, "r") as instream:
                        dummy, dumm1, restart_energy = basin1.read_config(instream)
                except Exception as e:
                    out_stream.write(f"{e}\n")
                    out_stream.flush()
                    traceback.print_exc()
                    MPI.COMM_WORLD.Abort()
                    exit()
                
                basin1.setup_configuration(spec)

                #create the restart file on each iteration
                basin1.config_write("", "restart.xyz", restart_energy, restart_iteration=iteration, restart_time=kmc_time)

                if job.kmc_parameters.recycle_events:
                    if job.kmc_parameters.recycle_style < 2:
                        event_list = self.recycle_saddle_search(event_list, basin1, restart_energy, ievent, job.kmc_parameters.regionRadius, job.kmc_parameters.recycle_distance, iteration, 
                                                                job.kmc_parameters.kmc_min_cap, job.kmc_parameters.recycle_prob, True, out_stream)
                    else:
                        event_list = self.recycle_saddle_limit(event_list, basin1, restart_energy, ievent, job.kmc_parameters.regionRadius, job.kmc_parameters.recycle_distance, iteration, 
                                                                job.kmc_parameters.kmc_min_cap, job.kmc_parameters.recycle_limit, True, out_stream)
                else:
                    event_list.clear()
        

            #if more than one mpi thread is being used then allcores must have a copy of the basin
            # basin1_energy is not required as a relaxation is done        
            if num_procs > 1:
                self.broadcast_positions(basin1)

            if iteration >= job.kmc_parameters.kmc_steps:
                converged = True
            
            
                
                #if job.kmc_parameters.super_basin:
                #    with open("excludedlist", 'w') as ex_list_stream:
                #        self.write_excluded_list(excluded_list, ex_list_stream)
                #        ex_list_stream.flush()

            ierr = MPI.COMM_WORLD.Barrier()

        if world_rank == 0:
            basin1.config_write("", "restart.xyz", restart_energy, restart_iteration=iteration, restart_time=kmc_time)
            

        ierr = MPI.COMM_WORLD.Barrier()

    def calculate_rate_constant(self, spec, fld, job: JobControl, basin1, basin1_energy: Energy, event_list, excluded_list, iteration, kmc_time, out_stream):
        comm = MPI.COMM_WORLD
        world_rank = comm.Get_rank()
        num_procs = comm.Get_size()

        num_events = len(event_list)
        search_number = -1

        if job.kmc_parameters.recycle_events == True:
            search_number = len(event_list) - 1

        act_energy = 0.0

        basin2_energy = Energy()
        saddle_energy = Energy()
        
        continue_search = False

        searcher = None
        status = Status()
        #ev = Event()
        attempted_search = 0
        successfull_search = 0
        num_busy = 0

        it_path = "iteration" + str(iteration)
        if world_rank == 0:
            if job.kmc_parameters.recycle_events == False:
                try:
                    os.mkdir(it_path)
                except Exception as e:
                    out_stream.write(f"{e}\n")
                    out_stream.flush()
                    MPI.COMM_WORLD.Abort()
                    exit()
            elif iteration == 1:
                try:
                    os.mkdir(it_path)
                except Exception as e:
                    out_stream.write(f"{e}\n")
                    out_stream.flush()
                    MPI.COMM_WORLD.Abort()
                    exit()

            basin1.config_write(it_path, "basin1.xyz", basin1_energy.get_total_energy(), iteration, kmc_time) # basin1 only needs to be written out once
        
        comm.Barrier()  # dont let processors get ahead and create a search directory before iteration

        if world_rank == 0:  ################################ parent process to initiate and keep track of events
            num_search = num_procs
            left_over = job.kmc_parameters.kmc_images - num_events
            if left_over < num_procs:
                num_search = left_over
            #send a command to all groups to start search
            for i in range(1, num_search):
                #out_stream.write(f"start  {i} {search_number} \n")
                #out_stream.flush()
                search_number += 1
                attempted_search +=1
                num_busy += 1
                self.send_command_to_child(search_number,i)


            while num_events < job.kmc_parameters.kmc_images:

                #wait to receive an event
                command, source, tag = self.check_for_message()
                
                num_busy -= 1
                #print("received message frm", source, "tag", tag, "command", command)

                #check to see whether it is a valid event
                if np.rint(command[0]) >= 0:
                    id = int(command[1])
                    #out_stream.write(f"command received by parent {source} {tag} {num_events} {id} \n")
                    #out_stream.flush()
                    basin2_energy.vdwEnergy = np.float64(command[3])
                    saddle_energy.vdwEnergy = np.float64(command[4])
                    act_energy = np.float64(command[5])
                    centre = np.zeros(3)
                    centre[0] = command[6]
                    centre[1] = command[7]
                    centre[2] = command[8]
                    ev = Event(basin1_energy.get_total_energy(), basin2_energy.get_total_energy(), saddle_energy.get_total_energy(), act_energy, centre, id)
                    ev.calculate_rate(job.kmc_parameters.kmcTemperature, job.kmc_parameters.kmcPreFactor)
                    event_list.append(ev)
                    num_events += 1
                    successfull_search += 1
                    #save the events
                    eventname = os.path.join(it_path, "events")
                    fp = open(eventname, "a")
                    #for e in event_list:   #make a copy of all valid events for restart
                    ev.write_event(fp)
                    fp.flush()
                    fp.close()

                #check to see if enough searches have been launched
                if job.kmc_parameters.kmc_images - (num_events + num_busy) >= 0:
                    #send the search number to the child to initiate search again
                    search_number += 1
                    attempted_search +=1
                    num_busy += 1
                    self.send_command_to_child(search_number, source)

            #have all the number of events so send signal to children for return
            self.send_signal_to_children()
                
            out_stream.write(f"\n\n attempted searches : {attempted_search} successfull searches : {successfull_search}")
        else:   ######################################### child processors
            command = 0
            while command >= 0:
                command, source, tag = self.check_for_message()
                
                if command < 0:
                    break

                if command >= 0:
                    search_number = command
                    flnm = "search" + str(search_number)
                    #create a new directory for this iteration and attemp for the search
                    new_path = os.path.join(it_path, flnm)
                    try:
                        os.mkdir(new_path)
                    except Exception as e:
                        sys.stderr.write(f"{e}\n")
                        sys.stderr.flush()
                        MPI.COMM_WORLD.Abort()
                        exit()
                    #create a stream for the log file to go into this directory
                    log_path = os.path.join(new_path, "logfile")
                    try:
                        log_stream = open(log_path, "w")
                    except Exception as e:
                        sys.stderr.write(f"{e}\n")
                        sys.stderr.flush()
                        MPI.COMM_WORLD.Abort()
                        exit()


                    start_time = time.time()
                    attempted_search += 1
                    searcher = Saddle(iteration, num_events, job.kmc_parameters.kmc_method)
                    #try:
                    basin2_energy, centre, status = searcher.find_saddle_point(spec, job, fld, basin1, saddle_energy, continue_search, new_path, log_stream)  
                    log_stream.write(f"\n\n finished saddle search {status.get_status()}\n")
                    #except Exception as e:
                    #    log_stream.write("\n kinetic_monte_carlo : ")
                    #    log_stream.write(f"{e}\n")
                    #    status.set_status_failed()

                    if status.get_status() == Status.FAILED:
                        log_stream.write("\n search failed\n")
                
                    if status.get_status() == Status.SUCCESS:
                        act_energy = saddle_energy.get_total_energy() - basin1_energy.get_total_energy()

                        if act_energy >= job.kmc_parameters.kmc_min_cap and act_energy <= job.kmc_parameters.kmc_max_cap:
                            log_stream.write(f"\n\n*** valid activation energy\n"
                                     f" basin 1 energy      : {basin1_energy.get_total_energy()}\n"
                                     f" basin 2 energy      : {basin2_energy.get_total_energy()}\n"
                                     f" saddle point energy : {saddle_energy.get_total_energy()}\n"
                                     f" activation energy   : {act_energy}\n"
                                     f" centre              : {centre[:]}\n"
                                     f" search number       : {search_number}\n")

                            finish_time = time.time()
                            elapsed_time = finish_time - start_time
                            log_stream.write(f" time to calculate event : {elapsed_time:.3f} seconds\n")

                
                
                        else:
                            log_stream.write(f"\n invalid activation energy {act_energy}\n")
                            status.set_status_failed()

                    #successful or otherwise write out the status for use later
                    status.write(new_path)
                    log_stream.close()
                    del searcher

                    #send an event back to parent process - always rank 0
                    #an event consists of basin1, basin2, saddle + act energies and search number
                    #also included is the status flag
                    results = np.zeros(9)
                    results[0] = float(status.get_status())
                    results[1] = float(search_number) + 0.1
                    if status.get_status() == SUCCESS:
                        results[2] = basin1_energy.get_total_energy()
                        results[3] = basin2_energy.get_total_energy()
                        results[4] = saddle_energy.get_total_energy()
                        results[5] = act_energy
                        results[6] = centre[0]
                        results[7] = centre[1]
                        results[8] = centre[2]
                    else:
                        results[2] = 0.0
                        results[3] = 0.0
                        results[4] = 0.0
                        results[5] = 0.0
                        results[6] = 0.0
                        results[7] = 0.0
                        results[8] = 0.0

                    #now send
                    self.send_message_to_parent(results)


        
    def calculate_rate_constant_experimental(self, spec, fld, job: JobControl, basin1, basin1_energy: Energy, event_list, excluded_list, iteration, kmc_time, out_stream):
        comm = MPI.COMM_WORLD
        world_rank = comm.Get_rank()
        num_procs = comm.Get_size()

        if world_rank == 0:
            num_events = len(event_list)
            search_number = -1

            if job.kmc_parameters.recycle_events == True:
                search_number = len(event_list) - 1
                print("the starting search number is")

        act_energy = 0.0

        basin2_energy = Energy()
        saddle_energy = Energy()
        
        continue_search = False

        searcher = None
        status = Status()
        #ev = Event()
        attempted_search = 0
        successfull_search = 0
        num_busy = 0

        it_path = "iteration" + str(iteration)
        if world_rank == 0:
            if job.kmc_parameters.recycle_events == False:
                try:
                    os.mkdir(it_path)
                except Exception as e:
                    out_stream.write(f"{e}\n")
                    out_stream.flush()
                    traceback.print_exc()
                    MPI.COMM_WORLD.Abort()
                    exit()
            elif iteration == 1:
                try:
                    os.mkdir(it_path)
                except Exception as e:
                    out_stream.write(f"{e}\n")
                    out_stream.flush()
                    traceback.print_exc()
                    MPI.COMM_WORLD.Abort()
                    exit()

            basin1.config_write(it_path, "basin1.xyz", basin1_energy.get_total_energy(), iteration, kmc_time) # basin1 only needs to be written out once
        
        comm.Barrier()  # dont let processors get ahead and create a search directory before iteration

        if world_rank == 0:  ################################ parent process to initiate and keep track of events
            num_search = num_procs
            i_am_busy = np.zeros(num_procs, dtype=int)
            
            #send a command to all groups to start search
            for i in range(1, num_search):
                search_number += 1
                attempted_search +=1
                num_busy += 1
                i_am_busy[i] = 1
                command_to_children = np.zeros(2, dtype=np.int32)
                command_to_children[0] = search_number
                command_to_children[1] = iteration
                out_stream.write(f"start  {i} {search_number} \n")
                out_stream.flush()
                self.send_command_to_child(command_to_children, i)


            while num_events < job.kmc_parameters.kmc_images:

                #wait to receive an event
                command, source, tag = self.check_for_message()
                i_am_busy[source] = 0
                num_busy -= 1
                #out_stream.write(f"command received by parent {source} {tag} {command[0]} {int(command[1])}  {np.rint(command[1])}\n")
                #out_stream.flush()

                #check to see whether it is a valid event
                if int(command[0]) >= 0:

                    recvd_iteration = int(command[2])
                    if recvd_iteration == iteration:
                        id = int(command[1])
                    
                        basin2_energy.vdwEnergy = np.float64(command[4])
                        saddle_energy.vdwEnergy = np.float64(command[5])
                        act_energy = np.float64(command[6])
                        centre = np.zeros(3, dtype=np.float64)
                        centre[0] = np.float64(command[7])
                        centre[1] = np.float64(command[8])
                        centre[2] = np.float64(command[9])
                        #out_stream.write(f"data received by parent {source} {tag} {command[:]}\n")
                        #out_stream.flush()
                        ev = Event(basin1_energy.get_total_energy(), basin2_energy.get_total_energy(), saddle_energy.get_total_energy(), act_energy, centre, id)
                        ev.calculate_rate(job.kmc_parameters.kmcTemperature, job.kmc_parameters.kmcPreFactor)
                        event_list.append(ev)
                        num_events += 1
                        successfull_search += 1
                        #save the events
                        eventname = os.path.join(it_path, "events")
                        fp = open(eventname, "a")
                        #for e in event_list:   #make a copy of all valid events for restart
                        ev.write_event(fp)
                        fp.flush()
                        fp.close()
                    else:
                        out_stream.write(f"WARNING:received data from the incorrect iteration {source} {iteration} {recvd_iteration} {command[:]}\n")
                        out_stream.flush()

                
                #check to see if enough searches have been launched
                if num_events < job.kmc_parameters.kmc_images:
                    #send the search number to the child to initiate search again
                    search_number += 1
                    attempted_search +=1
                    num_busy += 1
                    i_am_busy[source] = 1
                    #out_stream.write(f"restart search  {num_events} {job.kmc_parameters.kmc_images} { source} {search_number} \n")
                    #out_stream.flush()
                    command_to_children = np.zeros(2, dtype=np.int32)
                    command_to_children[0] = search_number
                    command_to_children[1] = iteration
                    self.send_command_to_child(command_to_children, source)

            #have all the number of events so send signal to children for return
            
            self.send_signal_to_children()
                
            out_stream.write(f"\n\n attempted searches : {attempted_search} successfull searches : {successfull_search}")
        else:   ######################################### child processors
            command = np.zeros(2, dtype=np.int32)
            command[0] = 0
            while command[0] >= 0:
                command, source, tag = self.check_for_message()
                search_number = command[0]
                sent_iteration = command[1]
                if search_number < 0:
                    break

                if search_number >= 0:
                    
                    flnm = "search" + str(search_number)
                    #create a new directory for this iteration and attemp for the search
                    new_path = os.path.join(it_path, flnm)
                    try:
                        os.mkdir(new_path)
                    except Exception as e:
                        sys.stderr.write(f"{e}\n")
                        sys.stderr.flush()
                        traceback.print_exc()
                        MPI.COMM_WORLD.Abort()
                        exit()
                    #create a stream for the log file to go into this directory
                    log_path = os.path.join(new_path, "logfile")
                    try:
                        log_stream = open(log_path, "w")
                    except Exception as e:
                        sys.stderr.write(f"{e}\n")
                        sys.stderr.flush()
                        traceback.print_exc()
                        MPI.COMM_WORLD.Abort()
                        exit()

                    log_stream.write(f"\n\n started new search {search_number} \n")

                    start_time = time.time()
                    attempted_search += 1
                    searcher = Saddle(iteration, search_number, job.kmc_parameters.kmc_method)
                    #try:
                    basin2_energy, centre, status = searcher.find_saddle_point(spec, job, fld, basin1, saddle_energy, continue_search, new_path, log_stream)  
                    log_stream.write(f"\n\n finished saddle search {status.get_status()}\n")
                    #except Exception as e:
                    #    log_stream.write("\n kinetic_monte_carlo : ")
                    #    log_stream.write(f"{e}\n")
                    #    status.set_status_failed()
                    new_command, new_source, new_tag = self.noninterrupt_check_for_message()
                    if new_command == status.TERMINATED:
                        status.set_status_terminated()
                

                    if status.get_status() == Status.TERMINATED:
                        log_stream.write("\n search terminated \n")
                        command[0] = status.TERMINATED

                    if status.get_status() == Status.FAILED:
                        log_stream.write("\n search failed\n")
                        

                    
                    if status.get_status() == Status.SUCCESS:
                        act_energy = saddle_energy.get_total_energy() - basin1_energy.get_total_energy()

                        if act_energy >= job.kmc_parameters.kmc_min_cap and act_energy <= job.kmc_parameters.kmc_max_cap:
                            log_stream.write(f"\n\n*** valid activation energy\n"
                                     f" basin 1 energy      : {basin1_energy.get_total_energy()}\n"
                                     f" basin 2 energy      : {basin2_energy.get_total_energy()}\n"
                                     f" saddle point energy : {saddle_energy.get_total_energy()}\n"
                                     f" activation energy   : {act_energy}\n"
                                     f" centre              : {centre[:]}\n"
                                     f" search number       : {search_number}\n")

                            finish_time = time.time()
                            elapsed_time = finish_time - start_time
                            log_stream.write(f" time to calculate event : {elapsed_time:.3f} seconds\n")

                
                
                        else:
                            log_stream.write(f"\n invalid activation energy {act_energy}\n")
                            status.set_status_failed()

                    #successful or otherwise write out the status for use later
                    status.write(new_path)
                    log_stream.close()
                    del searcher

                    #send an event back to parent process - always rank 0
                    #an event consists of basin1, basin2, saddle + act energies and search number
                    #also included is the status flag
                    if status.get_status() != Status.TERMINATED:
                        results = np.zeros(10, dtype=np.float64)
                        results[0] = np.float64(status.get_status())
                        results[1] = np.float64(search_number) + 0.1
                        results[2] = np.float64(sent_iteration) + 0.1
                        if status.get_status() == SUCCESS:
                            results[3] = np.float64(basin1_energy.get_total_energy())
                            results[4] = np.float64(basin2_energy.get_total_energy())
                            results[5] = np.float64(saddle_energy.get_total_energy())
                            results[6] = np.float64(act_energy)
                            results[7] = np.float64(centre[0])
                            results[8] = np.float64(centre[1])
                            results[9] = np.float64(centre[2])

                        #now send
                        self.send_message_to_parent(results)


        
    def calculate_rate_constant_single_core(self, spec, fld, job: JobControl, basin1: Config, basin1_energy: Energy, event_list, iteration, kmc_time, out_stream):
        num_events = len(event_list)
        search_number = -1
        attempted_search = 0
        successfull_search = 0

        if job.kmc_parameters.recycle_events == True:
            search_number = len(event_list) - 1

        act_energy = 0.0

        basin2_energy = Energy()
        saddle_energy = Energy()
        
        continue_search = False

        searcher = None
        status = Status()

        #create a folder for each iteration
        it_path = "iteration" + str(iteration)
        if job.kmc_parameters.recycle_events == False:
            try:
                os.mkdir(it_path)
            except Exception as e:
                out_stream.write(f"{e}\n")
                out_stream.flush()
                traceback.print_exc()
                MPI.COMM_WORLD.Abort()
                exit()
        elif iteration == 1:
            try:
                os.mkdir(it_path)
            except Exception as e:
                out_stream.write(f"{e}\n")
                out_stream.flush()
                traceback.print_exc()
                MPI.COMM_WORLD.Abort()
                exit()

        basin1.config_write(it_path, "basin1.xyz", basin1_energy.get_total_energy(), iteration, kmc_time) # basin1 only needs to be written out once

        while num_events < job.kmc_parameters.kmc_images:

            search_number += 1

            flnm = "search" + str(search_number)
            #create a new directory for this iteration and attemp for the search
            new_path = os.path.join(it_path, flnm)
            os.mkdir(new_path)

            #create a stream for the log file to go into this directory
            log_path = os.path.join(new_path, "logfile")
            try:
                log_stream = open(log_path, "w")
            except Exception as e:
                out_stream.write(f"{e}\n")
                out_stream.flush()
                traceback.print_exc()
                MPI.COMM_WORLD.Abort()
                exit()


            start_time = time.time()
            attempted_search += 1
            searcher = Saddle(iteration, num_events, job.kmc_parameters.kmc_method)
            #try:
            
            basin2_energy, centre, status = searcher.find_saddle_point(spec, job, fld, basin1, saddle_energy, continue_search, new_path, log_stream)  
            log_stream.write(f"\n\n finished saddle search {status.get_status()}\n")
            #except Exception as e:
            #    log_stream.write("\n kinetic_monte_carlo : ")
            #    log_stream.write(f"{e}\n")
            #    status.set_status_failed()

            if status.get_status() == Status.FAILED or status.get_status() == Status.TERMINATED:
                log_stream.write("\n search failed\n")
                continue

            if status.get_status() == Status.FAILED:
                continue

            act_energy = saddle_energy.get_total_energy() - basin1_energy.get_total_energy()

            if act_energy >= job.kmc_parameters.kmc_min_cap and act_energy <= job.kmc_parameters.kmc_max_cap:
                log_stream.write(f"\n\n*** valid activation energy\n"
                                     f" basin 1 energy      : {basin1_energy.get_total_energy()}\n"
                                     f" basin 2 energy      : {basin2_energy.get_total_energy()}\n"
                                     f" saddle point energy : {saddle_energy.get_total_energy()}\n"
                                     f" activation energy   : {act_energy}\n"
                                     f" centre              : {centre[:]}\n")

                finish_time = time.time()
                elapsed_time = finish_time - start_time
                log_stream.write(f" time to calculate event : {elapsed_time:.3f} seconds\n")

                
                ev = Event(basin1_energy.get_total_energy(), basin2_energy.get_total_energy(), saddle_energy.get_total_energy(), act_energy, centre, search_number)
                ev.calculate_rate(job.kmc_parameters.kmcTemperature, job.kmc_parameters.kmcPreFactor)
                event_list.append(ev)
                num_events += 1
                successfull_search += 1
                #save the events
                filename = os.path.join(it_path,"events")
                fp = open(filename, "w")
                for e in event_list:   #make a copy of all valid events for restart
                    e.write_event(fp)
                fp.flush()
                fp.close()
                
             
            else:
                log_stream.write(f"\n invalid activation energy {act_energy}\n")
                status.set_status_failed()
            
            #successful or otherwise write out the status for use later
            status.write(new_path)
            log_stream.close()
            del searcher

        out_stream.write(f"\n\n attempted searches : {attempted_search} successfull searches : {successfull_search}")

    @staticmethod
    def calculate_time(total_rate):
        ran_number = random.random()
        if abs(total_rate) > 1.0e-18:
            return -np.log(ran_number) / total_rate
        return 0.0
    
    def send_signal_to_children(self):
        comm = MPI.COMM_WORLD
        dummy = np.zeros(2, dtype = np.int32)
        dummy[0] = TERMINATED
        tag = 999
        for i in range(1, comm.Get_size()):
            #print("terminate signal", i)
            comm.send(dummy, dest=i, tag=tag)
        return 
    
    def send_command_to_child(self, command, child):
        comm = MPI.COMM_WORLD
        tag = 999
        comm.send(command, dest=child, tag=tag)

    def send_message_to_parent(self, command):
        comm = MPI.COMM_WORLD
        tag = 999
        comm.send(command, dest=0, tag=tag)

    def check_for_message(self):
        comm = MPI.COMM_WORLD
        status = MPI.Status()
        command = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        if status.Get_error() != MPI.SUCCESS:
            #print(f"\n MPI::check_for message error {status.Get_error()} {status.Get_error_class()} rank {comm.Get_rank()} \n")
            traceback.print_exc()
            MPI.COMM_WORLD.Abort(EXIT_FAILURE)
            exit(EXIT_FAILURE)
        return command, status.Get_source(), status.Get_tag()
    
    def noninterrupt_check_for_message(self):
        comm = MPI.COMM_WORLD
        status = MPI.Status()
        source = 0
        command = np.zeros(2, dtype=np.int32)

        flag = comm.Iprobe(source, tag=MPI.ANY_TAG, status=status)
        if flag:
            command = comm.recv(source=source, tag=MPI.ANY_TAG, status=status)

        if status.Get_error() != MPI.SUCCESS:
            print(f"\n ART::check_for message error {status.Get_error()} {status.Get_error_class()} rank {comm.Get_rank()} \n")
            traceback.print_exc()
            MPI.COMM_WORLD.Abort(EXIT_FAILURE)
            exit(EXIT_FAILURE)

        return command[0], status.Get_source(), status.Get_tag()

    
    def broadcast_positions(self, basin: Config):
        comm = MPI.COMM_WORLD
        #tag = 999
        buffer = np.zeros((basin.get_number_of_atoms(), 3))
        buffer[:,:] = basin.pos[:,:]
        comm.Bcast(buffer, root=0)
        basin.pos[:,:] = buffer[:,:]

    def recycle_saddle_search(self, event_list, new_basin, new_basin_energy, selected_event, active_distance, recycle_distance, iteration, min_cap, probability, create_dir, outstream):
        #stub for recycling - possiblt wont be good for ionic materials where polarisation of the atoms/ioons
        #extends over large distance
        radius = recycle_distance * recycle_distance
        new_search_number = 0
        num_images = len(event_list)
        new_event_list = []

        #create the new directory for updated searches
        it_path = "iteration" + str(iteration+1)
        outstream.write(f"\n new iteration being created {it_path} \n")
        outstream.flush()
        if create_dir:
            try:
                os.mkdir(it_path)
            except Exception as e:
                outstream.write(f"{e}\n")
                outstream.flush()
                traceback.print_exc()
                MPI.COMM_WORLD.Abort()
                exit()
        
        #obtain the centre of the new event the selected basin 
        select_centre = event_list[selected_event].get_centre()

        #make a list of all atoms within the active region of the new basin
        num_list, px, py, pz = new_basin.find_active_region(select_centre, active_distance)
        
        active_size = len(num_list)

        #try to find a basin outside of the active region and then update positions withing the active region
        latvector = new_basin.vectors
        rcpvec = np.linalg.inv(latvector)
        for i in range(num_images):

            if i == selected_event:
                continue

            #get the centre of the current event
            ev_centre = event_list[i].get_centre()
            #print("centre of event", i, ev_centre)
            #check the
            rx = select_centre[0] - ev_centre[0]
            ry = select_centre[1] - ev_centre[1]
            rz = select_centre[2] - ev_centre[2]

            #correct for PBC
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
            #print("distance ", np.sqrt(rsq))
            if rsq >= radius:
                #steps
                #1. read in the saddle point and basin2 corresponding to the event
                #2. replac the positions of the saddle and basin2 
                #3. put the new structures in a new directory for the next iteration 
                #4. update the evnt list
            
                search_number = event_list[i].get_search_number()
                tmp = "iteration" + str(iteration)
                tmp1 = "search" + str(search_number)
                file_name = os.path.join(tmp,tmp1,"saddle.xyz")
                
                #outstream.write(f"\n resetting basin1 from file  {file_name} \n")
                saddle = Config()
                #print("opening files in ", file_name)
                if os.path.exists(file_name) == False:
                    outstream.write(f"failed to find saddle file {file_name} \n")
                    continue

                with open(file_name, "r") as instream:
                    dummy, dumm1, saddle_energy = saddle.read_config(instream)
                file_name = os.path.join(tmp,tmp1,"basin2.xyz")
                #outstream.write(f"\n resetting basin1 from file  {file_name} \n")
                basin2 = Config()
                with open(file_name, "r") as instream:
                    dummy, dumm1, basin2_energy = basin2.read_config(instream)

                #check that the activation energy is still valid
                act_energy = saddle_energy - new_basin_energy

                if act_energy < min_cap:
                    continue

                if np.random.random() > probability:
                    continue

                #replace positions
                for j in range(len(num_list)):
                    atm = num_list[j]
                    saddle.pos[atm][0] = px[j]
                    basin2.pos[atm][0] = px[j]
                    saddle.pos[atm][1] = py[j]
                    basin2.pos[atm][1] = py[j]
                    saddle.pos[atm][2] = pz[j]
                    basin2.pos[atm][2] = pz[j]

                

                #save the updated positions in the new directory
                
                tmp1 = "search" + str(new_search_number)
                new_path = os.path.join(it_path, tmp1)
                try:
                    os.mkdir(new_path)
                except Exception as e:
                    outstream.write(f"{e}\n")
                    outstream.flush()
                    MPI.COMM_WORLD.Abort()
                    exit()
                basin2.config_write(new_path, "basin2.xyz", basin2_energy)   
                saddle.config_write(new_path, "saddle.xyz", saddle_energy) 

                #update the event, write out and append to the new list of events
               
                event_list[i].set_activation_energy(act_energy)
                event_list[i].setBasin1Energy(new_basin_energy)
                event_list[i].set_search_number(new_search_number)
                filename = os.path.join(it_path,"events")
                fs = open(filename,"a")
                event_list[i].write_event(fs)
                fs.close()
                new_event_list.append(event_list[i])

                new_search_number = new_search_number + 1


        outstream.write(f"\n the number of recycled saddle points {len(new_event_list)} \n")
        outstream.flush()

        del event_list[:]
        return new_event_list
    
    def recycle_saddle_limit(self, event_list, new_basin, new_basin_energy, selected_event, active_distance, recycle_distance, iteration, min_cap, act_limit, create_dir, outstream):
        #stub for recycling - possiblt wont be good for ionic materials where polarisation of the atoms/ioons
        #extends over large distance
        radius = recycle_distance * recycle_distance
        new_search_number = 0
        num_images = len(event_list)
        new_event_list = []

        #create the new directory for updated searches
        it_path = "iteration" + str(iteration+1)
        outstream.write(f"\n new iteration being created {it_path} \n")
        outstream.flush()
        if create_dir:
            try:
                os.mkdir(it_path)
            except Exception as e:
                outstream.write(f"{e}\n")
                outstream.flush()
                traceback.print_exc()
                MPI.COMM_WORLD.Abort()
                exit()
        
        #obtain the centre of the new event the selected basin 
        select_centre = event_list[selected_event].get_centre()

        
        #make a list of all atoms within the active region of the new basin
        num_list, px, py, pz = new_basin.find_active_region(select_centre, active_distance)
        
        active_size = len(num_list)

        #try to find a basin outside of the active region and then update positions withing the active region
        latvector = new_basin.vectors
        rcpvec = np.linalg.inv(latvector)
        for i in range(num_images):

            if i == selected_event:
                continue

            #get the centre of the current event
            ev_centre = event_list[i].get_centre()
            #print("centre of event", i, ev_centre)
            #check the
            rx = select_centre[0] - ev_centre[0]
            ry = select_centre[1] - ev_centre[1]
            rz = select_centre[2] - ev_centre[2]

            #correct for PBC
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
            #print("distance ", np.sqrt(rsq))
            if rsq >= radius:
                #steps
                #1. read in the saddle point and basin2 corresponding to the event
                #2. replac the positions of the saddle and basin2 
                #3. put the new structures in a new directory for the next iteration 
                #4. update the evnt list
            
                search_number = event_list[i].get_search_number()
                tmp = "iteration" + str(iteration)
                tmp1 = "search" + str(search_number)
                file_name = os.path.join(tmp,tmp1,"saddle.xyz")
                #outstream.write(f"\n resetting basin1 from file  {file_name} \n")
                saddle = Config()
                #print("opening files in ", file_name)
                with open(file_name, "r") as instream:
                    dummy, dumm1, saddle_energy = saddle.read_config(instream)
                file_name = os.path.join(tmp,tmp1,"basin2.xyz")
                #outstream.write(f"\n resetting basin1 from file  {file_name} \n")
                basin2 = Config()
                with open(file_name, "r") as instream:
                    dummy, dumm1, basin2_energy = basin2.read_config(instream)

                #check that the activation energy is still valid
                act_energy = saddle_energy - new_basin_energy

                if act_energy < min_cap:
                    continue

                if act_energy > act_limit:
                    continue

                #replace positions
                for j in range(len(num_list)):
                    atm = num_list[j]
                    saddle.pos[atm][0] = px[j]
                    basin2.pos[atm][0] = px[j]
                    saddle.pos[atm][1] = py[j]
                    basin2.pos[atm][1] = py[j]
                    saddle.pos[atm][2] = pz[j]
                    basin2.pos[atm][2] = pz[j]

                

                #save the updated positions in the new directory
                
                tmp1 = "search" + str(new_search_number)
                new_path = os.path.join(it_path, tmp1)
                try:
                    os.mkdir(new_path)
                except Exception as e:
                    outstream.write(f"{e}\n")
                    outstream.flush()
                    MPI.COMM_WORLD.Abort()
                    exit()
                basin2.config_write(new_path, "basin2.xyz", basin2_energy)   
                saddle.config_write(new_path, "saddle.xyz", saddle_energy) 

                #update the event, write out and append to the new list of events
               
                event_list[i].set_activation_energy(act_energy)
                event_list[i].setBasin1Energy(new_basin_energy)
                event_list[i].set_search_number(new_search_number)
                filename = os.path.join(it_path,"events")
                fs = open(filename,"a")
                event_list[i].write_event(fs)
                fs.close()
                new_event_list.append(event_list[i])

                new_search_number = new_search_number + 1


        outstream.write(f"\n the number of recycled saddle points {len(new_event_list)} \n")
        outstream.flush()

        del event_list[:]
        return new_event_list
    

    def recover_event_list(self, iteration, temperature, pre_factor, outstream):

        new_list = []

        try:
            it = "iteration" + str(iteration)
            filename = os.path.join(it, "events")
            instream = open(filename, 'r')
            lines = instream.readlines()
                
            idx = 0
            
            for _ in range(len(lines) // 8):
                line = lines[idx]
                
                words = line.split()
                
                if "EVENT" in words[0]:
                    centre = np.zeros(3)
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    basin_1_energy = float(words[1])
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    basin_2_energy = float(words[1])
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    transition_state_energy = float(words[1])
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    activation_energy = float(words[1])
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    rate  = float(words[1])
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    centre[0] = float(words[1])   
                    centre[1] = float(words[2])
                    centre[2] = float(words[3])
                    idx += 1
                    line = lines[idx]
                    words = line.split()
                    image_number  = int(words[1])

                    ev = Event(basin_1_energy, basin_2_energy, transition_state_energy, activation_energy, centre, image_number)
                    ev.calculate_rate(temperature, pre_factor)

                    new_list.append(ev)
                
                idx += 1
        except FileNotFoundError:
            outstream.write("\n*** could not find the events file\n")
            outstream.flush()
            traceback.print_exc()
            MPI.COMM_WORLD.Abort(EXIT_FAILURE)
            exit(EXIT_FAILURE)

        outstream.write(f"\n the number of recovered events {len(new_list)} \n")
        outstream.flush()
        return new_list
    
    def read_selected_event(self, old_iteration, outstream):
        it_path = "iteration" + str(old_iteration)
        filename = os.path.join(it_path,"selection")

        try:
            instream = open(filename, 'r')
            line = instream.readline()
           
            words = line.split()
            ievent = int(words[3])
        except FileNotFoundError:
            outstream.write("\n*** could not find the selected event file: {filename}\n")
            outstream.flush()
            traceback.print_exc()
            MPI.COMM_WORLD.Abort(EXIT_FAILURE)
            exit(EXIT_FAILURE)

        return ievent

                
