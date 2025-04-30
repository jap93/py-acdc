import os

import numpy as np
from scipy.linalg import eigh_tridiagonal, norm
from mpi4py import MPI

from energy import Energy
from job_control import JobControl
from search_control import SearchCntrl
from config import Config
from status import Status
from field import Field
from constants import EVTODL, PI

class ART(object):

    def __init__(self):
        self.vector_size = 0
        self.r_vector: np.float64 = None
        pass

    def _initialise(self, frozen, natoms, num_frozen, indeces, mask):
        self.vector_size = (natoms - num_frozen) * 3

    
    def run(self, saddle: Config, r_vector, displacement_vector, saddle_energy: Energy, fld: Field, job: JobControl, saddle_status: Status, search_peprpindicular: bool, out_stream):

        saddle_status.set_status_failed()
        indices = []
        mask = []
        self._initialise(saddle.frozen, saddle.natoms, saddle.num_frozen, indices, mask)

        ase_search = False
        flag = (os.environ.get("USE_ASE_SEARCH", ase_search) == "True")

        if flag == True:  # use ase dimer method
            #raise(NotImplementedError)
            out_stream.write("\n not implemented yet")
            out_stream.flush()
            exit(-1)
        else:
            #run search for the saddle point
            if job.search_parameters.min_method == "fire":
                self.fire_saddle_point(saddle, r_vector , saddle_energy, fld, job, saddle_status, search_peprpindicular, out_stream)
            elif job.search_parameters.min_method == "fire2":
                self.fire2_saddle_point(saddle, r_vector , saddle_energy, fld, job, saddle_status, search_peprpindicular, out_stream)
            else:
                out_stream.write("ART minimisation method not defined")
                exit(-1)
    

    def fire_saddle_point(self, bas: Config, r_vector, final_energy: Energy, fld, job: JobControl, status: Status, search_peprpindicular: bool, out_stream):

        max_iterations = 0
        n_step = 1

        alpha = job.search_parameters.alpha_start
        eigen_value = 0.0
        delta_e = 0.0

        total_energy_new = Energy()
        total_energy_old = Energy()

        converged_fire = False
        converged_lanczos = False

        status.set_status_failed()

        self.vector_size = len(r_vector)
    
        # Allocate work arrays
        velocity = np.zeros((int(bas.natoms),3), dtype=np.float64)
        self.r_vector = np.zeros(self.vector_size, dtype=np.float64)
        np.copyto(self.r_vector, r_vector)
        
       # try:
        force, eigen_value, converged_lanczos = self.calculate_lanczos_force(job, bas, fld, total_energy_old, search_peprpindicular, out_stream)

        # Set FIRE parameters
        if search_peprpindicular:
            max_iterations = job.search_parameters.max_perp_iter
            time_step = job.search_parameters.init_time_step_perp
            time_step_max = job.search_parameters.max_time_step_perp

            if eigen_value < job.search_parameters.max_eigenvalue:  #there is already a negative eigenvalue
                status.set_status_success()
                return
        else:
            max_iterations = job.search_parameters.max_par_iter
            time_step = job.search_parameters.init_time_step
            time_step_max = job.search_parameters.max_time_step
            

        if converged_lanczos == False:
            status.set_status_failed()
        
        if job.search_parameters.debug:
            out_stream.write(f"\nFIRE Saddle Search: initial energy {total_energy_old.get_total_energy()}\n")
            out_stream.write(f"\nFIRE Saddle Search: initial eigenvalue {eigen_value}\n")
            out_stream.flush()

        half_step = 0.5 * time_step * EVTODL

        it_art = 0
        while not converged_fire and it_art <= max_iterations:
            it_art += 1

            for i in range(bas.natoms):
                for j in range(3):
                    velocity[i,j] += half_step * force[i,j]

            f_norm, f_max = self.vector_norm_max(force)
            v_norm, v_max = self.vector_norm_max(velocity)

            P = np.dot(force.flatten(), velocity.flatten())
            OneLessAlpha = 1.0 - alpha

            for i in range(bas.natoms):
                for j in range(3):
                    velocity[i,j] = OneLessAlpha * velocity[i,j] + alpha * force[i,j] * v_norm / f_norm

            if P > 0:  # downhill move
                if n_step > job.search_parameters.N_min:
                    time_step = min(time_step * job.search_parameters.step_increase, time_step_max)
                    alpha *= job.search_parameters.alpha_decrease
                n_step += 1
            else:  # uphill move
                n_step = 0
                time_step *= job.search_parameters.step_decrease
                alpha = job.search_parameters.alpha_start
                for i in range(bas.natoms):
                    for j in range(3):
                        velocity[i,j] = 0.0

            # Advance coordinates
            max_displacement = -1.0e10
            for i in range(bas.natoms):
                if bas.frozen[i] == 0:
                    for j in range(3):
                        bas.pos[i][j] = bas.pos[i][j] + time_step * velocity[i][j]
                    
                    #bas.pos[i][:] += dpos[:]
                    #delta = np.sqrt(dx**2 + dy**2 + dz**2)
                    #if delta > max_displacement:
                    #    max_displacement = delta

            force, lowest_eigenvalue, converged_lanczos = self.calculate_lanczos_force(job, bas, fld, total_energy_new, search_peprpindicular, out_stream)

            velocity += half_step * force

            # Check convergence
            delta_e = abs(total_energy_old.get_total_energy() - total_energy_new.get_total_energy())

            if job.search_parameters.debug:
                out_stream.write(f"\n FIRE_ART: iteration {it_art} {total_energy_new.get_total_energy()} eigenvalue {lowest_eigenvalue} time_step {time_step}\n")
                out_stream.write(f" FIRE_ART: gnorm {f_norm} force max {f_max} energy difference {delta_e}\n")
                out_stream.flush()

            
            
            if it_art > 5:
                if lowest_eigenvalue < 0.0:
                    if (f_norm) < job.search_parameters.norm_tol or (f_max) < job.search_parameters.force_tol:
                        if job.search_parameters.debug:
                            out_stream.write(f"\nFIRE converged: delta E = {delta_e}\n")
                            #out_stream.write(f"\n lowest eigenvalue = {lowest_eigenvalue}")
                            converged_fire = True
                        if lowest_eigenvalue <= job.search_parameters.max_eigenvalue:
                            status.set_status_success()
                        else:
                            status.set_status_failed()

                        total_energy_new.copy_to(final_energy)
                        return

            half_step = 0.5 * time_step * EVTODL
            total_energy_new.copy_to(total_energy_old)

            if search_peprpindicular and lowest_eigenvalue < job.search_parameters.max_eigenvalue:  #there is a negative eigenvalue in perp search
                status.set_status_success()
                total_energy_new.copy_to(final_energy)
                out_stream.write(f"\n FIRE_ART: perpendicular iteration {it_art} a negative eigenvalue has been found {lowest_eigenvalue}\n")
                return
            
            if search_peprpindicular == False and job.search_parameters.stop_pos_eigenval == True:
                if lowest_eigenvalue > job.search_parameters.max_eigenvalue:
                    status.set_status_failed()
                    total_energy_new.copy_to(final_energy)
                    out_stream.write(f"\n FIRE_ART: iteration {it_art} positive eigenvalue has been found {lowest_eigenvalue}\n")
                    return
                
            if np.abs(lowest_eigenvalue) > 500.0:
                status.set_status_failed()
                total_energy_new.copy_to(final_energy)
                out_stream.write(f"\n FIRE_ART: iteration {it_art} something is going horribly wrong {lowest_eigenvalue}\n")
                return
            
            if abs(half_step) < 1.0e-10:
                status.set_status_failed()
                return

        total_energy_new.copy_to(final_energy)
        # Too many iterations
        status.set_status_failed()

    def fire2_saddle_point(self, bas: Config, r_vector, final_energy: Energy, fld, job: JobControl, status: Status, search_peprpindicular: bool, out_stream):
        
        nUpHill = 0
        nStep = 1
        alpha = job.search_parameters.alpha_start

        total_energy_new = Energy()
        total_energy_old = Energy()
        velocity = np.zeros((bas.natoms, 3), dtype=np.float64)
        force_new = np.zeros((bas.natoms, 3), dtype=np.float64)

        self.r_vector = np.zeros(self.vector_size, dtype=np.float64)
        np.copyto(self.r_vector, r_vector)
        
       # try:
        force_new, eigen_value, converged_lanczos = self.calculate_lanczos_force(job, bas, fld, total_energy_old, search_peprpindicular, out_stream)

        command, source, tag = self.check_for_message()
        if command == status.TERMINATED:
            status.set_status_terminated()
            return
        
        # Set FIRE parameters
        time_step_min = job.search_parameters.min_time_step
        if search_peprpindicular:
            max_iterations = job.search_parameters.max_perp_iter
            time_step = job.search_parameters.init_time_step_perp
            time_step_max = job.search_parameters.max_time_step_perp

            if eigen_value < job.search_parameters.max_eigenvalue:  #there is already a negative eigenvalue
                status.set_status_success()
                return
        else:
            max_iterations = job.search_parameters.max_par_iter
            time_step = job.search_parameters.init_time_step
            time_step_max = job.search_parameters.max_time_step
            

        if converged_lanczos == False:
            status.set_status_failed()        

        halfStep = 0.5 * time_step * EVTODL

        for it in range(1, max_iterations):
            command, source, tag = self.check_for_message()
            if command == status.TERMINATED:
                status.set_status_terminated()
                return

            P = np.dot(force_new.flatten(), velocity.flatten())

            OneLessAlpha = 1.0 - alpha
            fNorm, fMax = self.vector_norm_max(force_new)
            vNorm, vMax = self.vector_norm_max(velocity)

            alphaForce = alpha * vNorm / fNorm

            if P > 0:
                velocity = OneLessAlpha * velocity + alphaForce * force_new

                if nStep > job.search_parameters.N_min:
                    time_step = min(time_step * job.search_parameters.step_increase, time_step_max)
                    alpha *= job.search_parameters.alpha_decrease

                nStep += 1
                nUpHill = 0
            else:
                nStep = 0
                nUpHill += 1

                if nUpHill > job.search_parameters.N_max_uphill:
                    status.set_status_failed()
                    return

                if job.search_parameters.initial_delay:
                    if it >= job.search_parameters.N_delay_steps:
                        time_step *= job.search_parameters.step_decrease
                        if time_step < time_step_min:
                            time_step = time_step_min
                        alpha = job.search_parameters.alpha_start
                else:
                    time_step *= job.search_parameters.step_decrease
                    alpha = job.search_parameters.alpha_start

                bas.pos -= 0.5 * time_step * velocity

                velocity.fill(0.0)

            if job.search_parameters.verlet_integration == False:
                fact = time_step * EVTODL
                for i in range(bas.natoms):
                    if bas.frozen[i] == 0:
                        if P > 0.0:
                            velocity[i][:] = OneLessAlpha * velocity[i][:] + alphaForce * force_new[i][:]
                        
                        bas.pos[i][:] += time_step * velocity[i][:]
                    
                        velocity[i][:] += fact * force_new[i][:]
                    
            else:
                halfStep = 0.5 * time_step * EVTODL
                for i in range(bas.natoms):
                    if bas.frozen[i] == 0:
                        if P > 0.0:
                            velocity[i][:] = OneLessAlpha * velocity[i][:] + alphaForce * force_new[i][:]
                        
                        velocity[i][:] += halfStep * force_new[i][:]
                    
                bas.pos += time_step * velocity
            
            force_new, lowest_eigenvalue, converged_lanczos = self.calculate_lanczos_force(job, bas, fld, total_energy_new, search_peprpindicular, out_stream)

            if job.search_parameters.verlet_integration:
                for i in range(bas.natoms):
                    if bas.frozen[i] == 0:
                        velocity[i][:] += halfStep * force_new[i][:]
                    
            #checks
            if search_peprpindicular and lowest_eigenvalue < job.search_parameters.max_eigenvalue:  #there is a negative eigenvalue in perp search
                status.set_status_success()
                total_energy_new.copy_to(final_energy)
                out_stream.write(f"\n FIRE2_ART: perpendicular iteration {it} a negative eigenvalue has been found {lowest_eigenvalue}\n")
                return
            
            if search_peprpindicular == False and job.search_parameters.stop_pos_eigenval == True:
                if lowest_eigenvalue > job.search_parameters.max_eigenvalue:
                    status.set_status_failed()
                    total_energy_new.copy_to(final_energy)
                    out_stream.write(f"\n FIRE2_ART: iteration {it} positive eigenvalue has been found {lowest_eigenvalue}\n")
                    return
                
            if np.abs(lowest_eigenvalue) > 500.0:
                status.set_status_failed()
                total_energy_new.copy_to(final_energy)
                out_stream.write(f"\n FIRE2_ART: iteration {it} something is going horribly wrong {lowest_eigenvalue}\n")
                return
            #centreOfMass(velocity, natoms)

            fNorm, fMax = self.vector_norm_max(force_new)
            delta_e = abs(total_energy_old.get_total_energy() - total_energy_new.get_total_energy())

            #print(f"\n FIRE2_ART: iteration {it} {total_energy_new.get_total_energy()} lowest eigenvalue {lowest_eigenvalue} time_step {time_step}\n")
            #print(f" FIRE2_ART: gnorm {fNorm} force max {fMax} energy difference {delta_e}\n")
            if job.search_parameters.debug and it % job.search_parameters.print == 0:
                out_stream.write(f"\n FIRE2_ART: iteration {it} {total_energy_new.get_total_energy()} lowest eigenvalue {lowest_eigenvalue} time_step {time_step}\n")
                out_stream.write(f" FIRE2_ART: gnorm {fNorm} force max {fMax} energy difference {delta_e}\n")
                out_stream.flush()

            if fNorm < job.search_parameters.norm_tol or fMax < job.search_parameters.force_tol:
                if job.search_parameters.debug:
                    out_stream.write(f"\n FIRE2_ART converged : delta E {delta_e}\n")
                    out_stream.write(f" FIRE2_ART converged : max force {fMax} tolerance {job.search_parameters.force_tol}\n")
                    out_stream.write(f" FIRE2_ART converged : gnorm {fNorm} tolerance {job.search_parameters.norm_tol}\n")
                    out_stream.flush()
                total_energy_new.copy_to(final_energy)
                if lowest_eigenvalue <= job.search_parameters.max_eigenvalue:
                    status.set_status_success()
                else:
                    status.set_status_failed()
                    
                return

            total_energy_new.copy_to(total_energy_old)

        total_energy_new.copy_to(final_energy)
        status.set_status_failed()

    def vector_norm_max(self, vec):
        norm_val = np.linalg.norm(vec)
        max_val = np.max(np.abs(vec))
        return norm_val, max_val



    def calculate_lanczos_force(self, job: JobControl, bas: Config, fld: Field, total_energy: Energy, search_peprpindicular: bool, out_stream):
        
        alpha = 0.0
        beta = 0.0
        delta = 0.0
        eig_val_new = 0.0
        eig_val_old = 0.0

        converged = False

        lowest_eigenvalue = 1e6
        
        qx = np.zeros(self.vector_size, dtype=np.float64)
        q_vector = np.zeros(self.vector_size, dtype=np.float64)
        u_vector = np.zeros(self.vector_size, dtype=np.float64)
        #init_force = np.zeros(self.vector_size)
        #delta_force = np.zeros(self.vector_size)
        eig_vec_new = np.zeros(job.search_parameters.num_lanczos_vectors, dtype=np.float64)
        Q = np.zeros((job.search_parameters.num_lanczos_vectors, self.vector_size), dtype=np.float64)
        diagonal_t = np.zeros(job.search_parameters.num_lanczos_vectors, dtype=np.float64)
        sub_diag_t = np.zeros(job.search_parameters.num_lanczos_vectors - 1, dtype=np.float64)

        beta = np.linalg.norm(self.r_vector)

        init_force = self.calculate_forces(bas, fld, total_energy, out_stream)

        for i in range(job.search_parameters.num_lanczos_vectors): 
           
            for j in range(self.vector_size):
                q_vector[j] = self.r_vector[j] / beta #eqn 2
                Q[i, j] = q_vector[j]

            #do the finite difference part - eqns 3, 9 and 10
            #first create the vector q.r*
            for j in range(self.vector_size):
                qx[j] = q_vector[j] * job.search_parameters.delta_x
            bas.displace_atoms_mapped(qx)
            
            delta_force = self.calculate_forces( bas, fld, total_energy, out_stream)

            #put atoms back
            qx *= -1.0
            bas.displace_atoms_mapped(qx)

            for j in range(self.vector_size):
                u_vector[j] = -(delta_force[j] - init_force[j]) / job.search_parameters.delta_x

            if i == 0:  #eqn 4
                for j in range(self.vector_size):
                    self.r_vector[j] = u_vector[j]
            else:
                for j in range(self.vector_size):
                    self.r_vector[j]= u_vector[j] - beta * Q[i - 1, j]

            alpha = np.dot(q_vector, self.r_vector)  #eqn 5
            for j in range(self.vector_size):
                self.r_vector[j] = self.r_vector[j] - alpha * Q[i,j]  # eqn 6

            #diagonal part of matrix T
            diagonal_t[i] = alpha

            #sub diagonal parts of matrix T
            #ah blast must be greater than 0!
            if i > 0:
                sub_diag_t[i - 1] = beta

            beta = np.linalg.norm(self.r_vector)

            if i == 0:
                eig_val_new = eig_val_old = alpha
            else:
                eig_val_new, eig_vec_new = self.diagonalise(diagonal_t[:i + 1], sub_diag_t[:i])
                #get the lowest eigen value tolerance - eqn 8
                delta = abs((eig_val_new - eig_val_old) / eig_val_old)
                eig_val_old = eig_val_new
                lowest_eigenvalue = eig_val_new
                if delta < job.search_parameters.eig_val_tol:
                    converged = True

            if converged:
                break
            
        iterations = i
        eigen_vector = np.zeros(self.vector_size, dtype=np.float64)
        fvv = np.zeros(self.vector_size, dtype=np.float64)
        
        
        #v = Qv^t i.e. calculate the eigenvector of full Hessian the bit on page 9778 - bottom left para
        #for i in range(job.search_parameters.num_lanczos_vectors):
        for i in range(iterations):
            for j in range(self.vector_size):
                eigen_vector[j] += eig_vec_new[i] * Q[i, j]

        norm = np.linalg.norm(eigen_vector)
        eigen_vector /= norm

        #calculate forces using equation 18
        for i in range(self.vector_size):
            dot = np.dot(init_force, eigen_vector)
            for j in range(self.vector_size):
                fvv[j] = -dot * eigen_vector[j]

        forces = np.zeros(self.vector_size)
        if lowest_eigenvalue > job.search_parameters.max_eigenvalue or search_peprpindicular:
            forces[:] = fvv[:]  #force pependicular
        else:
            forces[:] = init_force[:] + 2.0 * fvv[:]  # force paralel
            

        np.copyto(self.r_vector[:], eigen_vector)

        # return forces on each atom even if it is zero/frozen
        atom_forces = np.zeros((bas.natoms,3), dtype=np.float64)
        j = 0
        for i in range(bas.natoms):
            if bas.frozen[i] == 0:
                atom_forces[i,0] = forces[j]
                atom_forces[i,1] = forces[j+1]
                atom_forces[i,2] = forces[j+2]
                j += 3
        #for i in range(bas.natoms):
        #    print("atom forces ", i, bas.symbol[i], atom_forces[i,:])
        #exit(-1)
        return atom_forces, lowest_eigenvalue, converged

    def calculate_forces(self, bas: Config, fld: Field, total_energy: Energy, out_stream):
        force_vector = np.zeros(self.vector_size, dtype=np.float64)
        force = np.zeros((bas.natoms, 3), dtype=np.float64)
        stress = np.zeros(6)
        rcp_vector = np.linalg.inv(bas.vectors)
        fld.calculate_forces(bas.pos, force, stress, bas.vectors, rcp_vector, bas.charge, bas.label, bas.symbol, bas.frozen, bas.natoms, total_energy)

        j = 0
        for i in range(bas.natoms):
            if bas.frozen[i] == 0:
                force_vector[j] = force[i][0]
                force_vector[j + 1] = force[i][1]
                force_vector[j + 2] = force[i][2]
                j += 3

        return force_vector

    def check_vector_size(self):
        if self.vector_size < num_lanczos_vectors:
            num_lanczos_vectors = self.vector_size

    def diagonalise(self, diagonal_t, sub_diag_t):
        eigvals, eigvecs = eigh_tridiagonal(diagonal_t, sub_diag_t)
        eig_val = eigvals[0]
        eig_vec = eigvecs[:, 0]
        return eig_val, eig_vec
    
    """def diagonalise(self, diagonal_t, sub_diag_t):
        eigvals, eigvecs = eigh_tridiagonal(diagonal_t, sub_diag_t)
        eig_val = eigvals[0]
        eig_vec = eigvecs[:, 0]
        return eig_val, eig_vec"""
    
    def check_for_message(self):
        comm = MPI.COMM_WORLD
        status = MPI.Status()
        source = 0
        command = np.zeros(2, dtype = np.int32)

        flag = comm.Iprobe(source, tag=MPI.ANY_TAG, status=status)
        if flag:
            command = comm.recv(source=source, tag=MPI.ANY_TAG, status=status)

        if status.Get_error() != MPI.SUCCESS:
            print(f"\n ART::check_for message error {status.Get_error()} {status.Get_error_class()} {self.grp_rank}\n")
            MPI.COMM_WORLD.Abort()
        return command[0], status.Get_source(), status.Get_tag()

