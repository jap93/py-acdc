import os
import numpy as np
import warnings

from mpi4py import MPI

from config import Config
from field import Field
from energy import Energy
from status import Status
from relax_control import RelaxCntrl
from constants import EVTODL


class Relax:
    def __init__(self):
        self.euler = False

    def relax_structure(self, cfg: Config, fld: Field, natoms, r_cntrl: RelaxCntrl,  out_stream):

        status = Status()
        eng = Energy()
            
        ase_relax = False
        flag = (os.environ.get("USE_ASE_RELAX", ase_relax) == "True")
        
        if flag == True:
            from ase import Atoms
            from ase.optimize import LBFGS, BFGS
            from ase.filters import UnitCellFilter

            basin1 = Atoms(positions=cfg.pos, cell=cfg.vectors, symbols=cfg.symbol, pbc=cfg.pbc)
            
            basin1.set_calculator(fld.get_calculator())
            
            flag = LBFGS(UnitCellFilter(basin1, mask=[1,1,1,0,0,0])).run(fmax=r_cntrl.force_tol, steps=r_cntrl.max_iter)  #constant volume at the moment
            
            np.copyto(cfg.pos, basin1.get_positions())
            np.copyto(cfg.vectors, basin1.cell[:])
            eng.vdwEnergy = basin1.get_potential_energy()

            if flag == True:
                status.set_status_success()
            else:
                status.set_status_failed()
                
        else:
            if r_cntrl.verlet_integration:
                self.euler = False

            if r_cntrl.method == "fire":
                self.fire_relax(cfg, fld, natoms, r_cntrl, status, eng, out_stream)
            elif r_cntrl.method == "fire2":
                self.fire2_relax(cfg, fld, natoms, r_cntrl, status, eng, out_stream)
            else:
                stress = np.zeros(6, dtype=np.float64)
                force_new = np.zeros((natoms, 3))
                rcp_vector = np.linalg.inv(cfg.vectors)
                fld.calculate_forces(cfg.pos, force_new, stress, cfg.vectors, rcp_vector, cfg.charge, cfg.label, cfg.symbol, cfg.frozen, natoms, eng)
                status.set_status_success()

        return eng, status

    def fire_relax(self, bas, fld, natoms, r_cntrl, status, eng, out_stream):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        n_step = 1
        alpha = r_cntrl.alpha_start

        time_step = r_cntrl.init_time_step
        time_step_max = r_cntrl.max_time_step

        stress = np.zeros(6, dtype=np.float64)
        velocity = np.zeros((natoms, 3))
        force_new = np.zeros((natoms, 3))
        rcp_vector = np.linalg.inv(bas.vectors)

        fld.calculate_forces(bas.pos, force_new, stress, bas.vectors, rcp_vector, bas.charge, bas.label, bas.symbol, bas.frozen, natoms, eng)

        total_energy_old = eng.get_total_energy()        

        f_norm, f_max = self.vector_norm_max(force_new[:][0], force_new[:][1], force_new[:][2])

        if rank == 0 and r_cntrl.debug:
            out_stream.write(f"\n FIRE: initial total energy\n")
            eng.print_energy(0, out_stream)
            out_stream.write(f" FIRE: gnorm {f_norm} force max {f_max}\n")
            out_stream.flush()

        half_step = 0.5 * time_step * EVTODL

        for it in range(1, r_cntrl.max_iter + 1):
            P = np.dot(force_new.flatten(), velocity.flatten())

            one_less_alpha = 1.0 - alpha

            f_norm, f_max = self.vector_norm_max(force_new[:][0], force_new[:][1], force_new[:][2])
            v_norm, v_max = self.vector_norm_max(velocity[:][0], velocity[:][1], velocity[:][2])

            with warnings.catch_warnings():
                warnings.filterwarnings('error')
    
                try:
                    alpha_force = alpha * v_norm / f_norm   # this gets caught and handled as an exception
                except Warning as e:
                    print(f"divide by zero {it } rank {rank} fnorm {f_norm} \n")
                    alpha_force = alpha
                        
            if P > 0:
                velocity = one_less_alpha * velocity + alpha_force * force_new
                if n_step > r_cntrl.N_min:
                    time_step = min(time_step * r_cntrl.step_increase, time_step_max)
                    alpha *= r_cntrl.alpha_decrease
                n_step += 1
            else:
                n_step = 0
                time_step *= r_cntrl.step_decrease
                alpha = r_cntrl.alpha_start
                velocity.fill(0)

            if self.euler:
                fact = time_step * EVTODL
                bas.pos[:, 0] += time_step * velocity[:, 0]
                bas.pos[:, 1] += time_step * velocity[:, 1]
                bas.pos[:, 2] += time_step * velocity[:, 2]
                velocity += fact * force_new
            else:
                velocity += half_step * force_new
                bas.pos[:, 0] += time_step * velocity[:, 0]
                bas.pos[:, 1] += time_step * velocity[:, 1]
                bas.pos[:, 2] += time_step * velocity[:, 2]

            fld.calculate_forces(bas.pos, force_new, stress, bas.vectors, rcp_vector, bas.charge, bas.label, bas.symbol, bas.frozen, natoms, eng)

            if not self.euler:
                velocity += half_step * force_new

            #self.centre_of_mass(velocity)

            total_energy_new = eng.get_total_energy()

            f_norm, f_max = self.vector_norm_max(force_new[:][0], force_new[:][1], force_new[:][2])
            delta_e = abs(total_energy_old - total_energy_new)

            if rank == 0 and r_cntrl.debug and it % r_cntrl.print == 0:
                out_stream.write(f"\n FIRE: iteration {it} {total_energy_new} timestep {time_step}\n")
                out_stream.write(f" FIRE: gnorm {f_norm} force max {f_max} energy difference {delta_e}\n")
                out_stream.flush()

            if f_norm < r_cntrl.norm_tol or f_max < r_cntrl.force_tol:
                if rank == 0 and r_cntrl.debug:
                    out_stream.write(f"\n FIRE converged: delta E {delta_e}\n")
                    out_stream.write(f" FIRE converged: max force {f_max} tolerance {r_cntrl.force_tol}\n")
                    out_stream.write(f" FIRE converged: gnorm {f_norm} tolerance {r_cntrl.norm_tol}\n")
                    out_stream.flush()
                status.set_status_success()
                return

            total_energy_old = total_energy_new

        status.set_status_failed()




    def fire2_relax(self, bas: Config, fld, natoms, rCntrl: RelaxCntrl, status: Status, eng: Energy, outStream):
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        nUpHill = 0
        nStep = 1

        timeStep = rCntrl.init_time_step
        timeStepMax = rCntrl.max_time_step
        timeStepMin = rCntrl.min_time_step
        alpha = rCntrl.alpha_start

        stress = np.zeros(6, dtype=np.float64)
        velocity = np.zeros((natoms, 3), dtype=np.float64)
        forceNew = np.zeros((natoms, 3), dtype=np.float64)
        rcp_vector = np.linalg.inv(bas.vectors)

        fld.calculate_forces(bas.pos, forceNew, stress, bas.vectors, rcp_vector, bas.charge, bas.label, bas.symbol, bas.frozen, natoms, eng)

        total_energy_old = eng.get_total_energy()

        fNorm, fMax = self.vector_norm_max(forceNew[:][0], forceNew[:][1], forceNew[:][2])

        if rank == 0 and rCntrl.debug:
            outStream.write(f"\n FIRE: initial total energy\n")
            eng.print_energy(0, outStream)
            outStream.write(f" FIRE: gnorm {fNorm} force max {fMax}\n")
            outStream.flush()

        halfStep = 0.5 * timeStep * EVTODL

        for it in range(1, rCntrl.max_iter + 1):
            P = self.dotProduct(forceNew[:][0], forceNew[:][1], forceNew[:][2], velocity[:][0], velocity[:][1], velocity[:][2])

            OneLessAlpha = 1.0 - alpha
            fNorm, fMax = self.vector_norm_max(forceNew[:][0], forceNew[:][1], forceNew[:][2])
            vNorm, vMax = self.vector_norm_max(velocity[:][0], velocity[:][1], velocity[:][2])

            
            if P > 0:
                alphaForce = alpha * vNorm / fNorm

                velocity = OneLessAlpha * velocity + alphaForce * forceNew

                if nStep > rCntrl.N_min:
                    timeStep = min(timeStep * rCntrl.step_increase, timeStepMax)
                    alpha *= rCntrl.alpha_decrease

                nStep += 1
                nUpHill = 0
            else:
                nStep = 0
                nUpHill += 1

                if nUpHill > rCntrl.N_maxUphill:
                    status.set_status_failed()
                    return

                if rCntrl.initialDelay:
                    if it >= rCntrl.N_delaySteps:
                        timeStep *= rCntrl.step_decrease
                        if timeStep < timeStepMin:
                            timeStep = timeStepMin
                        alpha = rCntrl.alpha_start
                else:
                    timeStep *= rCntrl.step_decrease
                    alpha = rCntrl.alpha_start

                bas.pos -= 0.5 * timeStep * velocity

                velocity.fill(0.0)

            if self.euler:
                fact = timeStep * EVTODL
                for i in range(natoms):
                    if bas.frozen[i] == 0:
                        if P > 0.0:
                            velocity[i][:] = OneLessAlpha * velocity[i][:] + alphaForce * forceNew[i][:]
                        
                        bas.pos[i][:] += timeStep * velocity[i][:]
                    
                        velocity[i][:] += fact * forceNew[i][:]
                    
            else:
                halfStep = 0.5 * timeStep * EVTODL
                for i in range(natoms):
                    if bas.frozen[i] == 0:
                        if P > 0.0:
                            velocity[i][:] = OneLessAlpha * velocity[i][:] + alphaForce * forceNew[i][:]
                        
                        velocity[i][:] += halfStep * forceNew[i][:]
                    
                bas.pos += timeStep * velocity
            

            fld.calculate_forces(bas.pos, forceNew, stress, bas.vectors, rcp_vector, bas.charge, bas.label, bas.symbol, bas.frozen, natoms, eng)

            if not self.euler:
                for i in range(natoms):
                    if bas.frozen[i] == 0:
                        velocity[i][:] += halfStep * forceNew[i][:]
                    

            #centreOfMass(velocity, natoms)

            total_energy_new = eng.get_total_energy()

            fNorm, fMax = self.vector_norm_max(forceNew[:][0], forceNew[:][1], forceNew[:][2])
            deltaE = abs(total_energy_old - total_energy_new)

            if rank == 0 and rCntrl.debug and it % rCntrl.print == 0:
                outStream.write(f"\n FIRE2: iteration {it} {total_energy_new} timestep {timeStep}\n")
                outStream.write(f" FIRE2: gnorm {fNorm} force max {fMax} energy difference {deltaE}\n")

            if fNorm < rCntrl.norm_tol or fMax < rCntrl.force_tol:
                if rank == 0 and rCntrl.debug:
                    outStream.write(f"\n FIRE2 converged : delta E {deltaE}\n")
                    outStream.write(f" FIRE2 converged : max force {fMax} tolerance {rCntrl.force_tol}\n")
                    outStream.write(f" FIRE2 converged : gnorm {fNorm} tolerance {rCntrl.norm_tol}\n")
                    outStream.flush()

                status.set_status_success()
                return

            total_energy_old = total_energy_new

        status.set_status_failed()

# Helper functions (placeholders, to be implemented)
    def vector_norm_max(self, x, y, z):
        norms = np.sqrt(x**2 + y**2 + z**2)
        fNorm = np.linalg.norm(norms)
        fMax = np.max(norms)
        return fNorm, fMax

    def dotProduct(self, x1, y1, z1, x2, y2, z2):
        return np.dot(x1, x2) + np.dot(y1, y2) + np.dot(z1, z2)

    def centreOfMass(self, velocityX, velocityY, velocityZ, natoms):
        # Placeholder implementation
        pass

"""
    @staticmethod
    def vector_norm_max(vector):
        norms = np.linalg.norm(vector, axis=1)
        max_norm = np.max(norms)
        norm_sum = np.sqrt(np.sum(norms ** 2))
        return norm_sum, max_norm

    @staticmethod
    def centre_of_mass(velocity):
        mass_center = np.mean(velocity, axis=0)
        velocity -= mass_center

"""
