import os

import constants as constants
class Status:
    SUCCESS = constants.SUCCESS
    FAILED = constants.FAILED
    INTERRUPTED = constants.INTERRUPTED
    TERMINATED = constants.TERMINATED

    def __init__(self):
        self.sim_status = self.FAILED

    def get_status(self):
        """Returns the status of the simulation"""
        return self.sim_status

    def set_status_failed(self):
        """Set the Status object to failed"""
        self.sim_status = self.FAILED

    def set_status_success(self):
        """Set the Status object to success"""
        self.sim_status = self.SUCCESS

    def set_status_interrupt(self):
        """Set the Status to indicate that it has been interrupted, but will continue, by MPI"""
        self.sim_status = self.INTERRUPTED

    def set_status_terminated(self):
        """Set the Status to indicate that it has been terminated by MPI"""
        self.sim_status = self.TERMINATED

    def write(self, new_path):
        file_name = os.path.join(new_path, "STATUS")
        with open(file_name, "w") as outstream:
            outstream.write(f"{self.sim_status}")