import numpy as np

class Energy:
    def __init__(self, src=None):
        if src:
            self.totalEnergy = src.totalEnergy
            self.rcpEnergy = src.rcpEnergy
            self.realEnergy = src.realEnergy
            self.vdwEnergy = src.vdwEnergy
            self.manyEnergy = src.manyEnergy
            self.selfEnergy = src.selfEnergy
            self.extEnergy = src.extEnergy
        else:
            self.totalEnergy: np.float64 = 0.0
            self.rcpEnergy: np.float64 = 0.0
            self.realEnergy: np.float64 = 0.0
            self.vdwEnergy: np.float64 = 0.0
            self.manyEnergy: np.float64 = 0.0
            self.selfEnergy: np.float64 = 0.0
            self.extEnergy: np.float64 = 0.0

    def __sub__(self, src):
        result = Energy()
        result.totalEnergy = self.totalEnergy - src.totalEnergy
        result.rcpEnergy = self.rcpEnergy - src.rcpEnergy
        result.realEnergy = self.realEnergy - src.realEnergy
        result.selfEnergy = self.selfEnergy - src.selfEnergy
        result.manyEnergy = self.manyEnergy - src.manyEnergy
        result.vdwEnergy = self.vdwEnergy - src.vdwEnergy
        result.extEnergy = self.extEnergy - src.extEnergy
        return result

    def __add__(self, src):
        result = Energy()
        result.totalEnergy = self.totalEnergy + src.totalEnergy
        result.rcpEnergy = self.rcpEnergy + src.rcpEnergy
        result.realEnergy = self.realEnergy + src.realEnergy
        result.selfEnergy = self.selfEnergy + src.selfEnergy
        result.manyEnergy = self.manyEnergy + src.manyEnergy
        result.vdwEnergy = self.vdwEnergy + src.vdwEnergy
        result.extEnergy = self.extEnergy + src.extEnergy
        return result

    def __iadd__(self, src):
        self.totalEnergy += src.totalEnergy
        self.manyEnergy += src.manyEnergy
        self.realEnergy += src.realEnergy
        self.vdwEnergy += src.vdwEnergy
        self.selfEnergy += src.selfEnergy
        self.extEnergy += src.extEnergy
        return self

    def zero(self):
        self.totalEnergy = 0.0
        self.manyEnergy = 0.0
        self.rcpEnergy = 0.0
        self.selfEnergy = 0.0
        self.realEnergy = 0.0
        self.vdwEnergy = 0.0
        self.extEnergy = 0.0

    def copy_to(self, c):
        c.totalEnergy = self.totalEnergy
        c.rcpEnergy = self.rcpEnergy
        c.realEnergy = self.realEnergy
        c.selfEnergy = self.selfEnergy
        c.manyEnergy = self.manyEnergy
        c.vdwEnergy = self.vdwEnergy
        c.extEnergy = self.extEnergy

    def print_energy(self, box, outStream):
        totalEnergy = self.rcpEnergy + self.realEnergy + self.vdwEnergy + self.manyEnergy + self.extEnergy
        outStream.write(f"\n\n energies of box {box}\n")
        outStream.write(f" total (internal) energy      {totalEnergy:25.15e}\n")
        outStream.write(f" total madelung energy        {(self.rcpEnergy + self.realEnergy):25.15e}\n")
        outStream.write(f" recip space energy           {self.rcpEnergy:25.15e}\n")
        outStream.write(f" real space energy            {self.realEnergy:25.15e}\n")
        outStream.write(f" non-bonded vdw energy        {self.vdwEnergy:25.15e}\n")
        outStream.write(f" many body energy             {self.manyEnergy:25.15e}\n")
        outStream.write(f" external potential energy    {self.extEnergy:25.15e}\n")

    def get_total_energy(self) -> np.float64:
        self.totalEnergy = self.rcpEnergy + self.realEnergy + self.vdwEnergy + self.manyEnergy + self.extEnergy
        return self.totalEnergy