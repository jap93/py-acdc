import numpy as np

from constants import BOLTZMANN

class Event (object):

    def __init__(self, basin1Energy = 1.0e8, basin2Energy = 1.0e8, saddleEnergy = 1.0e8, activationEnergy = 1.0e8, centre = None, num = -1):

        self.basin1Energy = basin1Energy
        self.basin2Energy = basin2Energy
        self.saddleEnergy = saddleEnergy
        self.activationEnergy = activationEnergy
        self.image = num
        self.rate = 1.0e6

        self.centre = np.zeros(3)
        np.copyto(self.centre, centre)

    def printEvent(self):

        print(" basin 1 energy {:15.8f}".format(self.basin1Energy))
        print(" basin 2 energy {:15.8f}".format(self.basin2Energy))
        print(" transition sate energy {:15.8f}".format(self.saddleEnergy))
        print(" activation energy {:15.8f}".format(self.activationEnergy))
        print(" rate              {:15.8f}".format(self.rate))
        print(" centre            {:15.8f} {:15.8f} {:15.8f}".format(self.centre[0], self.centre[1], self.centre[2]))
        print(" image number      {:8d}".format(self.image))

    def write_event(self, fp):

        fp.write("EVENT: \n")
        fp.write("    basin_1_energy {:15.8f} \n".format(self.basin1Energy))
        fp.write("    basin_2_energy {:15.8f} \n".format(self.basin2Energy))
        fp.write("    transition_sate_energy {:15.8f} \n".format(self.saddleEnergy))
        fp.write("    activation_energy {:15.8f} \n".format(self.activationEnergy))
        fp.write("    rate              {:15.8f} \n".format(self.rate))
        fp.write("    centre            {:15.8f} {:15.8f} {:15.8f} \n".format(self.centre[0], self.centre[1], self.centre[2]))
        fp.write("    image_number      {:8d} \n".format(self.image))

    def get_activation_energy(self):
        return self.activationEnergy
    
    def set_activation_energy(self, ea):
        self.activationEnergy = ea

    def getBasin1Energy(self):
        return self.basin1Energy

    def setBasin1Energy(self, eng):
        self.basin1Energy = eng

    def getBasin2Energy(self):
        return self.basin2Energy

    def getSaddleEnergy(self):
        return self.saddleEnergy

    def calculate_rate(self, temperature, prefactor):

        beta = 1.0 / (temperature * BOLTZMANN)
        #print("calculate rate ", temperature, prefactor, boltzmann, beta)
        #print("rate again", (-self.activationEnergy * beta), np.exp(-self.activationEnergy*beta))
        self.rate = prefactor * np.exp(-self.activationEnergy * beta)
        return self.rate

    def get_rate(self):
        return self.rate

    def get_search_number(self):
        return self.image
    
    def set_search_number(self, num):
        self.image = num
    
    def get_centre(self):
        return self.centre
