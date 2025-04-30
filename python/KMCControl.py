from constants import BOLTZMANN
def split(s):
    return s.split()

class KMCCntrl:
    def __init__(self):
        self.kmc_steps = 0
        self.kmc_min_cap = 0.001
        self.kmcWindow = 50.0
        self.kmcExcessSearches = 0
        self.kmcPreFactor = 0.0
        self.kmcTemperature = 300.0
        self.kmc_max_cap = 10.0
        self.kmc_images = 0
        self.basin2Delta = 0.02
        self.fullBasin2 = False
        self.kmc_method = ""
        self.kmcBasinRadius = 1.0
        self.regionStyle = 0
        self.regionRadius = 100.0
        self.artTypes1 = []
        self.artTypes2 = []
        self.coordDistance = []
        self.idealCoordNumber = []
        self.useGauss = False
        self.gaussWidth = 0.0
        self.min_defect_distance = 0.5
        self.recycle_events = False
        self.recycle_style = 0
        self.recycle_distance = 5.0
        self.recycle_prob = 1.1
        self.recycle_limit = 100.0

    def writeKMCDetails(self, outStream):

        outStream.write("\n\n *****************************************************************************************************")
        outStream.write("\n Kinetic Monte Carlo Details")
        outStream.write("\n ***************************************************************************************************** \n")

        
        numErrors = 0

        outStream.write(f" the number of KMC steps {self.kmc_steps}\n")
        outStream.write(f" the min activation energy for barriers {self.kmc_min_cap}\n")
        outStream.write(f" the energy window for the activation energy {self.kmcWindow}\n")
        outStream.write(f" the number of excess searches to be performed {self.kmcExcessSearches}\n")

        if self.kmcPreFactor < 0.0:
            outStream.write(" the pre-factor is less than zero and probably has not been set \n")
            numErrors += 1
        else:
            outStream.write(f" the KMC pre-factor {self.kmcPreFactor}\n")

        self.kmc_max_cap = self.kmc_min_cap + (self.kmcWindow * self.kmcTemperature * BOLTZMANN)
        outStream.write(f"\n maximum activation energy : {self.kmc_max_cap}\n")
        
        outStream.write(f" the temperature of the simulation {self.kmcTemperature} K \n")
        outStream.write(f" the number of KMC events per step {self.kmc_images}\n")
        outStream.write(f" distance of displacement from saddle prior to relaxation to basin 2 {self.basin2Delta} A\n")
        
        if self.fullBasin2:
            outStream.write(" the complete basin 2 will be relaxed at the end of the cycle \n")
        else:
            outStream.write(" the partial basin 2 will be relaxed at the end of the cycle \n")

        if self.kmc_method == "art":
            outStream.write(" the activation relaxation technique will be used to search for transition states\n")
        elif self.kmc_method == "dimer":
            outStream.write(" the dimer method will be used to search for transition states\n")
        else:
            outStream.write(" unrecognised activation method \n")
            numErrors += 1
        
        outStream.write(f" distance for new basin {self.kmcBasinRadius} A\n")
        
        outStream.write(" target atom details \n")
        if self.regionStyle == 1:
            outStream.write(" all atoms will be displaced\n")
        elif self.regionStyle == 2:
            outStream.write(f" a atom will be selected at random and those within {self.regionRadius} A will be relaxed \n")
        elif self.regionStyle == 3:
            outStream.write(" all atoms will be displaced but centered on types\n")
            for art_type in self.artTypes1:
                outStream.write(f" type {art_type}\n")
        elif self.regionStyle == 4:
            outStream.write(f" a atom will be selected at random and those within {self.regionRadius} A will be relaxed \n")
            for art_type in self.artTypes1:
                outStream.write(f" type {art_type}\n")
        elif self.regionStyle == 5:
            outStream.write(f" a atom will be selected at based on coordination and those within {self.regionRadius} A will be relaxed \n")
            for i in range(len(self.artTypes1)):
                outStream.write(f" type 1 {self.artTypes1[i]} and {self.artTypes2[i]}. Max distance {self.coordDistance[i]} ideal coordination {self.idealCoordNumber[i]}\n")
        elif self.regionStyle == 6:
             outStream.write(f" a atom will be selected at based on template and radius {self.min_defect_distance}. Atoms within {self.regionRadius} A will be relaxed \n")
        else:
            outStream.write(" unrecognised target style\n")
            numErrors += 1
            
        if self.useGauss:
            outStream.write(f" the displacement of atoms will be weighted by a Gaussian of width {self.gaussWidth}\n")

        if self.kmc_method != "art" and self.kmc_method != "dimer":
            outStream.write(" the kmc method is not recognised\n")
            numErrors += 1

        return numErrors

    def readKineticMonteCarlo(self, inStream, outStream):
        readMore = True

        while readMore:
            line = inStream.readline().strip().lower()
            if not line:
                break

            words = split(line)
            if not words:
                continue

            keyWord = words[0]

            if keyWord == "}":
                readMore = False
                break
            elif keyWord == "kmcsteps":
                self.kmc_steps = int(words[1])
            elif keyWord == "mincap":
                self.kmc_min_cap = float(words[1])
            elif keyWord == "prefactor":
                self.kmcPreFactor = float(words[1])
            elif keyWord == "window":
                self.kmcWindow = float(words[1])
            elif keyWord == "kmctemperature":
                self.kmcTemperature = float(words[1])
            elif keyWord == "kmcevents":
                self.kmc_images = int(words[1])
            elif keyWord == "excesssearches":
                self.kmcExcessSearches = int(words[1])
            elif keyWord == "basin2delta":
                self.basin2Delta = float(words[1])
            elif keyWord == "partialbasin2":
                self.fullBasin2 = False
            elif keyWord == "kmcmethod":
                self.kmc_method = words[1]
            elif keyWord == "kmcbasinradius":
                self.kmcBasinRadius = float(words[1])
            elif keyWord == "target":
                subWord = words[1]
                if subWord == "global":
                    self.regionStyle = 1
                elif subWord == "local":
                    self.regionStyle = 2
                    self.regionRadius = float(words[2])
                elif subWord == "globaltypes":
                    self.regionStyle = 3
                    num = int(words[2])
                    for _ in range(num):
                        line = inStream.readline().strip()
                        self.artTypes1.append(line.split()[0])
                elif subWord == "localtypes":
                    self.regionStyle = 4
                    num = int(words[2])
                    self.regionRadius = float(words[3])
                    for _ in range(num):
                        line = inStream.readline().strip()
                        self.artTypes1.append(line.split()[0])
                elif subWord == "coordination":
                    self.regionStyle = 5
                    num = int(words[2])
                    self.regionRadius = float(words[3])
                    for _ in range(num):
                        line = inStream.readline() #.strip()
                        parts = line.split()
                        self.artTypes1.append(parts[0])
                        self.artTypes2.append(parts[1])
                        self.coordDistance.append(float(parts[2]))
                        self.idealCoordNumber.append(int(parts[3]))
                elif subWord == "template":
                    self.regionStyle = 6
                    self.min_defect_distance = float(words[2])
                    self.regionRadius = float(words[3])
            elif keyWord == "usegaussian":
                self.useGauss = True
                self.gaussWidth = float(words[1])
            elif keyWord == "recycle":
                self.recycle_events = True

                if words[1] == "distance":
                    self.recycle_style = 0
                    self.recycle_distance = float(words[2])
                    self.recycle_prob = 1.1
                elif words[1] == "distprob":
                    self.recycle_style = 1
                    self.recycle_distance = float(words[2])
                    self.recycle_prob = float(words[3])
                elif words[1] == "limit":
                    self.recycle_style = 2
                    self.recycle_distance = float(words[2])
                    self.recycle_limit = float(words[3])
            else:
                print(f"\n KMC keyword not found {words[0]}")
                outStream.write(f"\n KMC keyword not found {words[0]}")
                outStream.flush()
                #MPI.COMM_WORLD.Abort()
                exit(1)