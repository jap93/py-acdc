import re

class Element:
    def __init__(self, name="", mass=0.0, charge=0.0, atomic_number=0):
        self.name = name
        self.mass = mass
        self.charge = charge
        self.atomic_number = atomic_number

class Species:
    def __init__(self):
        self.number_of_elements = 0
        self.ele_data = []

    def __del__(self):
        pass

    def load_species(self, in_stream, num):
        self.number_of_elements = num

        for _ in range(self.number_of_elements):
            line = in_stream.readline().strip()
            words = self.split(line)

            ele = Element()
            ele.name = words[0]
            ele.mass = float(words[1])
            ele.charge = float(words[2])
            if len(words) > 3:
                ele.atomic_number = int(words[3])

            self.ele_data.append(ele)

    def read_species(self, in_stream):

        for _ in range(10000):
            line = in_stream.readline().strip()
            words = self.split(line)
            if words[0].lower() == "species":
                self.number_of_elements = int(words[1])
                break

        for _ in range(self.number_of_elements):
            line = in_stream.readline().strip()
            words = self.split(line)

            ele = Element()
            ele.name = words[0]
            ele.mass = float(words[1])
            ele.charge = float(words[2])
            if len(words) > 3:
                ele.atomic_number = int(words[3])

            self.ele_data.append(ele)

    def print_species(self, out_stream):
        out_stream.write("\n \n *********************************************************************\n")
        out_stream.write("                         atom data \n")
        out_stream.write("\n \n *********************************************************************\n")
        out_stream.write("\n        element         charge           mass\n")

        for ele in self.ele_data:
            out_stream.write(f"{ele.name:15} {ele.charge:15.3f} {ele.mass:15.3f}\n")

    def get_species(self, i) -> Element:
        return self.ele_data[i]
    
    def get_num_species(self):
        return self.number_of_elements

    def get_mass(self, i):
        return self.ele_data[i].mass

    def split(self, s):
        return re.split(r'\s+', s.strip())

# Example usage:
# To load data from a file, you can use the following:
# with open('data.txt', 'r') as infile:
#     species = Species()
#     species.load_species(infile, num_elements)
# To print data to a file, you can use:
# with open('output.txt', 'w') as outfile:
#     species.print_species(outfile)

