import numpy as np
import os

def splitm(line):
    '''
    Correctly split pdb file
    '''
    return([line[0:6].strip(),line[6:11].strip(),line[12:16].strip(),line[16:17].strip(),line[17:20].strip(),
            line[21:22].strip(),line[22:26].strip(),line[26:27].strip(),line[30:38].strip(),line[38:46].strip(),
            line[46:54].strip(),line[54:60].strip(),line[60:66].strip(),line[76:78].strip(),line[78:80].strip()])


class Geometry:

    def __init__(self, path):
        
        
        name_ext = os.path.basename(path).split(".")
        self.name = name_ext[0]
        self.ext = name_ext[1]
        self.file = path
        self.atoms = []
        self.positions = []
        self.extra = []
        self.header = []
        self.footer = []

    def read(self):
        if self.ext == "in":
            with open(self.file, "r") as f:
                for line in f:
                    data = line.split()
                    if data[0] == "atom":
                        self.atoms.append(data[4])
                        self.positions.append(np.array(data[1:4], dtype = float))
        elif self.ext == "xyz":
            with open(self.file, "r") as f:
                for line in f:
                    data = line.split()
                    self.atoms.append(data[0])
                    self.positions.append(np.array(data[1:5], dtype = float))
        elif self.ext == 'pdb':
            first_atom = 1e10
            with open(self.file, "r") as f:
                for i, line in enumerate(f):
                    data = splitm(line)
                    if data[0] == 'ATOM':
                        self.atoms.append(data[2])
                        self.positions.append(np.array(data[8:11], dtype = float))
                        data[1] = int(data[1])
                        data[6] = int(data[6])
                        data[11] = float(data[11])
                        data[12] = float(data[12])
                        self.extra.append(data)
                        first_atom = i
                    elif i <= first_atom:
                        self.header.append(line)
                    else:
                        self.footer.append(line)


    def update_pos(self, new_pos):
        self.positions = new_pos

    def write(self, name_out, ext):
        if ext == "in":
            fmt = "atom{:12.6f}{:12.6f}{:12.6f}{:>3s}\n"
            with open(f"{name_out}.{ext}", "w") as f:
                for atom, pos in zip(self.atoms, self.positions):
                    f.write(fmt.format(*pos, atom))
        elif ext == 'xyz':
            fmt = "{:>3s}{:12.6f}{:12.6f}{:12.6f}\n"
            with open(f"{name_out}.{ext}", "w") as f:
                for atom, pos in zip(self.atoms, self.positions):
                    f.write(fmt.format(atom, *pos))
        elif ext == 'pdb':
            fmt = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
            with open(f"{name_out}.{ext}", "w") as f:
                for line in self.header:
                    f.write(line)
                for i in range(len(self.positions)):
                    f.write(fmt.format(*self.extra[i][0:2], self.atoms[i], *self.extra[i][3:8], *self.positions[i], *self.extra[i][11:]))
                for line in self.footer:
                    f.write(line)