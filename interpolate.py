from quaternions import Quaternion
from read_geo import Geometry
import numpy as np
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('--pdb', '-c', required=True, type=str, help='pdb with ligand')
args = parser.parse_args()

def com(names, positions, use):
    '''
    Center of mass for a given selection using certain selected atoms
    '''
    vec = np.zeros(3)
    mass_t = 0.0
    masses = {'C': 12.011, 'CA': 12.011, 'O': 15.999, 'N': 14.007}

    for i in range(len(positions)):
        if names[i] in use:
            vec += positions[i]*masses[names[i]]
            mass_t += masses[names[i]]
    return vec/mass_t

LIGAND = -39    # Ligand is last 39 atoms
LAST_A = 6557
geom = Geometry(args.pdb)
geom.read()
N = len(geom.positions) 

# Gets only chain A and computes center of mass of bb
chain_a_pos = geom.positions[:LAST_A]
chain_a_atoms = geom.atoms[:LAST_A]
com_chain_a = com(chain_a_atoms, chain_a_pos,  ['C', 'CA', 'N'])

# Transaltes so CoM is in the origin
new_pos = geom.positions - com_chain_a

# Gets only the ligand and computes com using C, N and O atoms
lig_pos = new_pos[LIGAND:]
lig_atoms = geom.atoms[LIGAND:]
com_lig = com(lig_atoms, lig_pos, ['C', 'N', 'O'])

# Converts CoM ligand to quaternion (x,y,z) -> (0, x i + y j + z k)
initial_q = Quaternion.vec2quat(com_lig)
# Creates the necessary quaternions to perform the rotation along x then z
q = Quaternion.angaxis(np.pi, [1, 1, 1])
rotated_q = initial_q.rotate(q)

for i in np.arange(0.0,1.01,0.01):
    print(i)
    tmp = initial_q.interpolate(rotated_q, i).normalize()
    rotated = np.zeros((N,3))
    for j in range(N):
        p = Quaternion.vec2quat(new_pos[j])
        rotated[j] = p.rotate(tmp).tovec()

    geom.update_pos(rotated)
    geom.write(f'ligand_interpolate', 'pdb', "a", int((i*100) + 1))