from quaternions import Quaternion # Custom library to perform rotations
from read_geo import Geometry   # Custom library to read pdb, xyz and aims
import numpy as np
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--pdb1', '-c1', required=True, type=str, help='pdb with ligand')
parser.add_argument('--pdb2', '-c2', required=True, type=str, help='pdb of tpr')
parser.add_argument('--tpr', '-s', required=True, type=str, help='tpr file')
parser.add_argument('--xtc', '-f', required=True, type=str, help='xtc file')
args = parser.parse_args()

################################################################################
#                               Useful functions                               #
################################################################################

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

def get_bb(A, names):
    '''
    Returns the indeces corresponding to atom names matching the backbone
    '''
    result = []
    for i in A:
        if names[i] in ['C', 'CA', 'N']:
            result.append(i)
    return np.array(result, dtype = int)

def angle(a, b):
    '''
    Angle between two vectors
    '''
    dot = np.dot(a, b)
    mag_a = np.linalg.norm(a)
    mag_b = np.linalg.norm(b)
    rad = np.arccos(dot/(mag_a*mag_b))
    return np.rad2deg(rad)

################################################################################
#  Read pdb with ligand figure our rotation for all pdbs based on placing the  #
#             ligand in the z axis and chain A CoM in the origin               #
################################################################################

LAST_A = 6557   # Last atom of chain A counting from 0 and non-inclusive
LIGAND = -39    # Ligand is last 39 atoms
BACKBONE = ['C', 'CA', 'N'] # Atoms to filter only the backbone
HASHES = 80*"#"
BETA_RES = list(range(466,472)) + list(range(476, 482))
print(HASHES)
print("Beta residues: ", BETA_RES)

geom = Geometry(args.pdb1)
geom.read()
N = len(geom.positions) 

# Gets only chain A and computes center of mass of bb
chain_a_pos = geom.positions[:LAST_A]
chain_a_atoms = geom.atoms[:LAST_A]
com_chain_a = com(chain_a_atoms, chain_a_pos, BACKBONE)
print("Center of mass chain A with ligand: ", com_chain_a)

# Transaltes so CoM is in the origin
new_pos = geom.positions - com_chain_a

# Gets only the ligand and computes com using C, N and O atoms
lig_pos = new_pos[LIGAND:]
lig_atoms = geom.atoms[LIGAND:]
com_lig = com(lig_atoms, lig_pos, ['C', 'N', 'O'])
print("Center of mass of ligand: ", com_lig)

# Converts CoM ligand to quaternion (x,y,z) -> (0, x i + y j + z k)
p = Quaternion.vec2quat(com_lig)

# Computes the necessary angles to rotate com of ligand unto x axis
theta1, theta2 = p.liex()

# Creates the necessary quaternions to perform the rotation along x then z
q1 = Quaternion.angaxis(theta1, [1, 0, 0])
q2 = Quaternion.angaxis(theta2, [0, 0, 1])
# Extra rotation to lie on positive z instead of x
q3 = Quaternion.angaxis(-np.pi/2, [0, 1, 0])

# For each position perform the 3 rotations and save to rotated
rotated = np.zeros((N,3))
for i in range(N):
    p = Quaternion.vec2quat(new_pos[i])
    rotated[i] = p.rotate(q3*q2*q1).tovec()

# Update positions and output the geometry to a file ligand_rotated.pdb
geom.update_pos(rotated)
geom.write('ligand_rotated', 'pdb')
print(HASHES)

################################################################################
#       Make index based on different criteria (top, bottom, betas) & bb       #
################################################################################

# Get z component of chain A
z_comp = rotated[:LIGAND,2]

# define top and bottom split by XY plane
top_ind = np.where(z_comp >= 0)[0]
bot_ind = np.where(z_comp < 0)[0]

# Get only the indeces corresponding to the bb for top, bottom and entire chain A
top_bb = get_bb(top_ind, chain_a_atoms)
bot_bb = get_bb(bot_ind, chain_a_atoms)
bb = get_bb(np.arange(0, LAST_A), chain_a_atoms)

################################################################################
#       Use gromacs to compute the CoM of the trajectory for the selections    #
################################################################################

# Read new geometry with both chains
geom = Geometry(args.pdb2)
geom.read()
N = len(geom.positions)

# Gets only chain A and computes center of mass of bb
chain_a_pos = geom.positions[:LAST_A]
chain_a_atoms = geom.atoms[:LAST_A]
com_chain_a = com(chain_a_atoms, chain_a_pos, BACKBONE)

# Translate based on com
new_pos = geom.positions - com_chain_a

# Perform rotation based on previosly defined quaternions with ligand pdb
rotated = np.zeros((N,3))
for i in range(N):
    p = Quaternion.vec2quat(new_pos[i])
    rotated[i] = p.rotate(q3*q2*q1).tovec()

# Update positions and output the geometry to a file tpr_rotate.pdb
geom.update_pos(rotated)
geom.write('tpr_rotated', 'pdb')

# Gets the indeces of the beta sheets only backbone
beta_A = []
beta_B = []
for i, ele in enumerate(geom.extra):
    if int(ele[1]) <= LAST_A:
        if ele[6] in BETA_RES:
            if ele[2] in BACKBONE:
                beta_A.append(i)
    else:
        if ele[6] in BETA_RES:
            if ele[2] in BACKBONE:
                beta_B.append(i) 

labels = ['[A]', '[B]', '[TOP_A]', '[TOP_B]', '[BOT_A]', '[BOT_B]', '[BETA_A]', '[BETA_B]']
selections = [bb + 1, bb + 1 + LAST_A, top_bb + 1, top_bb + 1 + LAST_A, 
              bot_bb + 1, bot_bb + 1 + LAST_A, beta_A, beta_B]

# Run gromacs CoM for all selections on the trajectory file
TRAJ = args.xtc
TPR = args.tpr
print("Running Gromacs")
for i in range(len(labels)):
    np.savetxt('index.ndx', selections[i], fmt='%6d', header = labels[i], comments='')
    os.system(f"echo 0 | gmx traj -f {TRAJ} -s {TPR} -n index.ndx -com -ox {labels[i][1:-1]} >> /dev/null 2>/dev/null")
    os.system("rm index.ndx")

print("# of atoms top: ", len(top_bb))
print("# of atoms bottom: ", len(bot_bb))
print(HASHES)

################################################################################
#           Calculate the collective variables CoM distance Tilt and Roll      #
################################################################################

time = np.loadtxt('A.xvg', comments = ['@', '#'])[:,0]

CoM_A = np.loadtxt('A.xvg', comments = ['@', '#'])[:,1:]
CoM_B = np.loadtxt('B.xvg', comments = ['@', '#'])[:,1:]

CoM_top_A = np.loadtxt('TOP_A.xvg', comments = ['@', '#'])[:,1:]
CoM_bot_A = np.loadtxt('BOT_A.xvg', comments = ['@', '#'])[:,1:]

CoM_top_B = np.loadtxt('TOP_B.xvg', comments = ['@', '#'])[:,1:]
CoM_bot_B = np.loadtxt('BOT_B.xvg', comments = ['@', '#'])[:,1:]

CoM_beta_A = np.loadtxt('BETA_A.xvg', comments = ['@', '#'])[:,1:]
CoM_beta_B = np.loadtxt('BETA_B.xvg', comments = ['@', '#'])[:,1:]

x_axis_A = CoM_top_A - CoM_bot_A
x_axis_B = CoM_top_B - CoM_bot_B

y_axis_A = CoM_beta_A - CoM_A
y_axis_B = CoM_beta_B - CoM_B

tilt = []
roll = []
for i in range(len(time)):
    tilt.append(angle(x_axis_A[i], x_axis_B[i]))
    roll.append(angle(y_axis_A[i], y_axis_B[i]))
dist = np.linalg.norm(CoM_A - CoM_B, axis=1)

np.savetxt('colvar.dat', list(zip(time, dist, tilt, roll)), fmt="%12.1f, %6.3f, %6.3f, %6.3f")
print(np.round(np.mean(dist),3), np.round(np.std(dist),3))
print(np.round(np.mean(tilt),3), np.round(np.std(tilt),3))
print(np.round(np.mean(roll),3), np.round(np.std(roll),3))
print("Finished collective varibles")
print(HASHES)

################################################################################
#     Read pdb with both chains and compute the CoM of all the selections      #
#           Create pymol script to visualize the newly defined axis            #
################################################################################

coms = []
COLOR_X = "1.0, 0.2, 0.2"
COLOR_Y = "0.5, 0.5, 1.0"
COLOR_A = "white"
COLOR_B = "bluewhite"
COLOR_x = "tv_red"
COLOR_y = "slate"
radius = 0.4
vdw = 1.5

com_A = com(geom.atoms[:LAST_A], geom.positions[:LAST_A], BACKBONE)
coms.append(com_A)

com_B = com(geom.atoms[LAST_A:], geom.positions[LAST_A:], BACKBONE)
coms.append(com_B)

pos = [geom.positions[i] for i in top_bb] 
atoms = [geom.atoms[i] for i in top_bb] 
com_top_A = com(atoms, pos, BACKBONE)
coms.append(com_top_A)

pos = [geom.positions[i] for i in bot_bb] 
atoms = [geom.atoms[i] for i in bot_bb] 
com_bot_A = com(atoms, pos, BACKBONE)
coms.append(com_bot_A)

pos = [geom.positions[i] for i in top_bb + LAST_A] 
atoms = [geom.atoms[i] for i in top_bb + LAST_A] 
com_top_B = com(atoms, pos, BACKBONE)
coms.append(com_top_B)

pos = [geom.positions[i] for i in bot_bb + LAST_A] 
atoms = [geom.atoms[i] for i in bot_bb + LAST_A] 
com_bot_B = com(atoms, pos, BACKBONE)
coms.append(com_bot_B)

pos = [geom.positions[i] for i in beta_A] 
atoms = [geom.atoms[i] for i in beta_A] 
com_beta_A = com(atoms, pos, BACKBONE)
coms.append(com_beta_A)

pos = [geom.positions[i] for i in beta_B] 
atoms = [geom.atoms[i] for i in beta_B] 
com_beta_B = com(atoms, pos, BACKBONE)
coms.append(com_beta_B)

pairs = [(0, 6), (1, 7), (2, 3), (4, 5)]

with open("vis.py", 'w') as f:
    f.write("from pymol.cgo import *\n")
    f.write("from pymol import cmd\n\n")
    f.write("cmd.load(\"tpr_rotated.pdb\")\n")

    for i in range(len(coms)):
        if i in [2,3,4,5]:
            f.write(f"cmd.pseudoatom(\"{labels[i][1:-1]}\", pos={list(coms[i])},\
                      chain = \"X\", vdw = {vdw})\n")
        else:
            f.write(f"cmd.pseudoatom(\"{labels[i][1:-1]}\", pos={list(coms[i])},\
                      chain = \"Z\", vdw = {vdw})\n")
    
    for i, pair in enumerate(pairs):
        p1 = coms[pair[0]]
        p2 = coms[pair[1]]
        if i in [0,1]:
            f.write(f"cmd.load_cgo([9.0,{p1[0]},{p1[1]},{p1[2]},{p2[0]},{p2[1]},{p2[2]},{radius},{COLOR_Y},{COLOR_Y}],\"cyl_{i}\")\n") 
        else:
            f.write(f"cmd.load_cgo([9.0,{p1[0]},{p1[1]},{p1[2]},{p2[0]},{p2[1]},{p2[2]},{radius},{COLOR_X},{COLOR_X}],\"cyl_{i}\")\n") 

    f.write("cmd.show(\"spheres\", selection = \"chain X\")\n")
    f.write("cmd.show(\"spheres\", selection = \"chain Z\")\n")
    f.write(f"cmd.color(\"{COLOR_x}\", selection = \"chain X\")\n")
    f.write(f"cmd.color(\"{COLOR_y}\", selection = \"chain Z\")\n")

    f.write(f"cmd.color(\"{COLOR_A}\", selection = \"chain A\")\n")
    f.write(f"cmd.color(\"{COLOR_B}\", selection = \"chain B\")\n")

    f.write(f"cmd.bg_color(\"white\")\n")
    f.write(f"cmd.set(\"ray_trace_mode\", 1)\n")

    f.write("cmd.set_view([\
    -1.0, 0.0, 0.0,\
     0.0, 0.0, 1.0,\
     0.0, 1.0, 0.0,\
     0.0, 0.0, -300.0,\
     22.0, 22.0, -22.0,\
    -100000.0, 10000.0, 20.0])\n")

print("Ceated PyMOL script")
print(HASHES)