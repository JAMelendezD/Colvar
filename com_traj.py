import numpy as np

def write_pdb(f,box,atoms,positions,connect,step):
	'''
	Function to write coordinates as a pdb to a file
	'''
	pdb_format = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"
	fixed_format = ['ATOM',
					'INDEX',
					'ATYPE',
					'',
					'',
					'A',
					'INDEX',
					'',
					'POSX',
					'POSY',
					'POSZ',
					1.00,
					1.00,
					'AELE',
					'']

	f.write('CRYST1 {:8.3f}{:8.3f}{:8.3f}{:8.2f}{:8.2f}{:8.2f}\n'.format(*box,90,90,90))
	f.write('MODEL {:>8d}\n'.format(step))
	for i in range(len(atoms)):
		fixed_format[1] = i+1
		fixed_format[2] = atoms[i]
		fixed_format[4] = atoms[i]
		fixed_format[6] = i+1
		fixed_format[8] = positions[i][0]
		fixed_format[9] = positions[i][1]
		fixed_format[10] = positions[i][2]
		fixed_format[13] = atoms[i]
		f.write(pdb_format.format(*fixed_format))
	for pair in connect:
		f.write(f"{'CONECT':6s}{pair[0]:5d}{pair[1]:5d}\n")
	f.write('TER\n')
	f.write('ENDMDL\n')
	return


time = np.loadtxt('A.xvg', comments = ['@', '#'])[:,0]

CoM_A = np.loadtxt('A.xvg', comments = ['@', '#'])[:,1:]
CoM_B = np.loadtxt('B.xvg', comments = ['@', '#'])[:,1:]

CoM_top_A = np.loadtxt('TOP_A.xvg', comments = ['@', '#'])[:,1:]
CoM_bot_A = np.loadtxt('BOT_A.xvg', comments = ['@', '#'])[:,1:]

CoM_top_B = np.loadtxt('TOP_B.xvg', comments = ['@', '#'])[:,1:]
CoM_bot_B = np.loadtxt('BOT_B.xvg', comments = ['@', '#'])[:,1:]

CoM_beta_A = np.loadtxt('BETA_A.xvg', comments = ['@', '#'])[:,1:]
CoM_beta_B = np.loadtxt('BETA_B.xvg', comments = ['@', '#'])[:,1:]


positions = np.concatenate((CoM_A,CoM_B,CoM_top_A,CoM_bot_A,CoM_top_B,CoM_bot_B,CoM_beta_A, CoM_beta_B), axis = 1)*10

fmt = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n"

atoms = ['C','C','O','O','O','O','C','C']
box = [150,150,150]
f = open("com_traj.pdb", 'a')

connect = [(1,7), (2,8), (3,4), (5,6)]

for i in range(len(time)):
    write_pdb(f, box, atoms, positions[i].reshape(len(atoms), 3), connect, i)