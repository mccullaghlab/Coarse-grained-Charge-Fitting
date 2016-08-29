#USAGE :  python opt_cg_charges.py [config file name]

#DEPENDCIES : numpy, MDAnalysis, sys

#CONFIG FILE FORMAT:
#atomtopfile = [atom topology file]
#atomcoordfile = [atom coordinate file]
#cgcoordfile = [CG coordinate file]
#cgtopfile = [CG topology file]
#outfile = [Output pdb trajectory file - will contain CG positions and charges for every minimization step]
#thresh = [minimization threshold]
#minsteps = [maximum number of minimization steps]
#
# NOTE ON FILE FORMATS: input topology and coordinate file formats are those allowed by MDAnalysis.  One can input the same pdb file for both.  The topolgy file of 
# of the atomistic system MUST contain charges.
#
# CITATION: P. McCullagh, P. T. Lake, M. McCullagh. Deriving Coarse-grained Charges from All-atom Systems: An Analytic Solution. J. Chem. Theory Comp., 2016

# Python Libraries
import numpy as np
import MDAnalysis
import sys

pi = 3.1415926535

##############################################################################################################
#######################################          SUBROUTINES          ########################################
##############################################################################################################

# read the configuration file and populate the global variables
def ParseConfigFile(cfg_file):
	global atom_top_file, atom_coord_file, cg_top_file, cg_coord_file, out_file, thresh, min_steps
	f = open(cfg_file)
	for line in f:
		# first remove comments
		if '#' in line:
			line, comment = line.split('#',1)
		if '=' in line:
			option, value = line.split('=',1)
			option = option.strip()
			value = value.strip()
			print "Option:", option, " Value:", value
			# check value
			if option.lower()=='atomtopfile':
				atom_top_file = value
			elif option.lower()=='atomcoordfile':
				atom_coord_file = value
			elif option.lower()=='cgtopfile':
				cg_top_file = value
			elif option.lower()=='cgcoordfile':
				cg_coord_file = value
			elif option.lower()=='outfile':
				out_file = value
			elif option.lower()=='thresh':
				thresh = float(value)
			elif option.lower()=='minsteps':
				min_steps = int(value)
			else :
				print "Option:", option, " is not recognized"
	f.close()

# compute the squared distance between two points 
def dist2(r1,r2):

	dist2 = 0

	for j in range(3):
		temp = r1[j]-r2[j]
		dist2 += temp*temp

	return dist2;

# compute the pair negative distance matrix between sites in a system given positions in an Nx3 array
def symPairDistMat(coord):

	rows = coord.shape[0]
	mat = np.empty((rows,rows),dtype=float)
		
	for i in range(rows-1):
		# set diagonal elements to zero
		mat[i,i] = 0.0
		for j in range(i+1,rows):
			# populate with negative distance between two sites
			mat[i,j] = - np.sqrt(dist2(coord[i,:],coord[j,:]))
			# symmetrize
			mat[j,i] = mat[i,j]
	# set final diagonal element to zero
	mat[rows-1,rows-1] = 0.0
	# return the symmetric pair distance matrix
	return mat;

# compute pair negative distance matrix between two sets of coordinates
def pairDistMat(coord1,coord2):

	rows = coord1.shape[0]
	columns = coord2.shape[0]
	mat = np.empty((rows,columns),dtype=float)
		
	for i in range(rows):
		for j in range(columns):
			# populate with negative distance between two sites
			mat[i,j] = - np.sqrt(dist2(coord1[i,:],coord2[j,:]))
	# return matrix
	return mat;

# fit the CG charges using procedure outlined in the paper
def fitCgCharges(A,B,atomic_charges):

	# get size of A = number of CG sites
	n_cg = B.shape[0]
	# get columns of B = number of atoms
	n_atoms = A.shape[1]

	# create D matrix
	D = np.zeros((n_cg-1,n_cg),dtype=float)
	for i in range(n_cg-1):
		D[i,i] = 1.0
		D[i,i+1] = -1.0
	
	# create temporary A and B matricies that we will manipulate
	A_temp = np.empty((n_cg,n_atoms),dtype=float)
	B_temp = np.empty((n_cg,n_cg),dtype=float)
	# the temporary matrices will be manipulated to remove infinite components
	A_temp[0:n_cg-1,:] = np.dot(D,A)
	A_temp[n_cg-1,:] = 1.0
	B_temp[0:n_cg-1,:] = np.dot(D,B)
	B_temp[n_cg-1,:] = 1.0

	solution = np.dot(A_temp,atomic_charges)
	cg_charges = np.linalg.solve(B_temp,solution)
#	print "Total charge of optimized CG charges:", sum(cg_charges)
	return cg_charges


# compute the ESP matching residual
def espResidual(A,B,atomic_charges,cg_charges,atomic_integral):

	rss = atomic_integral - 2*np.dot(cg_charges.T,np.dot(A,atomic_charges)) + np.dot(cg_charges.T,np.dot(B,cg_charges))
	return 2*pi*rss

# compute the spacial gradient of the ESP residual
def charge_gradient(cg_coord,cg_charges,atomic_coord,atomic_charges):

	n_cgs = cg_coord.shape[0]
	n_atoms = atomic_coord.shape[0]

	charge_force = np.zeros((n_cgs,3),dtype=float)

	for cg1 in range(n_cgs):
		for atom in range(n_atoms):
			temp = cg_coord[cg1,:]-atomic_coord[atom,:]
			temp /= np.linalg.norm(temp)
			temp = -4 * pi * cg_charges[cg1] * atomic_charges[atom] * temp
			charge_force[cg1,:] += temp
		for cg2 in range(n_cgs):
			if (cg1 != cg2) :
				temp = cg_coord[cg1,:]-cg_coord[cg2,:]
				temp /= np.linalg.norm(temp)
				temp = 4 * pi * cg_charges[cg1] * cg_charges[cg2] * temp
				charge_force[cg1,:] += temp

	return charge_force

# compute the RMSD between two coordinate sets
def rmsd(coord1, coord2):

	n_cgs = coord1.shape[0]
	rmsd_val = 0.0
	for atom in range(n_cgs):
		rmsd_val += dist2(coord1[atom,:],coord2[atom,:])
	rmsd_val = np.sqrt(rmsd_val/float(n_cgs))
	return rmsd_val


##############################################################################################################
#######################################         MAIN PROGRAM          ########################################
##############################################################################################################


# read in command line argument
cfg_file = sys.argv[1]

# read cfg file
ParseConfigFile(cfg_file)

print "Atomic topology file:", atom_top_file
print "Atomic coordinate file:", atom_coord_file
print "Coarse-grained topology file:", atom_coord_file
print "Initial coarse-grained coordinate file:", cg_coord_file
print "Coordinate trajectory will be written to:", out_file

# initiate MDAnalysis coordinate universe for atomic system
atomic = MDAnalysis.Universe(atom_top_file, atom_coord_file)
# select all atoms
atomic_sel = atomic.select_atoms("all")
# declare array that will contain the selected atom positions at each time frame
atom_positions = np.empty((atomic_sel.n_atoms,3),dtype=float)
# initiate MDAnalysis coordinate universe for coarse-grained system
cg = MDAnalysis.Universe(cg_top_file, cg_coord_file)
# select all CG particles
cg_sel = cg.select_atoms("all")
# declare array that will contain the CG positions
cg_positions = np.empty((cg_sel.n_atoms,3),dtype=float)

# print total charge
print "Total charge on atomic system:", atomic_sel.total_charge()

# center
#atomic_sel.translate(-atomic_sel.center_of_mass())
#cg_sel.translate(-atomic_sel.center_of_mass())
#atomic_sel.write("atomic_centered.pdb")

C = symPairDistMat(atomic_sel.coordinates())
atomic_integral = np.dot(atomic_sel.charges.T,np.dot(C,atomic_sel.charges))

A = pairDistMat(cg_sel.coordinates(),atomic_sel.coordinates())
B = symPairDistMat(cg_sel.coordinates())
# compute the optimal CG charges given these positions
cg_charges = fitCgCharges(A,B,atomic_sel.charges)
# compute the residual
rss = espResidual(A,B,atomic_sel.charges,cg_charges,atomic_integral)
print "Step:", 0, " Residual = ", rss

# open pdb file for writing output
pdbtrj = out_file
lbda = 0.01 # minimization step size
old_positions = np.empty((cg_sel.n_atoms,3),dtype=float)
with MDAnalysis.Writer(pdbtrj, multiframe=True, bonds=False, n_atoms=cg_sel.n_atoms) as PDB:
	cg_sel.set_bfactors(cg_charges)
	PDB.write(cg_sel)
	for min_step in range(min_steps):
		
		# compute charge gradient
		cg_charge_force = charge_gradient(cg_sel.coordinates(),cg_charges,atomic_sel.coordinates(),atomic_sel.charges)
		# save old positions
		old_positions = cg_sel.coordinates()
		# shift cg positions
		cg_sel.positions += lbda*cg_charge_force
		# recompute pair distance matrices
		A = pairDistMat(cg_sel.coordinates(),atomic_sel.coordinates())
		B = symPairDistMat(cg_sel.coordinates())
		# compute the optimal CG charges given these positions
		cg_charges = fitCgCharges(A,B,atomic_sel.charges)
		cg_sel.set_bfactors(cg_charges)
		PDB.write(cg_sel)
		
		# compute the residual
		previous_rss = rss
		rss = espResidual(A,B,atomic_sel.charges,cg_charges,atomic_integral)
		delta = rmsd(old_positions,cg_sel.positions)
		print "Step:", min_step+1, " Residual = ", rss, " position delta = ", delta
		if (delta < thresh):
			print "The residual has converged"
			sys.exit(0)


# print final positions and charges
sys.stdout.write("Current CG positions and charges:\n")
sys.stdout.write(("%20s%20s%20s%20s\n") % ("X","Y","Z","Charge"))
for i in range(cg_sel.n_atoms):
	sys.stdout.write(("%20.10f%20.10f%20.10f%20.10f\n") % (cg_sel.positions[i,0],cg_sel.positions[i,1],cg_sel.positions[i,2],cg_charges[i]) )
sys.stdout.write("Please cite: P. McCullagh, P. T. Lake, M. McCullagh. Deriving Coarse-grained Charges from All-atom Systems: An Analytic Solution. J. Chem. Theory Comp., 2016.\n")


