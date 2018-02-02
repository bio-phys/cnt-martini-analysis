# Lipid order per shell in Martini simulations
# Martin Voegele, 2017-01-10


## Import necessary modules
import argparse
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import MDAnalysis as mda
import geometry

    
def calculate_order_per_shell(reference,rn,d_range): 
    
    # Initialize the list for distance ranges
    order_d_A = [[],[],[],[],[],[]]
    order_d_B = [[],[],[],[],[],[]]
    
    # Select the lipids and the CNT
    lipids = u.select_atoms('resname '+rn)
    cnt    = u.select_atoms(reference)

    # Loop through the trajectory
    for ts in u.trajectory[::10]:
        
        # Box dimensions
        box_dim = ts.dimensions   
#        print box_dim
            
        # Calculate COM and principal axis of the CNT
        com_cnt = cnt.center_of_mass()
        pax_cnt = cnt.principal_axes(pbc=True)[0]
    
        # Go through all lipids
        for res in lipids.residues:  
            
            # Define the two tails of each lipid
            is_in_tail_A = [ res.atoms[i].name[-1]=='A' for i in xrange(len(res.atoms)) ]
            is_in_tail_B = [ res.atoms[i].name[-1]=='B' for i in xrange(len(res.atoms)) ]
            tail_A = res.atoms[is_in_tail_A]
            tail_B = res.atoms[is_in_tail_B]

            # Calculate the distance to the CNT
            com_A = tail_A.center_of_mass()
            com_B = tail_B.center_of_mass()
            distA = geometry.dist_point_line_pbc(com_A,com_cnt,pax_cnt,box_dim)
            distB = geometry.dist_point_line_pbc(com_B,com_cnt,pax_cnt,box_dim)
                
            # Calculate the local order parameter of each bond in tail A
            for bead in xrange(len(tail_A)-1):
                bond_vec = tail_A[bead+1].position - tail_A[bead].position
                bond_len = np.linalg.norm(bond_vec)
                costheta = bond_vec[2]/bond_len
                order = 0.5*(3.*costheta**2-1.)
                # Append the order parameter to the respective distance range
                if distA < d_range[0]:
                    order_d_A[0].append(order)
                elif distA < d_range[1]:
                    order_d_A[1].append(order)
                elif distA < d_range[2]:
                    order_d_A[2].append(order)
                elif distA < d_range[3]:
                    order_d_A[3].append(order)
                elif distA < d_range[4]:
                    order_d_A[4].append(order)
                else:
                    order_d_A[5].append(order)    
            
            # Calculate the local order parameter of each bond in tail B
            for bead in xrange(len(tail_B)-1):
                bond_vec = tail_B[bead+1].position - tail_B[bead].position
                bond_len = np.linalg.norm(bond_vec)
                costheta = bond_vec[2]/bond_len
                order = 0.5*(3.*costheta**2-1.)
                # Append the order parameter to the respective distance range
                if distB < d_range[0]:
                    order_d_B[0].append(order)
                elif distB < d_range[1]:
                    order_d_B[1].append(order)
                elif distB < d_range[2]:
                    order_d_B[2].append(order)
                elif distB < d_range[3]:
                    order_d_B[3].append(order)
                elif distB < d_range[4]:
                    order_d_B[4].append(order)
                else:
                    order_d_B[5].append(order)           
            
    mean_order_A = []
    std_order_A  = []
    sem_order_A  = []
    for l in order_d_A:
        mean_order_A.append(np.mean(l))
        std_order_A.append(np.std(l))
        sem_order_A.append(np.std(l)/np.sqrt(len(l)))
            
    mean_order_B = []
    std_order_B  = []
    sem_order_B  = []
    for l in order_d_B:
        mean_order_B.append(np.mean(l))
        std_order_B.append(np.std(l))
        sem_order_B.append(np.std(l)/np.sqrt(len(l)))
            
    return np.array(mean_order_A), np.array(std_order_A), np.array(sem_order_A), np.array(mean_order_B), np.array(std_order_B), np.array(sem_order_B)


    

## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p',   dest='prmfile',   default='martini_md.tpr',       help="topology file")
parser.add_argument('-t',   dest='trjfile',   default='martini_md_whole.xtc', help="trajectory file")
parser.add_argument('-o',   dest='outfile',   default='./system',             help="directory and prefix for output files")
parser.add_argument('-r',   dest='resname',   default='POPC', nargs='+',      help="residue name of the lipid tails")
parser.add_argument('-ref', dest='reference', default='resname CNT',          help="reference to which the distance is calculated")
args = parser.parse_args()


# Distance ranges
d_range = [13.0,17.5,22.0,26.5,31.0]

# Create universe from files
u = mda.Universe(args.prmfile,args.trjfile)

# Loop over all lipid tail species 
for rn in args.resname:
    
    # Calculate all values
    mean_order_A, std_order_A, sem_order_A, mean_order_B, std_order_B, sem_order_B = calculate_order_per_shell(args.reference,rn,d_range)
    order_A = [mean_order_A, std_order_A, sem_order_A]
    order_B = [mean_order_B, std_order_B, sem_order_B]
    
    # Save the results to files
    header = "Order parameters in the ranges between "+str(d_range[0])+", "+str(d_range[1])+", "+str(d_range[2])+", and "+str(d_range[3])+" Angstrom with Stdev. and SEM"
    np.savetxt( args.outfile+'_order_'+rn+'_A.dat', order_A, fmt='%.8e ', header=header )
    np.savetxt( args.outfile+'_order_'+rn+'_B.dat', order_B, fmt='%.8e ', header=header )
        
exit
