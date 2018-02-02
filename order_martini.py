# Lipid order per shell in Martini simulations
# Martin Voegele, 2017-01-10


## Import necessary modules
import argparse
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import MDAnalysis as mda
import geometry

    
def calculate_order(rn): 
    
    # Initialize the list for each tail
    order_d_A = []
    order_d_B = []
    
    # Select the lipids and the CNT
    lipids = u.select_atoms('resname '+rn)

    # Loop through the trajectory
    for ts in u.trajectory[::10]:
        
        # Box dimensions
        box_dim = ts.dimensions   
#        print box_dim
            
        # Go through all lipids
        for res in lipids.residues:  
            
            # Define the two tails of each lipid
            is_in_tail_A = [ res.atoms[i].name[-1]=='A' for i in xrange(len(res.atoms)) ]
            is_in_tail_B = [ res.atoms[i].name[-1]=='B' for i in xrange(len(res.atoms)) ]
            tail_A = res.atoms[is_in_tail_A]
            tail_B = res.atoms[is_in_tail_B]

            # Calculate the local order parameter of each bond in tail A
            for bead in xrange(len(tail_A)-1):
                bond_vec = tail_A[bead+1].position - tail_A[bead].position
                bond_len = np.linalg.norm(bond_vec)
                costheta = bond_vec[2]/bond_len
                order = 0.5*(3.*costheta**2-1.)
                # Append the order parameter to the respective distance range
                order_d_A.append(order)
            
            # Calculate the local order parameter of each bond in tail B
            for bead in xrange(len(tail_B)-1):
                bond_vec = tail_B[bead+1].position - tail_B[bead].position
                bond_len = np.linalg.norm(bond_vec)
                costheta = bond_vec[2]/bond_len
                order = 0.5*(3.*costheta**2-1.)
                # Append the order parameter to the respective distance range
                order_d_B.append(order)
            
        mean_order_A = np.mean(order_d_A)
        std_order_A  = np.std(order_d_A)
        sem_order_A  = np.std(order_d_A)/np.sqrt(len(order_d_A))
        
        mean_order_B = np.mean(order_d_B)
        std_order_B  = np.std(order_d_B)
        sem_order_B  = np.std(order_d_B)/np.sqrt(len(order_d_B))    
            
        return mean_order_A, std_order_A, sem_order_A, mean_order_B, std_order_B, sem_order_B


    

## MAIN ##

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', dest='prmfile', default='system.prmtop', help="topology file (prmtop)")
parser.add_argument('-t', dest='trjfile', default='trajectory.nc', help="trajectory file")
parser.add_argument('-o', dest='outfile', default='./system',      help="directory and prefix for output files")
parser.add_argument('-r', dest='resname', default='POPC', nargs='+', help="residue name of the lipid tails")
args = parser.parse_args()


# Create universe from files
u = mda.Universe(args.prmfile,args.trjfile)

# Loop over all lipid tail species 
for rn in args.resname:
    
    # Calculate all values
    mean_order_A, std_order_A, sem_order_A, mean_order_B, std_order_B, sem_order_B = calculate_order(rn)
    order_A = [mean_order_A, std_order_A, sem_order_A]
    order_B = [mean_order_B, std_order_B, sem_order_B]
    
    # Save the results to files
    header = "Order parameter with Stdev. and SEM"
    np.savetxt( args.outfile+'_order_'+rn+'_A.dat', order_A, fmt='%.8e ', header=header )
    np.savetxt( args.outfile+'_order_'+rn+'_B.dat', order_B, fmt='%.8e ', header=header )
        
exit
