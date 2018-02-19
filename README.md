# Analysis of Martini MD Simulations with Carbon Nanotube Porins (or membrane proteins)

This is a set of analysis tools for the Martini lipid model used in Gromacs. 
It was developed and optimized for use with carbon nanotubes, but can in many cases be used for membrane proteins as well.

It allows you to calculate the average order parameter of the lipids in any Martini lipid simulation. You can also compute the radial distribution function (RDF) of lipid tails around the CNT, its tilting angle and the motion of its center of mass.

Structure and topology of a CNT for Martini can be generated here: https://github.com/bio-phys/cnt-martini.


## Requirements
 - Python 2.7
 - Python packages: sys, argparse, numpy, scipy, MDAnalysis


## Installation
 - no installation needed. 


## Analysis of any Martini lipid simulations

### **Order parameter**: *order_martini.py* 
 calculates the coarse-grained order parameter of the lipids (works without CNT).
 
 Usage:
 
    python order_martini.py -p [topology file (.tpr)] -t [trajectory file (.trr/.xtc)] -o [directory and generic name for output files] -r [residue name of the lipid]
for example the following command calculates the order parameter of POPC and saves the results in the folder 'results':

    python order_martini.py -p topol.tpr -t traj.xtc -o ../results/order -r POPC



## Analysis of Martini simulations of a single CNT porin in a lipid membrane


### **Order by shell**: *order_by_shell_martini.py* 
 calculates the deuterium order parameter in each lipid shell around a carbon nanotube.
 
 Works like order_martini.py, but you need to provide the MDAnalysis selection command (via the option -ref) for the reference molecule (in our case a CNT), for example:

    python order_martini.py -p topol.tpr -t traj.xtc -o ../results/order -r POPC -ref 'resname CNT'


### **Tilt and center-of-mass motion**: *cntmotion_martini.py* 
 calculates the cosine of the tilting angle and the motion of the center of mass (COM) of a CNT (or any other molecule you select).
 
 Usage:
   
    python cntmotion_martini.py -p [topology file (.tpr)] -t [trajectory file (.trr/.xtc)] -ocom [output file for COM] -otil [output file for cos of tilt angle] -sel [selection command]
    
Example:

    python cntmotion_martini.py -p topol.tpr -t traj.xtc -ocom ../results/com.dat -otil ../results/tilt.dat -sel 'resname CNT'


### **Radial lipid density**: *rdf_martini.py* 
 calculates the radial distribution function of selected beads (e.g. acyl chain beads) around the central axis of the reference molecule (e.g. a CNT)
 
Usage:
   
    python rdf_martini.py -p [topology file (.tpr)] -t [trajectory file (.trr/.xtc)] -o [output file] -sel [selection of the molecules for the density] -ref [selection for the reference molecule]
    
Example:

    python rdf_martini.py -p topol.tpr -t traj.xtc -o ../results/rdf.dat -sel 'resname POPC and name C* D*' -ref 'resname CNT'


## Required geometry functions
 - geometry.py

## Literature
 - M. Vögele, J. Köfinger, G. Hummer: 
Simulations of Carbon Nanotube Porins in Lipid Bilayers.
Faraday Discussions, 2018 (to be submitted)
