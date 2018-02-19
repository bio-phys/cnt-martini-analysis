# Analysis of Martini MD Simulations with Carbon Nanotubes

This is a set of analysis tools for the Martini lipid model used in Gromacs. It allows you to calculate the order parameter of the lipids in a pure lipid membrane or with a carbon nanotube (CNT) in it. You can also compute the radial distribution function (RDF) of lipid tails around the CNT, its tilting angle and the motion of its center of mass.

Structure and topology of a CNT for Martini can be generated here: https://github.com/bio-phys/cnt-martini.

## This project contains:

### Analysis of Martini lipid simulations

#### **Order parameter**: *order_martini.py* 
 calculates the coarse-grained order parameter of the lipids (works without CNT).
 
 Usage:
 
    python order_martini.py -p [topology file (.tpr)] -t [trajectory file (.trr/.xtc)] -o [directory and generic name for output files] -r [residue name of the lipid]
for example the following command calculates the order parameter of POPC and saves the results in the folder 'results':

    python order_martini.py -p topol.tpr -t traj.xtc -o ../results/order -r POPC


### Analysis of Martini simulations of a single CNT porin in a lipid membrane

#### **Order by shell**: *order_by_shell_martini.py* 
 calculates the deuterium order parameter in each lipid shell around a carbon nanotube.
 
 Works like order_martini.py, but you need to provide the MDAnalysis selection command (via the option -ref) for the reference molecule (in our case a CNT), for example:

    python order_martini.py -p topol.tpr -t traj.xtc -o ../results/order -r POPC -ref 'resname CNT'

#### **Tilt and center-of-mass motion**: *cntmotion_martini.py* 
 calculates the tilting angle and the motion of the center of mass of a CNT.

#### **Radial lipid density**: *rdf_martini.py* 
 calculates the radial distribution function of acyl chain beads around the central axis of the CNT.


### Required geometry functions
 - geometry.py

## Requirements
 - Python 2.7
 - Python packages: sys, argparse, numpy, scipy, MDAnalysis

## Literature
 - M. Vögele, J. Köfinger, G. Hummer: 
Simulations of Carbon Nanotube Porins in Lipid Bilayers.
Faraday Discussions, 2018 (to be submitted)
