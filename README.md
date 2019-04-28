# ASTrA:
**A**tomistic **S**upercapacitor **Tr**ajectory **A**nalysis

File inputs:
1. input file containing system information and specifying analyses to run
2. xyz trajectory file
3. (optional) LAMMPS-formatted dump file containing charges of atoms
4. elegrps file specifying groups of electrode atoms

Outputs [filename]:
- Radial distribution functions [rdf]
- Average degree of confinement [DoC]: as described in Merlet et al.[2]
- Individual ion degrees of confinement + indices [DoCIndices]
- Coordination numbers [coordnum]: as in Merlet et al. [3]
- Atom counts [atoms]: binned number of atoms (separated by type) throughout production run
- Density [density]: binned density throughout production run
- Ions [ions]: binned number of ions (by type, based on COM) throughout production run
- Ions in layers [layers]: average number of ions in layers (anode, cathode, bulk) throughout production run
- Ions in layers as time series [numionslayers]: number of ions in layers over time throughout full simulation
- Collective variables [ionCV]: charging parameter
- Number of ions in electrode slices [numionselecslices]
- Charge in electrode slices [elecchargeslices]

References:

[1] Liu, Y. M. et al. Mechanisms of Charge Storage in Zeolite-Templated Carbon Supercapacitors. In prep.

[2] Merlet, C. et al. Highly confined ions store charge more efficiently in supercapacitors. Nat. Commun. 4, 2701. issn: 2041-1723 (Dec. 2013).

[3] Merlet, C. et al. On the molecular origin of supercapacitance in nanoporous carbon electrodes. Nat. Mater. 11, 306–310. issn: 1476-1122 (2012).

