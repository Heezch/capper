# capper
"I have pdb file with multiple chains. Any easy to run script that can go over chains and cap each end with ACE/NME without causing substantial clashes?"
"Is there a scripted tool to add neutral terminal caps (e.g. NME, ACE) to peptide chains for use with AMBER forcefields?" 

Yes! Look no further.


Requires the following packages:
- mdtraj: https://mdtraj.org/1.9.3/installation.html
- biopython: https://biopython.org/wiki/Download
- numpy
- copy


Example:

`traj = md.load('./pdbs/uncapped_protein.pdb')`
`capper = Capper(copy.deepcopy(traj),leading_chainid=2,N_term=True,C_term=False)`

Uses backbone of terminal residues and those of termini to fit the caps on the protein using SVD on atoms.
