import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import nglview as nv
import copy

# https://biopython.org/docs/1.74/api/Bio.SVDSuperimposer.html
from Bio.SVDSuperimposer import SVDSuperimposer

class Capper:

    def __init__(self,traj,leading_chainid, N_term=True, C_term=False):

        """

        Example Usage:

        traj = md.load('./pdbs/FI_HNS_Atr.pdb')
        capper = Capper(copy.deepcopy(traj),leading_chainid=2,N_term=True,C_term=False)
        print(capper.capped_traj.top)
        for c in capper.capped_traj.top.chains:
        print([res for res in c.residues])


        # Caps have been made as follows:
        peptide = md.load('./capped_peptide.pdb')
        ACE_indices = list(peptide.top.select('(resSeq 0) or (resSeq 1 and backbone)'))
        ACE = peptide.atom_slice(ACE_indices)
        ACE.save('ACE_cap.pdb')
        NME_indices = peptide.top.select('(resid 4) or (resid 3 and backbone)')
        NME = peptide.atom_slice(NME_indices)
        NME.save('NME_cap.pdb')
        """

        self.traj = traj
        self.top = traj.top
        self.leading_chainid = leading_chainid
        self.N_term = N_term
        self.C_term = C_term
    
        self.analyse_leading_chain()
        self.load_caps()
        self.fit_caps()
        self.add_caps()
        self.rebuild_system()

    
    def get_rot_and_trans(self,subtraj_A,subtraj_B):
        """ fit only works now on a single frame (mdtraj returns xyz with shape (n_frames, atoms, xyz) 
            even for single frame trajs so hence the xyz[0]"""
        
        # load super imposer
        sup = SVDSuperimposer()

        # Set the coords, y will be rotated and translated on x
        x = subtraj_A.xyz[0]
        y = subtraj_B.xyz[0]
        sup.set(x, y)

        # Do the leastsquared fit
        sup.run()

        # Get the rms
        rms = sup.get_rms()

        # Get rotation (right multiplying!) and the translation
        rot, tran = sup.get_rotran()
        
        # now we have the instructions to rotate B on A
        return rot,tran,rms

    def apply_superimposition(self,traj, rot, tran):
        
        # get xyz coordinates
        xyz = traj.xyz[0]
        
        # rotate subject on target
        new_xyz = np.dot(xyz, rot) + tran

        # replace coordinates of traj
        traj.xyz = new_xyz
        return traj
    
    def fit_B_on_A(self,A, B, selection_A, selection_B):
        
        # create trajs containing only the selections
        subtraj_A = A.atom_slice(selection_A)
        subtraj_B = B.atom_slice(selection_B)

        # obtain instructions to rotate and translate B on A based on substraj structures
        rot, tran, rms = self.get_rot_and_trans(subtraj_A,subtraj_B)
        
        # do the superimposition of B on A and subsitute old with new xyz of B
        sup_B = self.apply_superimposition(B, rot, tran)

        # remove overlapping backbone atoms from B 
        selection_to_keep = [at.index for at in  B.top.atoms if at.index not in selection_B]

        new_B = sup_B.atom_slice(selection_to_keep)
        return new_B, rms

    def get_fit_indices(self,residue):
        res_indices = [at.index  for at in residue.atoms if at.name in ['N','CA','C','O']]
        return res_indices
            
    def analyse_leading_chain(self):

        # Get terminal residues of leading chain
        chain_top = self.top.chain(self.leading_chainid)
        self.leading_chain_residues = list(chain_top.residues)
        print('First and last residue of leading chain: ', self.leading_chain_residues[0], self.leading_chain_residues[-1])
        self.first_res_fit_indices = self.get_fit_indices(self.leading_chain_residues[0]) 
        self.last_res_fit_indices = self.get_fit_indices(self.leading_chain_residues[-1]) 

    def load_caps(self):
        # Load N and C terminal caps with backbone of the residues one after ACE and one before NME
        if self.N_term:
            self.ACE = md.load('./pdbs/ACE_cap.pdb')
            self.ACE_fit_indices = list(self.ACE.top.select('backbone and not resname ACE'))
        if self.C_term:
            self.NME = md.load('./pdbs/NME_cap.pdb')
            self.NME_fit_indices = list(self.NME.top.select('backbone and not resname NME'))

        if not self.N_term and not self.C_term:
            print('Please specify N_term and/or C_term')

    def fit_caps(self):
        if self.N_term:
            self.fitted_ACE, rms_ACE = self.fit_B_on_A(A=self.traj, B=self.ACE, selection_A=self.first_res_fit_indices, selection_B=self.ACE_fit_indices)
            print(f'RMS of fit ACE: {rms_ACE}')
        if self.C_term:
            self.fitted_NME, rms_NME = self.fit_B_on_A(A=self.traj, B=self.NME, selection_A=self.last_res_fit_indices, selection_B=self.NME_fit_indices)
            print(f'RMS of fit NME: {rms_NME}')
        
    def add_caps(self):

        protein = self.traj.atom_slice(self.traj.top.select(f'chainid {self.leading_chainid}'))
        cap_ace = copy.deepcopy(self.fitted_ACE) if self.N_term else None
        cap_nme = copy.deepcopy(self.fitted_NME) if self.C_term else None

        new_top = md.Topology()
        # set resSeq counter to one residue backwards from first residue of protein
        first_res = self.leading_chain_residues[0]
        s = first_res.resSeq - 1 if self.N_term else first_res.resSeq
        
        # Create empty chain
        c = new_top.add_chain()
        for chain in protein.top.chains:

            # add ACE cap at the beginning if add_ace is True
            if self.N_term:
                for residue in list(cap_ace.top.residues):
                    r = new_top.add_residue(residue.name, c, s, residue.segment_id)
                    s += 1
                    for atom in residue.atoms:
                        new_top.add_atom(atom.name, atom.element, r)
            
            # add the chain residues
            for residue in list(chain.residues):
                r = new_top.add_residue(residue.name, c, s, residue.segment_id)
                s += 1
                for atom in residue.atoms:
                    new_top.add_atom(atom.name, atom.element, r)

            # add NME cap at the end if add_nme is True
            if self.C_term:
                for residue in list(cap_nme.top.residues):
                    r = new_top.add_residue(residue.name, c, s, residue.segment_id)
                    s += 1
                    for atom in residue.atoms:
                        new_top.add_atom(atom.name, atom.element, r)

        # Create bonds
        new_top.create_standard_bonds()

        # Create a list of the xyz attributes to merge, only if the corresponding cap is not None
        xyzs_to_merge = []
        if cap_ace is not None:
            xyzs_to_merge.append(cap_ace.xyz)
        xyzs_to_merge.append(protein.xyz)
        if cap_nme is not None:
            xyzs_to_merge.append(cap_nme.xyz)
        
        # Merge the xyzs and create a new trajectory
        merged_xyz = np.concatenate(xyzs_to_merge, axis=1)
        self.capped_protein = md.Trajectory(merged_xyz, new_top)

    def rebuild_system(self):
        # Create a new empty topology
        new_top = md.Topology()

        # Initialize a list to store new coordinates
        new_coords = []

        # Iterate over the chains in the original trajectory
        for chain in self.traj.top.chains:
            new_chain = new_top.add_chain()
            # If it's the chain to cap, add the capped protein
            if chain.index == self.leading_chainid:
                for residue in self.capped_protein.top.residues:
                    new_res = new_top.add_residue(residue.name, new_chain, resSeq=residue.resSeq, segment_id=residue.segment_id)
                    for atom in residue.atoms:
                        new_top.add_atom(atom.name, atom.element, new_res)
                new_coords.append(self.capped_protein.xyz)
            else:
                # If it's not the chain to cap, copy it over
                sub_traj = self.traj.atom_slice(self.traj.top.select(f'chainid {chain.index}'))
                for residue in chain.residues:
                    new_res = new_top.add_residue(residue.name, new_chain, resSeq=residue.resSeq, segment_id=residue.segment_id)
                    for atom in residue.atoms:
                        new_top.add_atom(atom.name, atom.element, new_res)
                new_coords.append(sub_traj.xyz)
        
        # Create bonds
        new_top.create_standard_bonds()

        # Concatenate all the new coordinates
        new_xyz = np.concatenate(new_coords, axis=1)

        # This will become the new trajectory
        self.capped_traj = md.Trajectory(new_xyz, new_top)