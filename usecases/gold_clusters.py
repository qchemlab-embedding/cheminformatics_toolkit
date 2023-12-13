import os, sys
from pathlib import Path
import glob
import numpy as np

# if testing ctk at the version in directory up, uncomment the line below:
#sys.path.append('../')
import cheminformatics_toolkit as ctk

class gold_clusters():

    def __init__(self, moldict):
        self.moldict = moldict

    def get_ref_for_alignment(self, filepath):
        ref_inp = filepath
        ref_mol = Chem.AddHs(Chem.MolFromMolFile(ref_inp))
        mol_smiles = '[Au]SC.O'
        core_smiles = 'SC'
        return core_smiles, ref_mol

    def find_atoms(self):
        """
        WARNING: 
        this function is specific to systems studied in the project
        and will not be useful in other contexts!
        """
        # in the "environment", we will always have only O and H
        O_env  = {}
        H_env  = {}
        # in the "active" molecule (Au complex), we will always have only Au, S, C, H   
        Au_mol = {}
        S_mol  = {}
        C_mol  = {}
        H_mol  = {}
    
        for k, m in self.moldict.items():
            # we create lists, as we might have more atoms of each type in a molecular system
            O_e = []
            H_e = [] 
            Au_m = []
            S_m = []
            C_m = []
            H_m = []
            
            # now, we collect atom IDs;
            # this atom numbering starts from 0, and (!) will be different than in original xyz/sdf files
            for atom in m.GetAtoms():
                # 1. oxygen atoms - only in environment (H2O):
                if atom.GetAtomicNum() == 8:
                    O_e.append(atom.GetIdx())
    
                # 2. hydrogen atoms:
                if atom.GetAtomicNum() == 1:
                    # if covalently bonded to oxygen, then they "belong" to the environment (H2O)
                    if all(x.GetAtomicNum() == 8 for x in atom.GetNeighbors()):
                        H_e.append(atom.GetIdx())
                    # if covalently bonded to carbon, then they "belong" to the "active" molecule (Au complex)                    
                    elif all(x.GetAtomicNum() == 6 for x in atom.GetNeighbors()):
                        H_m.append(atom.GetIdx())
    
                # gold atoms - only in the "active" molecule (Au complex):
                if atom.GetAtomicNum() == 79:
                    Au_m.append(atom.GetIdx()) 
    
                # sulfur atoms - only in the "active" molecule (Au complex):
                if atom.GetAtomicNum() == 16:
                    S_m.append(atom.GetIdx())
    
                # carbon atoms - only in the "active" molecule (Au complex):
                if atom.GetAtomicNum() == 6:
                    C_m.append(atom.GetIdx())                  
    
                O_env[k] = O_e
                H_env[k] = H_e
                Au_mol[k] = Au_m
                S_mol[k] = S_m
                C_mol[k] = C_m
                H_mol[k] = H_m
     
        return O_env, H_env, Au_mol, S_mol, C_mol, H_mol


verbose=True

# example 1
#test_file='test_crest_pyadf.inp'
#args = ctk.cli.read_input(finp=test_file, verbose=verbose)
#setup = ctk.cli.input_data(args)
#setup.process_options()
#calc = ctk.calculate.work(setup.options)
#calc.run(verbose=True)


# example 2:
test_file='test_analysis.inp'
args = ctk.cli.read_input(finp=test_file, verbose=verbose)
setup = ctk.cli.input_data(args)
setup.process_options()

removeHs=False                                                                                   
name_from_filename=True  
calc = ctk.calculate.work(setup.options)
moldict = ctk.create_moldict_from_sdf_files(calc.coordinates_dir, removeHs=False, name_from_filename=True)
gc=gold_clusters(moldict)
O_env, H_env, Au_mol, S_mol, C_mol, H_mol = gc.find_atoms()
sa = ctk.utils_struct_analysis.struct(moldict)
dist = sa.get_distances(Au_mol, H_env)
#dist = sa.get_distances(Au_mol, H_env)





