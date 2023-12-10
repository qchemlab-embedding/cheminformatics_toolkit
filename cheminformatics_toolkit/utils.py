import os
from pathlib import Path
import glob
import shutil
import argparse
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdMolTransforms
from rdkit import rdBase

from rdfpy import rdf, rdf3d

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

try:
    import py3Dmol
except ImportError:
    py3Dmol_available = False
    pass
else:
    py3Dmol_available = True



def plot3Dinteractive(moldict, p=None, removeHs=False):
    """
    return an interactive visualization in 3D
    which works in jupyter notebook
    """

    if p is None:
        p = py3Dmol.view(width=400, height=400)
    p.removeAllModels()

    for key, mol in moldict.items(): 
        #mb = Chem.MolToMolBlock(mol, removeHs=removeHs)
        mb = Chem.MolToMolBlock(mol)
        p.addModel(mb, 'sdf')

    p.setStyle({'stick':{'radius':'0.15'}})
    p.setBackgroundColor('0xeeeeee')
    p.zoomTo()
    return p


class alignment():

    def __init__(self, moldict):
        self.moldict = moldict

    def align_structures_to_lowest_energy_and_show(self, energy_dict, core_smiles):
        """
        align all structures in "moldict" to the one of the lowest energy
        """
        energy_sorted = sorted(energy_dict.items(), key=lambda x: x[1])
        lowest = energy_sorted[0][0]
        core_lowest = moldict[lowest].GetSubstructMatch(Chem.MolFromSmiles(core_smiles))
        
        for key, mol in self.moldict.items():
            match_mol_to_core = mol.GetSubstructMatch(Chem.MolFromSmiles(core_smiles))
            AllChem.AlignMol(mol,self.moldict[lowest],atomMap=list(zip(match_mol_to_core,core_lowest)))
            
        if py3Dmol_available:
            p = plot3Dinteractive(self.moldict)
            return p.show()
    
    
    def align_and_show(self, core_smiles, ref_mol):
        """
        align all structures in "moldict" to a reference structure ("ref_mol"), 
        "core_smiles" provides a molecular pattern to prioritize in this alignment ("core_smiles")
        """ 
        match_ref_to_core = ref_mol.GetSubstructMatch(Chem.MolFromSmiles(core_smiles))
        for key, mol in self.moldict.items():    
            match_mol_to_core = mol.GetSubstructMatch(Chem.MolFromSmiles(core_smiles))
            aligned = AllChem.AlignMol(mol,ref_mol,atomMap=list(zip(match_mol_to_core,match_ref_to_core)))
        
        if py3Dmol_available:
            p = plot3Dinteractive(self.moldict)
            return p.show()


class structural_analysis():
    def __init__(self, moldict):
        self.moldict = moldict

    def get_distances(self, atom_type1, atom_type2):
        dist = {}
        for k, m in self.moldict.items():
            dm = Chem.Get3DDistanceMatrix(m)
            dl=[]
            for at1_ind in atom_type1[k]:
                for at2_ind in atom_type2[k]:
                    dl.append(dm[at1_ind, at2_ind])
            dist[k]=sorted(dl)
        return dist
    
    def get_angles(self, atom_type1, atom_type2, atom_type3):
        angles = {}
        for k, m in self.moldict.items():
            conf=m.GetConformer(0)
            a=[]
            for at1_ind in atom_type1[k]:
                for at2_ind in atom_type2[k]:
                    for at3_ind in atom_type3[k]:
                        angle = rdMolTransforms.GetAngleDeg(conf,at1_ind,at2_ind,at3_ind)
                        # TODO avoid double-counting!
                        a.append(angle)
            angles[k]=sorted(a)
        return angles
    
    def get_rdf(self):
        "uses RDKit"
        rdf={}
        for k, m in self.moldict.items():
            rdf[k]=rdMolDescriptors.CalcRDF(m)
        return rdf
    
    
    def get_rdf_from_rdfpy(self, atom_type1, atom_type2):
        "uses rdfpy"
        rdf={}
        all_at_coords  = []
        for k, m in self.moldict.items():
            # rdfpy works on coordinates passed as a numpy array,
            # so we prepare one that collects coordinates
            # of two selected atom types
            at  = m.GetAtoms()
            at_symbols = [a.GetSymbol() for a in at]
            for i, a in enumerate(at):
                for at1_ind in atom_type1[k]:
                    if at1_ind == i:
                        coords = m.GetConformer().GetAtomPosition(i)
                        all_at_coords.append([coords.x, coords.y, coords.z])
                for at2_ind in atom_type2[k]:
                    if at2_ind == i:
                        coords = m.GetConformer().GetAtomPosition(i)
                        all_at_coords.append([coords.x, coords.y, coords.z])
                #print('TEST k, at c1.x, c1.y, c1.z ', k, at_symbols[i], coords.x, coords.y, coords.z)
            
        coords_array = np.array(all_at_coords)
        g_r, radii = rdf3d(coords_array, dr=0.1)
        return coords_array, g_r, radii 


class dataout():

    """
    This class collects examples of how to save data on files
    """

    def __init__(self, data_dir=None):
        self.data_dir = data_dir


    def data_dict_to_dataframe(self, data_dict):
        """
        entering data_dict is a dict of dict
        """
        result = []
        for p, d in data_dict.items():
            for m, v in d.items():
                for s in v:
                    data = {'parameter':p, 'mol':m, 'value':s}
                    result.append(data)
        data_array = pd.DataFrame(result)
        return data_array

    def write_to_csvfile(self, data_filename=None, dataframe=None):
        if data_filename is not None:
            Path(self.data_dir).mkdir(parents=True, exist_ok=True)
            f=Path(self.data_dir, data_filename).resolve()
            if dataframe is not None:
                dataframe.to_csv(f, index=False)


class plots():

    """
    This class collects examples of plots done with matplotlib
    """

    def __init__(self, plot_dir=None):
        self.plot_dir = plot_dir


    def plot_atoms_as_points(self, g_r, radii, coords=None, plot_filename=None, plot_title=None):
        if coords is not None:
            fig = plt.figure(figsize=(9, 9))
            ax = fig.add_subplot(111, projection='3d')
            ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], color='aqua', alpha=0.3, edgecolors='k', s=25)
            plt.axis('off')
            if plot_filename is not None:
                Path(self.plot_dir).mkdir(parents=True, exist_ok=True)
                f=Path(self.plot_dir, plot_filename).resolve()
                plt.savefig(f, bbox_inches='tight', pad_inches=0.0)
                plt.close()
            else:
                return plt


    def plot_rdf_and_atoms(self, g_r, radii, coords=None, plot_filename=None, plot_title=None):
        if coords is not None:
            tempfile='atoms.png'
            self.plot_atoms_as_points(g_r, radii, coords=coords, plot_filename=tempfile)
            f = Path(self.plot_dir, tempfile).resolve()
            image = plt.imread(f)

            fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={'width_ratios': [1, 2]})

            axes[0].imshow(image)
            axes[0].axis('off')

            axes[1].plot(radii, g_r, color='k', alpha=0.75)
            axes[1].hlines(y=1.0, xmin=0.0, xmax=max(radii), color='r', linestyle='--', alpha=0.4)
            axes[1].set_ylabel(r'g(r)')
            axes[1].set_xlabel(r'r')
            axes[1].set_xlim(0.0, max(radii))
            if plot_title is not None:
                axes[0].set_title(plot_title)
            if plot_filename is not None:
                Path(self.plot_dir).mkdir(parents=True, exist_ok=True)
                f = Path(self.plot_dir, plot_filename).resolve()
                plt.savefig(f, bbox_inches='tight', pad_inches=0.0)
                plt.close()
            else:
                plt.show()
                return plt





