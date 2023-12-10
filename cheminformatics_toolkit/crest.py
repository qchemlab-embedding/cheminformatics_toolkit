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



class crest_analysis():
    def __init__(self, crest_parentdirname, qm_calcdirname, qm_tpldirname):
        self.crest_parentdirname = crest_parentdirname
        self.qm_parentdirname    = qm_calcdirname
        self.qm_templatedirname  = qm_tpldirname

    def ignore_files(self, dirname, files):
        return [f for f in files if os.path.isfile(os.path.join(dirname, f))]

    def setup_qm_dirstructure(self, dir_template=None, dir_result=None):
        if dir_template is None:
            dir_template = Path(self.crest_parentdirname)
        if dir_result is None:
            dir_result = self.qm_parentdirname
            Path(dir_result).mkdir(parents=True, exist_ok=True)
        shutil.copytree(dir_template, dir_result, ignore=self.ignore_files, dirs_exist_ok=True)

    def get_single_struc(self, finp):
        data = {}
        with open(finp, 'r') as f:
            lines = f.readlines()
            nr_struct = 0
            new_struct_line = 0
            for i, line in enumerate(lines):
                if new_struct_line < len(lines):
                    nr_at = int(lines[new_struct_line])
                    energy = float(lines[new_struct_line+1])
                    coords = [x.strip() for x in lines[new_struct_line+2:new_struct_line+2+nr_at]]
                    nr_struct += 1
                    new_struct_line = new_struct_line + (2 + nr_at)*nr_struct
                    data[nr_struct] = [str(nr_at)] + [str(energy)] + [x for x in coords]
        return data


    def parse_structures(self, input_info=None, output_dir=None):
        files = Path(self.crest_parentdirname).rglob("*.xyz")
        output_files_list = []
        if output_dir is not None:
            qm_parent = output_dir[0]
            qm_coords = output_dir[1]
        else:
            qm_parent = self.qm_parentdirname
            qm_coords = 'coordinates'
        for f in files:
            if f.name == input_info[1]:
                if input_info[0] in f.parts:
                    # TODO this is error prone, be careful
                    molgroup_name = f.parts[1]
                    data = self.get_single_struc(str(f.resolve()))
                    for k, v in data.items():
                        mol_name = 'struc_'+str(k)
                        fileout_name = molgroup_name + '_' + mol_name + '.xyz'
                        fileout_dir  = Path('/'.join([x.replace(self.crest_parentdirname, qm_parent) for x in f.parts[:-1]]), qm_coords)
                        fileout_dir.mkdir(parents=True, exist_ok=True)
                        fileout = Path(fileout_dir, fileout_name).resolve()
                        output_files_list.append(fileout)
                        with open(str(fileout), 'w') as fout:
                            for l in v:
                                fout.write('{}\n'.format(l))
        return output_files_list
                
    def prep_for_pyadf(self, inpfiles, tpl_files):
        for inps in inpfiles:
            calcdir = inps.resolve().parents[1]
            print(calcdir)
            for tpl in tpl_files:
                shutil.copy(Path(self.qm_templatedirname, tpl).resolve(), calcdir)


