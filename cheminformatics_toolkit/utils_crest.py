import os
from pathlib import Path
import glob
import shutil
import argparse
import numpy as np
import pandas as pd


class crest_analysis():

    def __init__(self, options):

        self.options = options

        self.rundir      = self.options['runinp_dir']
        self.workdir     = self.options['work_pathdir']

        self.out_filepath  = self.options['outfile']
        self.log_filepath  = self.options['logfile']

        self.inpxyzf  = self.options['inpxyzf']
        self.inpmoln  = self.options['molname']
        self.nr_parsed_structures = 0

        self.qm_tmpl_inp = self.options['qm_input_template_fullpath']
        self.qm_tmpl_run = self.options['qm_runscript_template_fullpath']
         

    def ignore_files(self, dirname, files):
        return [f for f in files if os.path.isfile(os.path.join(dirname, f))]


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



    def parse_structures(self):
        fileout_dir= Path(self.workdir, 'tmpdir').absolute()
        fileout_dir.mkdir(parents=True, exist_ok=True)
        output_xyz_files = []
        for fxyz in self.inpxyzf:
            data = self.get_single_struc(str(Path(fxyz).resolve()))
            for k, v in data.items():
                mol_name = self.inpmoln+'_struc_'+str(k)
                fileout  = Path(fileout_dir, mol_name+'.xyz').absolute()
                output_xyz_files.append(fileout)
                with open(str(fileout), 'w') as fout:
                    for l in v:
                        fout.write('{}\n'.format(l))
        return output_xyz_files


    def prep_for_pyadf(self, xyzfiles=None):
        if xyzfiles is None:
            xyzfiles = self.parse_structures()
        self.setup_qm_dirstructure(xyzfiles)


    def setup_qm_dirstructure(self, xyzfiles):
        hs = self.options['hamiltonians']
        ds = self.options['dftfuns']
        bs = self.options['basis_sets']
        for mf in xyzfiles:
            molname = mf.stem
            for h in hs:
                for d in ds:
                    for b in bs:
                        if 'pyadf' in self.options['job_type']:
                            qm_dir = Path(self.workdir, 'qm_calcs', 'pyadf', molname, h, d, b).absolute()
                            Path(qm_dir).mkdir(parents=True, exist_ok=True)
                            shutil.copy(self.qm_tmpl_inp, qm_dir)
                            shutil.copy(self.qm_tmpl_run, qm_dir)


