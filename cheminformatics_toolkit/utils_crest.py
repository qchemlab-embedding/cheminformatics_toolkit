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
         
        if self.options['nr_atoms_in_active_sub'] is not None:
            self.nr_atoms_in_active = int(self.options['nr_atoms_in_active_sub'])
        else:
            self.nr_atoms_in_active = 0

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



    def parse_structures(self, to_subsystems=False):
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
                if to_subsystems == True:
                    #warning - error-prone:
                    sub1f = Path(fileout_dir, mol_name+'_sub1.xyz').absolute()
                    sub2f = Path(fileout_dir, mol_name+'_sub2.xyz').absolute()
                    with open(str(sub1f), 'w') as fout:
                        fout.write('{}\n'.format(self.nr_atoms_in_active))
                        fout.write('sub1\n')
                        for l in v[2:self.nr_atoms_in_active-1+3]:
                            fout.write('{}\n'.format(l))
                    with open(str(sub2f), 'w') as fout:
                        fout.write('{}\n'.format(int(v[0])-self.nr_atoms_in_active))
                        fout.write('sub2\n')
                        for l in v[2+self.nr_atoms_in_active:]:
                            fout.write('{}\n'.format(l))
        return output_xyz_files


    def cleanup(self, f):
        print('warning: removing dir: ', f)
        shutil.rmtree(f) 

    def prep_for_pyadf(self, xyzfiles=None):
        if xyzfiles is None:
            if self.nr_atoms_in_active == 0:
                xyzfiles = self.parse_structures(to_subsystems=False)
            else:
                xyzfiles = self.parse_structures(to_subsystems=True)
        self.setup_qm_dirstructure(xyzfiles)
        self.cleanup(Path(self.workdir, 'tmpdir').absolute())


    def setup_qm_dirstructure(self, xyzfiles):
        hs = self.options['hamiltonians']
        ds = self.options['dftfuns']
        bs = self.options['basis_sets']
        for mf in xyzfiles:
            molname = mf.stem
            for h in hs:
                h = h.replace(" ", "")
                for d in ds:
                    d = d.replace(" ", "")
                    for b in bs:
                        b = b.replace(" ", "")
                        if 'pyadf' in self.options['job_type']:
                            subdirname=self.options['username']+'-space/qm_calcs/pyadf/' + 'molname/' + h + '/' + d + '/' + b
                            qm_dir = Path(self.workdir, 'qm_calcs', 'pyadf', molname, h, d, b).absolute()
                            coor_dir = Path(self.workdir, 'qm_calcs', 'pyadf', molname, h, d, b, 'coordinates').absolute()
                            Path(qm_dir).mkdir(parents=True, exist_ok=True)
                            Path(coor_dir).mkdir(parents=True, exist_ok=True)
                            shutil.copy(mf, coor_dir)
                            sub1f=str(mf).replace('.xyz', '_sub1.xyz')
                            sub2f=str(mf).replace('.xyz', '_sub2.xyz')
                            print('BUBA: ', sub1f)
                            if Path(sub1f).exists(): 
                                print('BUBA2: ', sub1f)
                                shutil.copy(Path(sub1f), coor_dir)
                            if Path(sub2f).exists(): 
                                shutil.copy(Path(sub2f), coor_dir)
                            shutil.copy(self.qm_tmpl_inp, qm_dir)
                            shutil.copy(self.qm_tmpl_run, qm_dir)
                            # now, string replacement:
                            final_runf = Path(qm_dir, self.qm_tmpl_run.name).absolute()
                            final_inpf = Path(qm_dir, self.qm_tmpl_inp.name).absolute()
                            with open(self.qm_tmpl_run, "r") as f:
                                lines = f.readlines()
                                with open(final_runf, "w") as g:
                                    self.modify_lines(g, lines, \
                                                      data_dir_in_scratch=subdirname, data_dir_in_storage=subdirname, \
                                                      project_name=self.qm_tmpl_inp.stem, \
                                                      molfilename=mf.name, envfilename=mf.name, molcharge=0, \
                                                      cluster_ntasks=self.options['cluster_ntasks'], \
                                                      cluster_timeh=self.options['cluster_timeh'], \
                                                      cluster_part=self.options['cluster_part'], \
                                                      basis=b, hamiltonian=h, dftfun=d)
                            with open(self.qm_tmpl_inp, "r") as f:
                                lines = f.readlines()
                                with open(final_inpf, "w") as g:
                                    self.modify_lines(g, lines, \
                                                      data_dir_in_scratch=subdirname, data_dir_in_storage=subdirname, \
                                                      project_name=self.qm_tmpl_inp.stem, \
                                                      molfilename=mf.name, envfilename=mf.name, molcharge=0, \
                                                      cluster_ntasks=self.options['cluster_ntasks'], \
                                                      cluster_timeh=self.options['cluster_timeh'], \
                                                      cluster_part=self.options['cluster_part'], \
                                                      basis=b, hamiltonian=h, dftfun=d)




    def modify_lines(self, g, lines, \
                     data_dir_in_scratch, data_dir_in_storage, project_name, \
                     molfilename, envfilename, molcharge, \
                     cluster_ntasks=None, cluster_timeh=None, cluster_part=None, \
                     basis=None, hamiltonian=None, dftfun=None):
    
        patterns={}
    
        patterns["cluster_ntasks"] = str(cluster_ntasks)
        patterns["cluster_timeh"] = str(cluster_timeh)
        patterns["cluster_part"] = cluster_part

        patterns["data_dir_in_scratch"] = data_dir_in_scratch
        patterns["data_dir_in_storage"] = data_dir_in_storage
        patterns["project_name"] = project_name
    
        patterns["molfilename"] = molfilename
        patterns["envfilename"] = envfilename
        patterns["molcharge"] = molcharge
    
        patterns["choice_of_basis"] = basis
        patterns["choice_of_dftfun"] = dftfun

        if 'ZORA' in hamiltonian:
            if 'spinorbit' in hamiltonian:
                patterns["choice_of_hamiltonian"] = 'ZORA=True, SpinOrbit=True'
                patterns["choice_of_unrestricted"] = 'remove'
                patterns["choice_of_noncollinear"] = 'True'
            else:
                patterns["choice_of_hamiltonian"] = 'ZORA=True, SpinOrbit=False'
                patterns["choice_of_unrestricted"] = 'False'
                patterns["choice_of_noncollinear"] = 'remove'
        else:
            patterns["choice_of_hamiltonian"] = 'remove'
            patterns["choice_of_unrestricted"] = 'remove'

    
        for l in lines:
            for k, v in patterns.items():
                if v:
                    l = l.replace(k, v)
            if 'remove' not in l:
                g.write(l)
    
