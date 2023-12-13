import os
import sys
from pprint import pprint
from pathlib import Path

from .utils_crest import *
from .utils_struct_analysis import *

class work():

    def __init__(self, options):

        self.options  = options

        self.rundir      = self.options['runinp_dir']
        self.workdir     = self.options['work_pathdir']

        self.out_filepath  = self.options['outfile']
        self.log_filepath  = self.options['logfile']

        if self.options['inpsdfdir'] is not None:
            self.coordinates_dir = self.options['inpsdfdir']
        else:
            self.coordinates_dir = None
 

    def run(self, moldict=None, verbose=True):
        if self.options['job_type'] == 'inptest':
            self.input_test(verbose=verbose)
        else:
            self.hello(verbose=verbose)
            if self.options['job_type'] == 'from_crest_to_pyadf':
                self.from_crest_to_pyadf(verbose=verbose)
            if self.options['job_type'] == 'structural_analysis':
                self.structural_analysis(moldict=moldict, verbose=verbose)


    def input_test(self, verbose=True):
        print('Hello!')
        print('You are testing an input (print setup and stop).')
        with open(self.out_filepath, 'w') as f:
            for k, v in self.options.items():
                f.write('{} : {}\n'.format(k, v))
        pprint(self.options)


    def hello(self, verbose=True):
        print('Hello!')
        print('You are running cheminformatics_toolkit with the following arguments:')
        pprint(self.options)


    def from_crest_to_pyadf(self, verbose=True):
        print('Hello!')
        print('You are running cheminformatics_toolkit: input data from CREST, generated output for pyADF.')
        with open(self.out_filepath, 'w') as f:
            for k, v in self.options.items():
                f.write('{} : {}\n'.format(k, v))
        pprint(self.options)
        ca = crest_analysis(self.options)
        ca.prep_for_pyadf()


    def structural_analysis(self, moldict=None, verbose=True):
        print('Hello!')
        print('You are running cheminformatics_toolkit: input data - molecular structure(s), output data - structural information.')
        with open(self.out_filepath, 'w') as f:
            for k, v in self.options.items():
                f.write('{} : {}\n'.format(k, v))
        pprint(self.options)
     
        #for inpxyzf in self.options['inpxyzf']:
        #    coordinates_dir = Path(inpxyzf).resolve().parent

        if moldict is None:
            removeHs=False
            name_from_filename=True 
            if self.coordinates_dir is not None:
                moldict = create_moldict_from_sdf_files(self.coordinates_dir, removeHs=removeHs, name_from_filename=name_from_filename)

        # WIP; thus far see "usecases"
        #sa = struct(moldict) #self.options)
        #for p in self.options['structural_parameters']:
        #    if p == 'distances':
        #        dist = sa.get_distances(self.options['at1'], self.options['at1'])




