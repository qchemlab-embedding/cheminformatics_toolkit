import os
import sys
from pprint import pprint
from pathlib import Path

from .utils_crest import *

class work():

    def __init__(self, options):

        self.options  = options

        self.rundir      = self.options['runinp_dir']
        self.workdir     = self.options['work_pathdir']

        self.out_filepath  = self.options['outfile']
        self.log_filepath  = self.options['logfile']

    def run(self, verbose=True):
        if self.options['job_type'] == 'inptest':
            self.input_test(verbose=verbose)
        else:
            self.hello(verbose=verbose)
            if self.options['job_type'] == 'from_crest_to_pyadf':
                self.from_crest_to_pyadf(verbose=verbose)


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





