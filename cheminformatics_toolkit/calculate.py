import os
import sys
from pprint import pprint
from pathlib import Path

class work():

    def __init__(self, options, fout=None, flog=None):
        self.options  = options
        if fout is not None:
            self.out_filepath  = fout
        else:
            self.out_filepath  = Path(options['work_pathdir'], options['out_filename']).absolute()
        if flog is not None:
            self.log_filepath  = flog
        else:
            self.log_filepath  = Path(options['work_pathdir'], 'log').absolute()

    def run(self, verbose=True):
        if self.options['job_type'] == 'inptest':
            # print to log and stop
            self.input_test(verbose=verbose)
        else:
            self.hello(verbose=verbose)

    def input_test(self, verbose=True):
        print('Input test')
        with open(self.out_filepath, 'w') as f:
            for k, v in self.options.items():
                f.write('{} : {}\n'.format(k, v))
        pprint(self.options)


    def hello(self, verbose=True):
        print('Hello!')
        print('You are running {} with the following arguments:'.format(None))
        pprint(self.options)






