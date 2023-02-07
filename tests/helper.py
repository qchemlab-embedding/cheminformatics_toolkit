import os
from pathlib import Path

class helper:

    def __init__(self):
        self.test_space = None               # the root of 'tests' directory
        self.testdata_dir = None             # the directory with data used for tests
        self.testinp_dir = None              # the directory with one testcase
        self.scratch_dirname = "scratch"     # name of 'scratch' directory
        self.scratch_dir = None              # path to 'scratch' directory
        self.test_space_is_set = False
        self.scratch_space_is_set = False


    def set_test_space(self, verbose=True):
    
        if not self.test_space_is_set:
            self.test_space = Path(__file__).resolve().parent
            self.testdata_dir = Path(self.test_space, "testdata")
            self.test_space_is_set = True

        if verbose:
            print('test space: ')
            print('  - the root of tests directory:                    ', self.test_space)
            print('  - the directory with data used for tests:         ', self.testdata_dir)
    


    def set_scratch_space(self, test_dir_name, verbose=True):
    
        if not self.scratch_space_is_set:
            self.scratch_dir = Path(self.test_space, self.scratch_dirname, test_dir_name)

            os.makedirs(self.scratch_dir, exist_ok=True)
            self.scratch_space_is_set = True

        if verbose:
            print('scratch space: ')
            print('  - the name of scratch directory: ', self.scratch_dirname)
            print('  - the path to scratch directory: ', self.scratch_dir)
    
 



