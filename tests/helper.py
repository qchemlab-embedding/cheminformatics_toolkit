import os
from pathlib import Path
from deepdiff import DeepDiff
import filecmp

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
            self.scratch_dir = Path(self.test_space, self.scratch_dirname, test_dir_name).resolve()

            os.makedirs(self.scratch_dir, exist_ok=True)
            self.scratch_space_is_set = True

        if verbose:
            print('scratch space: ')
            print('  - the name of scratch directory: ', self.scratch_dirname)
            print('  - the path to scratch directory: ', self.scratch_dir)
    
 
    def same_files(self, f1, f2):
        p1 = Path(f1).resolve()
        p2 = Path(f1).resolve()
        if p1.exists() and p1.stat().st_size > 0 and p2.exists() and p2.stat().st_size > 0:
            check = filecmp.cmp(f1,f2)
            if check:
                print("{} and {} are the same".format(f1,f2))
            else:
                print("{} and {} are different".format(f1,f2))
            return check





