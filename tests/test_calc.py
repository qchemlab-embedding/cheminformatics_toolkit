import os
import sys
import shutil
import pytest
from pathlib import Path

import cheminformatics_toolkit as ctk
from .helper import *


def run_test_generic(testdirs, debug=False, verbose=False):

    verbose=True

    th = helper()

    #
    # 1. set paths:
    #    th.test_space = path to the root of test directory
    #    th.testdata_dir = path to the directory of test data
    #
    if not th.test_space_is_set:
        th.set_test_space(verbose=verbose)

    for testdir in testdirs:

        testdir_path = Path(testdir).absolute()

        #
        # 2. set paths to the scratch space for this test
        #
        if not th.scratch_space_is_set:
            th.set_scratch_space(testdir, verbose=verbose)

        os.chdir(th.scratch_dir)

        #
        # 3. find input file
        #
        for tf in os.listdir(testdir_path):
            if tf.endswith('.inp'):
                test_file = Path(testdir_path, tf).absolute()

        #
        # 4. run test
        #
        args = ctk.cli.read_input(finp=test_file, verbose=verbose)

        setup = ctk.cli.input_data(args)
        setup.parse_options()
        setup.process_options()
        setup.print_options()

        calc = ctk.calculate.work(setup.options)
        calc.run(verbose=True)

        #
        # 5. compare output files with reference files
        #
        same = []
        supported_extensions = ['.csv', '.txt']
        for f in os.listdir(Path(testdir_path, 'reference').absolute()):
            for e in supported_extensions:
                if f.endswith(e):
                    f_ref  = Path(testdir_path, 'reference', f).absolute()
                    f_test = Path(th.scratch_dir, f).absolute()
                    if f_ref.exists() and f_test.exists():
                        same.append(th.same_files(f_test, f_ref))

        os.chdir(th.test_space)

        assert (all(x==True for x in same))


def cleanup(testdirs):

    verbose=True

    th = helper()

    if not th.test_space_is_set:
        th.set_test_space()

    for testdir in testdirs:

        if os.path.isdir(testdir):
            testdir_path = Path(testdir).absolute()
            print('CLEANUP: testdir_path = ', testdir_path)

            if not th.scratch_space_is_set:
                th.set_scratch_space(testdir)

            print('CLEANUP: th.scratch_dir = ', th.scratch_dir)
            shutil.rmtree(th.scratch_dir)

            supported_extensions = ['.csv', '.vti']
            for e in supported_extensions:
                for f in os.listdir(testdir_path):
                    if f.endswith(e):
                        f_test = Path(testdir_path, f)
                        os.remove(f_test)


#
#
#def test_test1():
#
#    testdirs = [
#        "test1",
#        ]
#    
#    run_test_generic(testdirs, debug=True)
#    #cleanup(testdirs)
#

def test_inptest():

    testdirs = [
        "inptest",
        ]
    
    run_test_generic(testdirs, debug=True)
    #cleanup(testdirs)


def test_crest():

    testdirs = [
        "crest_test1",
        ]
    
    run_test_generic(testdirs, debug=True)
    #cleanup(testdirs)


#def test_rdkit():
#
#    testdirs = [
#        "rdkit_test1",
#        ]
#    
#    run_test_generic(testdirs, debug=True)
#    #cleanup(testdirs)
