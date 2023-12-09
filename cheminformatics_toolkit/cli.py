import os
import sys
import argparse
import pprint
from pathlib import Path
from .global_data import *

class input_data:

    def __init__(self, args_list):

        self.args_list = args_list[1:]
        if args_list[0]:
            self.runinp_dir = Path(args_list[0]).resolve().parent 
        else:
            self.runinp_dir = os.getcwd()

        self.options         = {}

        self.data = global_data.dataset1


    def parse_options(self):

        parser = argparse.ArgumentParser()

        required_args = parser.add_argument_group('required arguments')

        required_args.add_argument('--arg1',
                                   dest='arg1',
                                   action='store',
                                   required=True,
                                   help='''
                                        help msg
                                        ''')


        optional_args = parser.add_argument_group('optional arguments')

        optional_args.add_argument('--arg2',
                                   dest='arg2',
                                   action='append',
                                   choices=self.data,
                                   required=False,
                                   help='''
                                        help msg
                                        ''')

        optional_args.add_argument('--work_pathdir',
                                   dest='work_pathdir',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')

        optional_args.add_argument('--inp_xyz_fullpath',
                                   dest='inp_xyz_fullpath',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--molname',
                                   dest='molname',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--molcharge',
                                   dest='molcharge',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--qm_runscript_template_fullpath',
                                   dest='qm_runscript_template_fullpath',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--qm_input_template_fullpath',
                                   dest='qm_input_template_fullpath',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')

        optional_args.add_argument('--cluster_ntasks',
                                   dest='cluster_ntasks',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--cluster_timeh',
                                   dest='cluster_timeh',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--cluster_part',
                                   dest='cluster_part',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--hamiltonians',
                                   dest='hamiltonians',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--dftfuns',
                                   dest='dftfuns',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        optional_args.add_argument('--basis_sets',
                                   dest='basis_sets',
                                   action='store',
                                   metavar='writeme',
                                   required=False,
                                   help='''
                                        help msg
                                        ''')


        args = parser.parse_args(self.args_list)
        self.options = vars(args)


    def print_options(self):

        if self.options == {}:
            self.parse_options()

        print("Input options:")
        #pprint(self.options)
        print(self.options)


def read_input(finp=None, verbose=False):

    """
    read input options, either from a command line or from a file
    """

    args = []

    if finp is None:
        try:
            args = sys.argv[1:]
        except:
            sys.exit(1)
        if verbose:
            msg1 = 'Arguments for {} are read from command line'.format(sys.argv[0])
            msg2 = '\n'.join(x for x in args[1:])
            print(msg1+'\n'+msg2)
        finp_path = None
    else:
        finp_path = Path(finp)
        with open(finp_path, 'r') as f:
            args = [line.strip() for line in f if line[0] != '#' and line != '\n']
        if verbose:
            msg1 = 'Arguments for {} are read from file {}'.format(sys.argv[0], finp_path)
            msg2 = '\n'.join(x for x in args[1:])
            print(msg1+'\n'+msg2)

    args = [finp_path] + args
            
    return args



if __name__ == '__main__':
    args = read_input()
    data = input_data(args)
    data.parse_options()
    data.print_options()




