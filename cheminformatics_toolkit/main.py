from .cli import *
from .calculate import *

def run(args_file=None, verbose=False):

    """
    main engine of this package
    """

    args = read_input(finp=args_file, verbose=verbose)
    setup = input_data(args)
    calc = work(setup.options)
    calc.run()

if __name__ == "__main__":
    run()
