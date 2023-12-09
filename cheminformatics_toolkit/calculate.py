import os
import sys
from pprint import pprint
from pathlib import Path

class work():

    def __init__(self, options):
        self.options  = options


    def hello(self, verbose=True):
        print('Hello!')
        print('You are running {} with the following arguments:'.format(None))
        pprint(self.options)






