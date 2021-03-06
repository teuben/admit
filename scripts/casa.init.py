#
#  this is a piece of code MAC users will need to place in
#  ~/.casa/init.py
#
import os
import sys
from os.path import join

try:
    admit_path = os.environ['ADMIT']
    sys.path.append(admit_path)
    os.environ["PATH"] += os.pathsep + join(admit_path,'bin')
    os.environ["PATH"] += os.pathsep + '/usr/local/bin/'
except KeyError:
    print("ADMIT path not defined. If you wish to use ADMIT, source the admit_start.[c]sh file.")
