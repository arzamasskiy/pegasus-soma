"""
Insert an addintional output in restart file

Last update  2016-03-29  leva@astro.princetone.edu
"""


import numpy as np
import struct
import math
import sys
import getopt
import os
import time
from multiprocessing import Pool


n = 16000 # number of processors (equal to the number of folders)
m = 40    # number of processors to use for changing the restart files

arguments = np.arange(m)

def worker(num):
  for i in range(n//m):
    j = num * (n//m) + i
    filename = "".join(["id",str(j),"/heat-id",str(j),".0002.rst"]) # name of the restart file. change 0000 to apropriate number
    if j == 0:
      filename = "".join(["id0/heat.0002.rst"]) # name of the restart file. change 0000 to apropriate number

    # read the contents of the restart file
    f = open(filename, "r")
    contents = f.readlines()
    f.close()

    # 10th element in contents contains maxout. need to change it
    contents[8] = 'maxout      = 8   \n'
  
    # insert information about the new output (start from line number 88)
    start = 88
    contents.insert(start   , '<output8> \n')
    contents.insert(start+ 1, 'out_fmt = heat                 \n')
    contents.insert(start+ 2, 'dt      = 1.0                  \n')
    contents.insert(start+ 3, 'time    = 2.401000000000000e+03 # Next Output Time\n')
    contents.insert(start+ 4, 'num     = 0                     # Next Output Number\n')
    contents.insert(start+ 5, 'level   = -1                    # Default Value\n')
    contents.insert(start+ 6, 'domain  = -1                    # Default Value\n')
    contents.insert(start+ 7, 'id      = out8                  # Default Value\n')
    contents.insert(start+ 8, 'out     = all                   # Default Value\n')
    contents.insert(start+ 9, 'pargrid = 0                     # Default Value\n')
    contents.insert(start+10, '\n')
  
    # combine and rewrite the restart file

    f = open(filename, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close()

pool = Pool(m)
pool.map(worker,arguments.tolist())