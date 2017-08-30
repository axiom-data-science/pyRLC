#!python
# utf-8
"""Methods for reading .rlc files and converting to .png and .geotiff."""
import os

import numpy as np


class RLCfile:
    """All methods for manipulating .rcl files."""
    def __init__(self, fpath):
        file_path = os.path.abspath(fpath) 

#        with open(file_path, 'rb') as f:
#            (filesig, rlctype, _) = read_header(f)

    def read_header(self, f):
        """Reads and returns header."""
        filesig = np.frombuffer(f.read(4), np.uint8)
        rlctype = np.frombuffer(f.read(1), np.uint8)
        three_zeros = np.frombuffer(f.read(3), np.uint8)

        return (filesig, rlctype, three_zeros)
