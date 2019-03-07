# Copyright (c) 2013-2017 LSST Dark Energy Science Collaboration (DESC)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# The version is stored in _version.py as recommended here:
# http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
from ._version import __version__, __version_info__
import os,cffi,glob

# Also let coord.version show the version.
version = __version__

# Set module level attributes for the include directory and the library file name.
coord_dir = os.path.dirname(__file__)
include_dir = os.path.join(coord_dir,'include')

lib_file = os.path.join(coord_dir,'_coord.so')
# Some installation (e.g. Travis with python 3.x) name this e.g. _coord.cpython-34m.so,
# so if the normal name doesn't exist, look for something else.
if not os.path.exists(lib_file): # pragma: no cover
    alt_files = glob.glob(os.path.join(coord_dir,'_coord*.so'))
    if len(alt_files) == 0:
        raise IOError("No file '_coord.so' found in %s"%coord_dir)
    if len(alt_files) > 1:
        raise IOError("Multiple files '_coord*.so' found in %s: %s"%(coord_dir,alt_files))
    lib_file = alt_files[0]

# Load the C functions with cffi
_ffi = cffi.FFI()
for file_name in glob.glob(os.path.join(include_dir,'*_C.h')):
    with open(file_name) as f:
        _ffi.cdef(f.read())
_lib = _ffi.dlopen(lib_file)

# Explicitly import things that we want to be in the coord namespace
from .angleunit import AngleUnit, arcsec, arcmin, degrees, hours, radians
from .angle import Angle, _Angle
from .celestial import CelestialCoord, _CelestialCoord

# This isn't imported to the coord namespace.  You need to do coord.util.blah
from . import util
