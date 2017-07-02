from __future__ import print_function
import sys,os,glob,re

import setuptools
from setuptools import setup, Extension

print('Python version = ',sys.version)
py_version = "%d.%d"%sys.version_info[0:2]  # we check things based on the major.minor version.

with open('requirements.txt') as f:
    required = f.read().splitlines()

sources = glob.glob(os.path.join('src','*.cpp'))
print('sources = ',sources)
headers = glob.glob(os.path.join('include','*.h'))
print('headers = ',headers)
ext=Extension("coord._coord", sources, depends=headers, include_dirs=['include'])

with open('README.rst') as file:
    long_description = file.read()

# Read in the coord version from coord/_version.py
# cf. http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
version_file=os.path.join('coord','_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    coord_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (version_file,))
print('Coord version is %s'%(coord_version))

dist = setup(name="LSSTDESC.Coord",
      version=coord_version,
      author="LSST DESC (contact: Mike Jarvis)",
      author_email="michael@jarvis.net",
      description="Python module for handling angles and celestial coordinates",
      long_description=long_description,
      license = "MIT License",
      url="https://github.com/LSSTDESC/Coord",
      download_url="https://github.com/LSSTDESC/Coord/releases/tag/v%s.zip"%coord_version,
      packages=['coord'],
      package_data={'coord' : headers },
      ext_modules=[ext],
      install_requires=required)

# setup.py doesn't put the .so file in the coord directory, so this bit makes it possible to
# import coord from the root directory.  Not really advisable, but everyone does it at some
# point, so might as well facilitate it.
lib = os.path.join('coord','_coord.so')
if os.path.lexists(lib): os.unlink(lib)
os.link(glob.glob(os.path.join('build','*',lib))[0], lib)

