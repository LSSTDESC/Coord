from __future__ import print_function
import os
import re

from setuptools import setup

with open('requirements.txt') as f:
    required = f.read().splitlines()

with open('README.rst') as file:
    long_description = file.read()

# Read in the coord version from coord/_version.py
# cf. http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package  # noqa
version_file = os.path.join('coord', '_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    coord_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % version_file)
print('Coord version is %s' % coord_version)

download_url = (
    "https://github.com/LSSTDESC/Coord/releases/tag/v%s.zip" % coord_version
)
dist = setup(
    name="LSSTDESC.Coord",
    version=coord_version,
    author="LSST DESC (contact: Mike Jarvis)",
    author_email="michael@jarvis.net",
    description="Python module for handling angles and celestial coordinates",
    long_description=long_description,
    license="MIT License",
    url="https://github.com/LSSTDESC/Coord",
    download_url=download_url,
    packages=['coord'],
    install_requires=required,
)
