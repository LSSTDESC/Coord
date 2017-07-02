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

# This file tests the equivalence of the CelestialCoord functionality with the corresponding
# functionality in the astropy package.

from __future__ import print_function
import numpy as np
import math
import coord
import astropy.coordinates
from astropy import units

from helper_util import *

# We'll use these a lot, so just import them.
from math import pi

@timer
def test_angle():
    """Tests comparing coord.Angle to astropy.coordinates.Angle
    """
    pass


@timer
def test_basic():
    """Basic tests of CelestialCoord construction and corresponding SkyCoord construction.
    """
    pass


@timer
def test_distance():
    """Test calculations of distances on the sphere.
    """
    pass


@timer
def test_precess():
    """Test precession between epochs.
    """
    pass


@timer
def test_galactic():
    """Test the conversion from equatorial to galactic coordinates and back.
    """
    pass


@timer
def test_ecliptic():
    """Test the conversion from equatorial to ecliptic coordinates and back.
    """
    pass


if __name__ == '__main__':
    test_angle()
    test_basic()
    test_distance()
    test_precess()
    test_galactic()
    test_ecliptic()
