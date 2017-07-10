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
import time
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
    # radians
    theta_coord = 45. * coord.degrees
    theta_astro = astropy.coordinates.Angle(pi/4., units.radian)

    # degrees
    np.testing.assert_almost_equal(theta_coord.rad, theta_astro.rad, decimal=12)
    np.testing.assert_almost_equal(theta_coord / coord.degrees, theta_astro.degree, decimal=12)
    np.testing.assert_almost_equal(theta_coord / coord.hours, theta_astro.hour, decimal=12)
    np.testing.assert_almost_equal(theta_coord / coord.arcmin, theta_astro.arcminute, decimal=12)
    np.testing.assert_almost_equal(theta_coord / coord.arcsec, theta_astro.arcsec, decimal=12)


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
    t0 = time.time()
    c1 = coord.CelestialCoord(0.234 * coord.radians, 0.342 * coord.radians)
    t1 = time.time()
    c2 = c1.precess(2000, 1950)
    c3 = c2.precess(1950, 1900)
    t2 = time.time()

    a1 = astropy.coordinates.SkyCoord(0.234, 0.342, unit=units.radian,
                                      frame=astropy.coordinates.FK5(equinox='J2000'))
    t3 = time.time()
    a2 = a1.transform_to(astropy.coordinates.FK5(equinox='J1950'))
    a3 = a2.transform_to(astropy.coordinates.FK5(equinox='J1900'))
    t4 = time.time()

    np.testing.assert_allclose(c1.rad, [a1.ra.rad, a1.dec.rad], rtol=1.e-5,
                               err_msg='starting coords different')
    np.testing.assert_allclose(c2.rad, [a2.ra.rad, a2.dec.rad], rtol=1.e-5,
                               err_msg='coord/astropy different after 2000->1950')
    np.testing.assert_allclose(c3.rad, [a3.ra.rad, a3.dec.rad], rtol=1.e-5,
                               err_msg='coord/astropy different after 1950->1900')

    print('Compare times for precession calculations:')
    print('  Make CelestialCoords: t = ',t1-t0)
    print('  Make SkyCoords: t = ',t3-t2)
    print('  Precess with Coord: t = ',t2-t1)
    print('  Precess with Astropy: t = ',t4-t3)
    # On my laptop, these times are
    #   Make CelestialCoords: t =  9.10758972168e-05
    #   Make SkyCoords: t =  0.00361394882202
    #   Precess with Coord: t =  0.000560998916626
    #   Precess with Astropy: t =  0.0377740859985

    # Make sure we don't get slow like astropy.  ;)
    # (Travis is a bit slower than the above times, but these limits should still be safe.)
    assert t1-t0 < 0.001, 'Building CelestialCoord is too slow'
    assert t2-t1 < 0.01, 'CelestialCoord.precess is too slow'


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
