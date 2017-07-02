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
    # astropy has more ways to create an Angle, but some are quite analogous to coord.Angle.
    a1 = astropy.coordinates.Angle(10.2345 * units.deg)
    c1 = 10.2345 * coord.degrees
    np.testing.assert_almost_equal(c1.rad, a1.rad, decimal=12)

    a2 = astropy.coordinates.Angle(23.09, units.arcsec)
    c2 = coord.Angle(23.09, coord.arcsec)
    np.testing.assert_almost_equal(c2.rad, a2.rad, decimal=12)

    a3 = astropy.coordinates.Angle(-0.17, unit='rad')
    c3 = coord._Angle(-0.17)
    np.testing.assert_almost_equal(c3.rad, a3.rad, decimal=12)

    # astropy wrapping uses a different convention than we do.  Their argument is
    # the upper end of the target range, not the center.
    a4 = a3.wrap_at(360 * units.deg)
    c4 = c3.wrap(180 * coord.degrees)
    np.testing.assert_almost_equal(c4.rad, a4.rad, decimal=12)

    a5 = a3.wrap_at(-100 * units.deg)
    c5 = c3.wrap(-280 * coord.degrees)
    np.testing.assert_almost_equal(c5.rad, a5.rad, decimal=12)

    a6 = astropy.coordinates.Angle('03:34:12', unit='hourangle')
    c6 = coord.Angle.from_hms('03:34:12')
    np.testing.assert_almost_equal(c6.rad, a6.rad, decimal=12)

    a7 = astropy.coordinates.Angle('03:34:12', unit='deg')
    c7 = coord.Angle.from_dms('03:34:12')
    np.testing.assert_almost_equal(c7.rad, a7.rad, decimal=12)

    # Their default arguments to to_string are different from ours, but can make them compatible.
    print('a6.hms = ',a6.to_string(sep=':', pad=True))
    print('c6.hms = ',c6.hms())
    assert c6.hms() == a6.to_string(sep=':', pad=True)

    print('a7.dms = ',a7.to_string(sep=':', pad=True))
    print('c7.dms = ',c7.dms())
    assert c7.dms() == a7.to_string(sep=':', pad=True)

    print('a6.hms = ',a6.to_string())
    print('c6.hms = ',c6.hms(sep='hms', pad=False))
    assert c6.hms(sep='hms', pad=False) == a6.to_string()

    print('a7.hms = ',a7.to_string())
    print('c7.hms = ',c7.dms(sep='dms', pad=False))
    assert c7.dms(sep='dms', pad=False) == a7.to_string()


@timer
def test_basic():
    """Basic tests of CelestialCoord construction and corresponding SkyCoord construction.
    """
    a1 = astropy.coordinates.SkyCoord(10.625 * units.degree, 41.2 * units.degree)
    c1 = coord.CelestialCoord(10.625 * coord.degrees, 41.2 * coord.degrees)
    np.testing.assert_almost_equal(c1.ra.rad, a1.ra.rad, decimal=12)
    np.testing.assert_almost_equal(c1.dec.rad, a1.dec.rad, decimal=12)

    a2 = astropy.coordinates.SkyCoord('00:42.5 +41:12', unit=(units.hourangle, units.deg))
    c2 = coord.CelestialCoord(coord.Angle.from_hms('00:42:30'), coord.Angle.from_dms('41:12:00'))
    np.testing.assert_almost_equal(c2.ra.rad, a2.ra.rad, decimal=12)
    np.testing.assert_almost_equal(c2.dec.rad, a2.dec.rad, decimal=12)

    # astropy made different choices for the HMS and DMS string representations, but we can
    # make it match our choice without too much trouble.
    astr = a2.to_string(style='hmsdms', sep=':')
    print('astr = ',astr)
    astr_ra, astr_dec = astr.split()
    cstr_ra = c2.ra.hms()
    cstr_dec = c2.dec.dms(plus_sign=True)
    assert astr_ra == cstr_ra
    assert astr_dec == cstr_dec


@timer
def test_distance():
    """Test calculations of distances on the sphere.
    """
    t0 = time.time()
    c1 = coord.CelestialCoord(0.234 * coord.radians, 0.342 * coord.radians)
    c2 = coord.CelestialCoord(0.234 * coord.radians, -1.093 * coord.radians)
    c3 = coord.CelestialCoord((pi + 0.234) * coord.radians, -0.342 * coord.radians)
    c4 = coord.CelestialCoord((pi + 0.234) * coord.radians, 0.832 * coord.radians)
    c5 = coord.CelestialCoord(1.832 * coord.radians, -0.723 * coord.radians)
    c6 = coord.CelestialCoord((0.234 + 2.3e-9) * coord.radians, (0.342 + 1.2e-9) * coord.radians)
    t1 = time.time()

    a1 = astropy.coordinates.SkyCoord(0.234 * units.radian, 0.342 * units.radian)
    a2 = astropy.coordinates.SkyCoord(0.234 * units.radian, -1.093 * units.radian)
    a3 = astropy.coordinates.SkyCoord((pi + 0.234) * units.radian, -0.342 * units.radian)
    a4 = astropy.coordinates.SkyCoord((pi + 0.234) * units.radian, 0.832 * units.radian)
    a5 = astropy.coordinates.SkyCoord(1.832 * units.radian, -0.723 * units.radian)
    a6 = astropy.coordinates.SkyCoord(0.234 + 2.3e-9, 0.342 + 1.2e-9, unit=units.radian)
    t2 = time.time()

    coord_dist = [c1.distanceTo(c).rad for c in [c2,c3,c4,c5,c6]]
    t3 = time.time()
    astropy_dist = [a1.separation(a).rad for a in [a2,a3,a4,a5,a6]]
    t4 = time.time()

    np.testing.assert_almost_equal(coord_dist, astropy_dist, decimal=12)
    # For the last one, the distance is rather small in radians, so test in arcsec
    np.testing.assert_almost_equal(coord_dist[-1] * (coord.radians/coord.arcsec),
                                   astropy_dist[-1] * (coord.radians/coord.arcsec), decimal=10)

    print('Compare times for distance calculations:')
    print('  Make CelestialCoords: t = ',t1-t0)
    print('  Make SkyCoords: t = ',t2-t1)
    print('  Calculate distances with Coord: t = ',t3-t2)
    print('  Calculate distances with Astropy: t = ',t4-t3)


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

    np.testing.assert_almost_equal(c1.ra.rad, a1.ra.rad, decimal=6)
    np.testing.assert_almost_equal(c1.dec.rad, a1.dec.rad, decimal=6)
    np.testing.assert_almost_equal(c2.ra.rad, a2.ra.rad, decimal=6)
    np.testing.assert_almost_equal(c2.dec.rad, a2.dec.rad, decimal=6)
    np.testing.assert_almost_equal(c3.ra.rad, a3.ra.rad, decimal=6)
    np.testing.assert_almost_equal(c3.dec.rad, a3.dec.rad, decimal=6)

    print('Compare times for precession calculations:')
    print('  Make CelestialCoords: t = ',t1-t0)
    print('  Make SkyCoords: t = ',t3-t2)
    print('  Precess with Coord: t = ',t2-t1)
    print('  Precess with Astropy: t = ',t4-t3)


@timer
def test_galactic():
    """Test the conversion from equatorial to galactic coordinates and back.
    """
    center = coord.CelestialCoord(coord.Angle.from_hms('17:45:37.1991'),
                                  coord.Angle.from_dms('-28:56:10.2207'))
    north = coord.CelestialCoord(coord.Angle.from_hms('12:51:26.27549'),
                                 coord.Angle.from_dms('27:07:41.7043'))
    south = coord.CelestialCoord(coord.Angle.from_hms('00:51:26.27549'),
                                 coord.Angle.from_dms('-27:07:41.7043'))
    anticenter = coord.CelestialCoord(coord.Angle.from_hms('05:45:37.1991'),
                                      coord.Angle.from_dms('28:56:10.2207'))
    random = coord.CelestialCoord(0.234 * coord.radians, 0.342 * coord.radians)

    for c1 in [center, north, south, anticenter, random]:
        a1 = astropy.coordinates.SkyCoord(c1.ra.rad, c1.dec.rad, unit=units.rad, frame='fk5')
        #print('c1.galactic() = ',c1.galactic())
        #print('a1.galactic = ',a1.galactic)
        el, b = c1.galactic()
        # Wrap el to the same phase as a1
        el = el.wrap(a1.galactic.l.rad * coord.radians)
        if c1 not in [north, south]:
            np.testing.assert_almost_equal(el.rad, a1.galactic.l.rad, decimal=6)
        np.testing.assert_almost_equal(b.rad, a1.galactic.b.rad, decimal=6)

        c2 = coord.CelestialCoord.from_galactic(el,b)
        a2 = astropy.coordinates.SkyCoord(el.rad, b.rad, unit=units.radian, frame='galactic')
        a2 = a2.transform_to('fk5')
        c2_ra = c2.ra.wrap(a2.ra.rad * coord.radians)
        np.testing.assert_almost_equal(c2_ra.rad, a2.ra.rad, decimal=6)
        np.testing.assert_almost_equal(c2.dec.rad, a2.dec.rad, decimal=6)

@timer
def test_ecliptic():
    """Test the conversion from equatorial to ecliptic coordinates and back.
    """
    # Note: These don't have the same level of agreement as galactic coords did.
    # Possibly this is related to a caveat in the astropy documentation:
    #     "In the current version of astropy, the ecliptic frames do not yet have stringent
    #      accuracy tests."
    # Or possibly, the Coord conversions aren't super accurate.  I'm not really sure.
    north = coord.CelestialCoord(coord.Angle.from_hms('18:00:00.00'),
                                 coord.Angle.from_dms('66:33:38.55'))
    south = coord.CelestialCoord(coord.Angle.from_hms('06:00:00.00'),
                                 coord.Angle.from_dms('-66:33:38.55'))
    vernal = coord.CelestialCoord(0.*coord.radians, 0.*coord.radians)
    autumnal = coord.CelestialCoord(pi*coord.radians, 0.*coord.radians)
    random = coord.CelestialCoord(0.234 * coord.radians, 0.342 * coord.radians)

    for c1 in [north, south, vernal, autumnal, random]:
        a1 = astropy.coordinates.SkyCoord(c1.ra.rad, c1.dec.rad, unit=units.rad, frame='fk5')
        eclip =  a1.transform_to('geocentrictrueecliptic')
        #print('c1.ecliptic() = ',c1.ecliptic())
        #print('a1.ecliptic = ',eclip)
        el, b = c1.ecliptic()
        el = el.wrap(eclip.lon.rad * coord.radians)
        if c1 not in [north, south]:
            np.testing.assert_almost_equal(el.rad, eclip.lon.rad, decimal=4)
        np.testing.assert_almost_equal(b.rad, eclip.lat.rad, decimal=4)

        c2 = coord.CelestialCoord.from_ecliptic(el,b)
        a2 = astropy.coordinates.SkyCoord(el.rad, b.rad, unit=units.radian,
                                          frame='geocentrictrueecliptic')
        a2 = a2.transform_to('fk5')
        c2_ra = c2.ra.wrap(a2.ra.rad * coord.radians)
        np.testing.assert_almost_equal(c2_ra.rad, a2.ra.rad, decimal=3)
        np.testing.assert_almost_equal(c2.dec.rad, a2.dec.rad, decimal=3)

if __name__ == '__main__':
    test_angle()
    test_basic()
    test_distance()
    test_precess()
    test_galactic()
    test_ecliptic()
