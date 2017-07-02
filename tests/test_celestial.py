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

from __future__ import print_function
import numpy as np
import math
import coord

from helper_util import *

# We'll use these a lot, so just import them.
from math import pi, sin, cos, tan, acos, sqrt
from coord import radians, degrees, hours, arcmin, arcsec

@timer
def test_basic():
    """Basic tests of CelestialCoord construction. etc.
    """
    pass


@timer
def test_pickle():
    """Check picklability of CelestialCoords
    """


@timer
def test_invalid():
    """Check some invalid initializations of CelestialCoord
    """


@timer
def test_eq():
    """Check that equal coords are equal, but unequal ones are not.
    """
    pass


@timer
def test_distance():
    """Test calculations of distances on the sphere.
    """
    pass


@timer
def test_xyz():
    """Test get_xyz and from_xyz functions
    """


@timer
def test_angleBetween():
    """Test calculations of angles between positions on the sphere.
    """
    pass


@timer
def test_gnomonic_projection():
    """Test the gnomonic projection.
    """
    pass


@timer
def test_stereographic_projection():
    """Test the stereographic projection.
    """
    pass


@timer
def test_lambert_projection():
    """Test the lambert projection.
    """
    # Test that a small triangle has the correct properties.
    # Lambert projections preserve area, but not angles.
    center = coord.CelestialCoord(0.234 * radians, 0.342 * radians)
    cA = coord.CelestialCoord(-0.193 * radians, 0.882 * radians)
    cB = coord.CelestialCoord((-0.193 + 1.7e-6) * radians,
                               (0.882 + 1.2e-6) * radians)
    cC = coord.CelestialCoord((-0.193 - 2.4e-6) * radians,
                               (0.882 + 3.1e-6) * radians)

    a = cB.distanceTo(cC).rad
    b = cC.distanceTo(cA).rad
    c = cA.distanceTo(cB).rad
    A = cA.angleBetween(cB,cC).rad
    B = cB.angleBetween(cC,cA).rad
    C = cC.angleBetween(cA,cB).rad
    E = cA.area(cB,cC)

    uA, vA = center.project(cA, projection='lambert')
    uB, vB = center.project(cB, projection='lambert')
    uC, vC = center.project(cC, projection='lambert')
    uA = uA / arcsec  # Easier to just deal with these in arcsec
    vA = vA / arcsec
    uB = uB / arcsec
    vB = vB / arcsec
    uC = uC / arcsec
    vC = vC / arcsec

    # The shoelace formula gives the area of a triangle given coordinates:
    # A = 1/2 abs( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) )
    area = 0.5 * abs( (uB-uA) * (vC-vA) - (uC-uA) * (vB-vA) )
    area *= (arcsec / radians)**2
    print('lambert area = ',area,E)
    np.testing.assert_almost_equal(area / E, 1, decimal=5)

    # Check that project_rad does the same thing
    uA2, vA2 = center.project_rad(cA.ra.rad, cA.dec.rad, projection='lambert')
    np.testing.assert_array_almost_equal(uA, uA2, decimal=10)
    np.testing.assert_array_almost_equal(vA, vA2, decimal=10)

    # Check the deprojection
    cA2 = center.deproject(uA*arcsec, vA*arcsec, projection='lambert')
    np.testing.assert_almost_equal(cA.ra.rad, cA2.ra.rad, decimal=12)
    np.testing.assert_almost_equal(cA.dec.rad, cA2.dec.rad, decimal=12)
    cA3 = center.deproject_rad(uA, vA, projection='lambert')
    np.testing.assert_array_almost_equal([cA.ra.rad, cA.dec.rad], cA3, decimal=12)

    # The angles are not preserved
    a = sqrt( (uB-uC)**2 + (vB-vC)**2 )
    b = sqrt( (uC-uA)**2 + (vC-vA)**2 )
    c = sqrt( (uA-uB)**2 + (vA-vB)**2 )
    cosA = ((uB-uA)*(uC-uA) + (vB-vA)*(vC-vA)) / (b*c)
    cosB = ((uC-uB)*(uA-uB) + (vC-vB)*(vA-vB)) / (c*a)
    cosC = ((uA-uC)*(uB-uC) + (vA-vC)*(vB-vC)) / (a*b)

    print('lambert cosA = ',cosA,cos(A))
    print('lambert cosB = ',cosB,cos(B))
    print('lambert cosC = ',cosC,cos(C))

    # The deproject jacobian should tell us how the area changes
    dudx, dudy, dvdx, dvdy = center.jac_deproject(uA*arcsec, vA*arcsec, 'lambert').ravel()
    jac_area = abs(dudx*dvdy - dudy*dvdx)
    np.testing.assert_almost_equal(jac_area, E/area, decimal=5)

    dudx, dudy, dvdx, dvdy = center.jac_deproject_arcsec(uA, vA, 'lambert').ravel()
    np.testing.assert_almost_equal(jac_area, abs(dudx*dvdy - dudy*dvdx))

    # center projects to 0,0 with unit area
    u, v = center.project(center, 'lambert')
    np.testing.assert_almost_equal(u.rad, 0., decimal=12)
    np.testing.assert_almost_equal(v.rad, 0., decimal=12)
    c2 = center.deproject(u,v, 'lambert')
    np.testing.assert_almost_equal(c2.ra.rad, center.ra.rad, decimal=12)
    np.testing.assert_almost_equal(c2.dec.rad, center.dec.rad, decimal=12)
    np.testing.assert_almost_equal(np.linalg.det(center.jac_deproject(u, v, 'lambert')), 1.)


@timer
def test_postel_projection():
    """Test the postel projection.
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


@timer
def test_ecliptic_date():
    """Test the date option of the ecliptic and from_ecliptic functions.
    """
    pass


if __name__ == '__main__':
    test_basic()
    test_eq()
    test_distance()
    test_angleBetween()
    test_gnomonic_projection()
    test_stereographic_projection()
    test_lambert_projection()
    test_postel_projection()
    test_precess()
    test_galactic()
    test_ecliptic()
    test_ecliptic_date()
