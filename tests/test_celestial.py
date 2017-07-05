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
from numpy import pi, sin, cos, tan, arccos, sqrt
from coord import radians, degrees, hours, arcmin, arcsec

@timer
def test_init():
    """Basic tests of CelestialCoord construction.
    """
    pass


@timer
def test_invalid():
    """Check some invalid initializations of CelestialCoord
    """


@timer
def test_pickle():
    """Check picklability of CelestialCoords
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
    # First, let's test some distances that are easy to figure out without any spherical trig.
    eq1 = coord.CelestialCoord(0. * radians, 0. * radians)  # point on the equator
    eq2 = coord.CelestialCoord(1. * radians, 0. * radians)  # 1 radian along equator
    eq3 = coord.CelestialCoord(pi * radians, 0. * radians)  # antipode of eq1
    north_pole = coord.CelestialCoord(0. * radians, pi/2. * radians)  # north pole
    south_pole = coord.CelestialCoord(0. * radians, -pi/2. * radians) # south pole

    np.testing.assert_almost_equal(eq1.distanceTo(eq2).rad, 1., decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(eq1).rad, 1., decimal=12)
    np.testing.assert_almost_equal(eq1.distanceTo(eq3).rad, pi, decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(eq3).rad, pi-1., decimal=12)

    np.testing.assert_almost_equal(north_pole.distanceTo(south_pole).rad, pi, decimal=12)

    np.testing.assert_almost_equal(eq1.distanceTo(north_pole).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(north_pole).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal(eq3.distanceTo(north_pole).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal(eq1.distanceTo(south_pole).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(south_pole).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal(eq3.distanceTo(south_pole).rad, pi/2., decimal=12)

    # Some random point
    c1 = coord.CelestialCoord(0.234 * radians, 0.342 * radians)
    # Same meridian
    c2 = coord.CelestialCoord(0.234 * radians, -1.093 * radians)
    # Antipode
    c3 = coord.CelestialCoord((pi + 0.234) * radians, -0.342 * radians)
    # Different point on opposide meridian
    c4 = coord.CelestialCoord((pi + 0.234) * radians, 0.832 * radians)

    for c, d in zip( (c1,c2,c3,c4), (0., 1.435, pi, pi-1.174) ):
        np.testing.assert_almost_equal(c1.distanceTo(c).rad, d, decimal=12)

    # Now some that require spherical trig calculations.
    # Importantly, this uses the more straightforward spherical trig formula, the cosine rule.
    # The CelestialCoord class uses a different formula that is more stable for very small
    # distances, which are typical in the correlation function calculation.
    # Some other random point:
    c5 = coord.CelestialCoord(1.832 * radians, -0.723 * radians)
    # The standard formula is:
    # cos(d) = sin(dec1) sin(dec2) + cos(dec1) cos(dec2) cos(delta ra)
    d = arccos(sin(c1.dec) * sin(c5.dec) + cos(c1.dec) * cos(c5.dec) * cos(c1.ra-c5.ra))
    np.testing.assert_almost_equal(c1.distanceTo(c5).rad, d, decimal=12)

    # Tiny displacements should have dsq = (dra^2 cos^2 dec) + (ddec^2)
    c6 = coord.CelestialCoord(c1.ra + 1.7e-9 * radians, c1.dec)
    c7 = coord.CelestialCoord(c1.ra, c1.dec + 1.9e-9 * radians)
    c8 = coord.CelestialCoord(c1.ra + 2.3e-9 * radians, c1.dec + 1.2e-9 * radians)

    # Note that the standard formula gets these wrong.  d comes back as 0.
    d = arccos(sin(c1.dec) * sin(c6.dec) + cos(c1.dec) * cos(c6.dec) * cos(c1.ra-c6.ra))
    print('d(c6) = ',1.7e-9 * cos(0.342), c1.distanceTo(c6), d)
    d = arccos(sin(c1.dec) * sin(c7.dec) + cos(c1.dec) * cos(c7.dec) * cos(c1.ra-c7.ra))
    print('d(c7) = ',1.9e-9, c1.distanceTo(c7), d)
    d = arccos(sin(c1.dec) * sin(c8.dec) + cos(c1.dec) * cos(c8.dec) * cos(c1.ra-c8.ra))
    true_d = sqrt( (2.3e-9 * cos(0.342))**2 + 1.2e-9**2)
    print('d(c7) = ',true_d, c1.distanceTo(c8), d)
    np.testing.assert_allclose(c1.distanceTo(c6).rad, 1.7e-9 * c1.dec.cos(), rtol=1.e-7)
    np.testing.assert_allclose(c1.distanceTo(c7).rad, 1.9e-9, rtol=1.e-7)
    np.testing.assert_allclose(c1.distanceTo(c8).rad, true_d, rtol=1.e-7)

    # Near antipodes, the formula we usually use becomes somewhat inaccurate.
    # Check a variety of antipodes.
    eq1 = coord.CelestialCoord(0. * radians, 0. * radians)
    eq2 = coord.CelestialCoord(1. * radians, 0. * radians)
    eq3 = coord.CelestialCoord(pi * radians, 0. * radians)
    north_pole = coord.CelestialCoord(0. * radians, pi/2. * radians)
    south_pole = coord.CelestialCoord(0. * radians, -pi/2. * radians)
    for c in [c1, c2, c3, c4, c5, c6, c7, c8, eq1, eq2, eq3, north_pole, south_pole]:
        antipode = coord.CelestialCoord(c.ra + pi * radians, -c.dec)
        np.testing.assert_almost_equal(c.distanceTo(antipode).rad, pi, decimal=12)
        np.testing.assert_almost_equal(antipode.distanceTo(c).rad, pi, decimal=12)

    # Also some near, but not quite antipodes
    # Note: this range crosses the point where the formula in distanceTo changes.
    for delta in range(500):
        eq4 = coord.CelestialCoord(delta * arcmin, 0 * radians)
        np.testing.assert_allclose(pi - eq3.distanceTo(eq4).rad, eq4.ra.rad, rtol=1.e-7)


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
    cB = coord.CelestialCoord((-0.193 + 1.7e-8) * radians,
                               (0.882 + 1.2e-8) * radians)
    cC = coord.CelestialCoord((-0.193 - 2.4e-8) * radians,
                               (0.882 + 3.1e-8) * radians)

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

    # The shoelace formula gives the area of a triangle given coordinates:
    # A = 1/2 abs( (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1) )
    area = 0.5 * abs( (uB.rad-uA.rad) * (vC.rad-vA.rad) - (uC.rad-uA.rad) * (vB.rad-vA.rad) )
    print('lambert area = ',area,E)
    np.testing.assert_allclose(area, E, err_msg="lambert didn't preserve area")

    # Check that project_rad does the same thing
    uA2, vA2 = center.project_rad(cA.ra.rad, cA.dec.rad, projection='lambert')
    np.testing.assert_allclose([uA2,vA2], [uA.rad,vA.rad], rtol=1.e-8,
                               err_msg="project_rad not equivalent")

    # Check the deprojection
    cA2 = center.deproject(uA, vA, projection='lambert')
    np.testing.assert_allclose(cA2.rad, cA.rad, err_msg="deproject didn't return to orig")
    cA3 = center.deproject_rad(uA.rad, vA.rad, projection='lambert')
    np.testing.assert_allclose(cA3, cA.rad, err_msg="deproject_rad not equivalent")

    # The angles are not preserved
    a = sqrt( (uB.rad-uC.rad)**2 + (vB.rad-vC.rad)**2 )
    b = sqrt( (uC.rad-uA.rad)**2 + (vC.rad-vA.rad)**2 )
    c = sqrt( (uA.rad-uB.rad)**2 + (vA.rad-vB.rad)**2 )
    cosA = ((uB.rad-uA.rad)*(uC.rad-uA.rad) + (vB.rad-vA.rad)*(vC.rad-vA.rad)) / (b*c)
    cosB = ((uC.rad-uB.rad)*(uA.rad-uB.rad) + (vC.rad-vB.rad)*(vA.rad-vB.rad)) / (c*a)
    cosC = ((uA.rad-uC.rad)*(uB.rad-uC.rad) + (vA.rad-vC.rad)*(vB.rad-vC.rad)) / (a*b)

    print('lambert cosA = ',cosA,cos(A))
    print('lambert cosB = ',cosB,cos(B))
    print('lambert cosC = ',cosC,cos(C))

    # The deproject jacobian should tell us how the area changes
    dudx, dudy, dvdx, dvdy = center.jac_deproject(uA, vA, 'lambert').ravel()
    jac_area = abs(dudx*dvdy - dudy*dvdx)
    np.testing.assert_allclose(jac_area, E/area, err_msg='jac_deproject gave wrong area')

    dudx, dudy, dvdx, dvdy = center.jac_deproject_rad(uA.rad, vA.rad, 'lambert').ravel()
    np.testing.assert_allclose(jac_area, abs(dudx*dvdy - dudy*dvdx),
                               err_msg='jac_deproject_rad not equivalent')

    # center projects to 0,0 with unit area
    u, v = center.project(center, 'lambert')
    np.testing.assert_allclose([u.rad, v.rad], 0., err_msg='center did not project to (0,0)')
    c2 = center.deproject(u,v, 'lambert')
    np.testing.assert_allclose(c2.rad, center.rad, err_msg='(0,0) did not deproject to center')
    np.testing.assert_allclose(np.linalg.det(center.jac_deproject(u, v, 'lambert')), 1.,
                               err_msg='determinant of jac_deproject matrix != 1 at center')


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
    test_init()
    test_invalid()
    test_pickle()
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
