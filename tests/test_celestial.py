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
