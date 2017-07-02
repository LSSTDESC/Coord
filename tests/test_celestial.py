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
    c1 = coord.CelestialCoord(0. * radians, 0. * radians)
    np.testing.assert_almost_equal(c1.ra.rad, 0., decimal=12)
    np.testing.assert_almost_equal(c1.dec.rad, 0., decimal=12)

    # These next three are the same place
    c2 = coord.CelestialCoord(11. * hours, -37. * degrees)
    np.testing.assert_almost_equal(c2.ra / hours, 11., decimal=12)
    np.testing.assert_almost_equal(c2.dec / degrees, -37., decimal=12)

    c3 = coord.CelestialCoord(35. * hours, -37. * degrees)
    np.testing.assert_almost_equal(c3.ra / hours, 35., decimal=12)
    np.testing.assert_almost_equal(c3.dec / degrees, -37., decimal=12)

    c4 = coord.CelestialCoord(-13. * hours, -37. * degrees)
    np.testing.assert_almost_equal(c4.ra / hours, -13., decimal=12)
    np.testing.assert_almost_equal(c4.dec / degrees, -37., decimal=12)

    # We'll test distance later, but for check that these last 3 have distances = 0.
    np.testing.assert_almost_equal(c2.distanceTo(c3).rad, 0., decimal=12)
    np.testing.assert_almost_equal(c2.distanceTo(c4).rad, 0., decimal=12)

    # Check picklability
    do_pickle(c1)
    do_pickle(c2)
    do_pickle(c3)
    do_pickle(c4)


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
    c1 = coord.CelestialCoord(11. * hours, -37. * degrees)
    c2 = coord.CelestialCoord(165. * degrees, -37. * degrees)
    assert c2 == c1

    # These should all test as unequal.  Note some non-Angles in the list.
    diff_list = [ c1,
                  coord.CelestialCoord(0 * radians, 0 * radians),
                  coord.CelestialCoord(13. * hours, 37. * degrees),
                  coord.CelestialCoord(11. * hours, 37. * degrees),
                  pi, coord.CelestialCoord, None ]
    all_obj_diff(diff_list)


@timer
def test_distance():
    """Test calculations of distances on the sphere.
    """
    # First, let's test some distances that are easy to figure out without any spherical trig.
    eq1 = coord.CelestialCoord(0 * radians, 0 * radians)  # point on the equator
    eq2 = coord.CelestialCoord(1 * radians, 0 * radians)  # 1 radian along equator
    eq3 = coord.CelestialCoord(pi * radians, 0 * radians)  # antipode of eq1
    north_pole = coord.CelestialCoord(0 * radians, pi/2 * radians)  # north pole
    south_pole = coord.CelestialCoord(0 * radians, -pi/2 * radians) # south pole

    np.testing.assert_almost_equal(eq1.distanceTo(eq2).rad, 1, decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(eq1).rad, 1, decimal=12)
    np.testing.assert_almost_equal(eq1.distanceTo(eq3).rad, pi, decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(eq3).rad, pi-1, decimal=12)

    np.testing.assert_almost_equal(north_pole.distanceTo(south_pole).rad, pi, decimal=12)

    np.testing.assert_almost_equal(eq1.distanceTo(north_pole).rad, pi/2, decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(north_pole).rad, pi/2, decimal=12)
    np.testing.assert_almost_equal(eq3.distanceTo(north_pole).rad, pi/2, decimal=12)
    np.testing.assert_almost_equal(eq1.distanceTo(south_pole).rad, pi/2, decimal=12)
    np.testing.assert_almost_equal(eq2.distanceTo(south_pole).rad, pi/2, decimal=12)
    np.testing.assert_almost_equal(eq3.distanceTo(south_pole).rad, pi/2, decimal=12)

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

    # Note that the standard formula gets these wrong.  d is 0.0
    d = arccos(sin(c1.dec) * sin(c6.dec) + cos(c1.dec) * cos(c6.dec) * cos(c1.ra-c6.ra))
    print('d(c6) = ',1.7e-9 * cos(c1.dec), c1.distanceTo(c6), d)
    d = arccos(sin(c1.dec) * sin(c7.dec) + cos(c1.dec) * cos(c7.dec) * cos(c1.ra-c7.ra))
    print('d(c7) = ',1.9e-9, c1.distanceTo(c7), d)
    d = arccos(sin(c1.dec) * sin(c8.dec) + cos(c1.dec) * cos(c8.dec) * cos(c1.ra-c8.ra))
    true_d = sqrt( (2.3e-9 * cos(c1.dec))**2 + 1.2e-9**2)
    print('d(c7) = ',true_d, c1.distanceTo(c8), d)
    np.testing.assert_allclose(c1.distanceTo(c6).rad, 1.7e-9 * cos(c1.dec), rtol=1.e-7)
    np.testing.assert_allclose(c1.distanceTo(c7).rad, 1.9e-9, rtol=1.e-7)
    np.testing.assert_allclose(c1.distanceTo(c8).rad, true_d, rtol=1.e-7)

    # Near antipodes, the formula we usually use becomes somewhat inaccurate.
    # Check a variety of antipodes.
    for c in [c1, c2, c3, c4, c5, c6, c7, c8, eq1, eq2, eq3, north_pole, south_pole]:
        antipode = coord.CelestialCoord(c.ra + pi*radians, -c.dec)
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
    eq1 = coord.CelestialCoord(0. * radians, 0. * radians)  # point on the equator
    eq2 = coord.CelestialCoord(1. * radians, 0. * radians)  # 1 radian along equator
    eq3 = coord.CelestialCoord(pi * radians, 0. * radians) # antipode of eq1
    north_pole = coord.CelestialCoord(0. * radians, pi/2. * radians)  # north pole
    south_pole = coord.CelestialCoord(0. * radians, -pi/2. * radians) # south pole
    c1 = coord.CelestialCoord(0.234 * radians, 0.342 * radians)
    c2 = coord.CelestialCoord(0.234 * radians, -1.093 * radians)
    c3 = coord.CelestialCoord((pi + 0.234) * radians, -0.342 * radians)
    c4 = coord.CelestialCoord((pi + 0.234) * radians, 0.832 * radians)
    c5 = coord.CelestialCoord(1.832 * radians, -0.723 * radians)

    for c in [c1, c2, c3, c4, c5, eq1, eq2, eq3, north_pole, south_pole]:
        x, y, z = c.get_xyz()
        np.testing.assert_almost_equal(x, c.dec.cos() * c.ra.cos(), decimal=12)
        np.testing.assert_almost_equal(y, c.dec.cos() * c.ra.sin(), decimal=12)
        np.testing.assert_almost_equal(z, c.dec.sin(), decimal=12)

        cc1 = coord.CelestialCoord.from_xyz(x,y,z)
        np.testing.assert_almost_equal(cc1.ra.rad, c.ra.rad, decimal=12)
        np.testing.assert_almost_equal(cc1.dec.rad, c.dec.rad, decimal=12)

        # Works in x,y,z are scaled arbitrarily
        cc2 = coord.CelestialCoord.from_xyz(x*17,y*17,z*17)
        np.testing.assert_almost_equal(cc2.ra.rad, c.ra.rad, decimal=12)
        np.testing.assert_almost_equal(cc2.dec.rad, c.dec.rad, decimal=12)

        cc3 = coord.CelestialCoord.from_xyz(x*1.e-9,y*1.e-9,z*1.e-9)
        np.testing.assert_almost_equal(cc3.ra.rad, c.ra.rad, decimal=12)
        np.testing.assert_almost_equal(cc3.dec.rad, c.dec.rad, decimal=12)

    # constructing from x,y,z = 0,0,0 is undefined.
    np.testing.assert_raises(ValueError, coord.CelestialCoord.from_xyz, 0., 0., 0.)


@timer
def test_xyz():
    """Test get_xyz and from_xyz functions
    """


@timer
def test_angleBetween():
    """Test calculations of angles between positions on the sphere.
    """
    # Again, let's start with some answers we can get by inspection.
    eq1 = coord.CelestialCoord(0. * radians, 0. * radians)  # point on the equator
    eq2 = coord.CelestialCoord(1. * radians, 0. * radians)  # 1 radian along equator
    eq3 = coord.CelestialCoord(pi * radians, 0. * radians) # antipode of eq1
    north_pole = coord.CelestialCoord(0. * radians, pi/2. * radians)  # north pole
    south_pole = coord.CelestialCoord(0. * radians, -pi/2. * radians) # south pole

    np.testing.assert_almost_equal(north_pole.angleBetween(eq1,eq2).rad, -1., decimal=12)
    np.testing.assert_almost_equal(north_pole.angleBetween(eq2,eq1).rad, 1., decimal=12)
    np.testing.assert_almost_equal(north_pole.angleBetween(eq2,eq3).rad, 1.-pi, decimal=12)
    np.testing.assert_almost_equal(north_pole.angleBetween(eq3,eq2).rad, pi-1., decimal=12)
    np.testing.assert_almost_equal(south_pole.angleBetween(eq1,eq2).rad, 1., decimal=12)
    np.testing.assert_almost_equal(south_pole.angleBetween(eq2,eq1).rad, -1., decimal=12)
    np.testing.assert_almost_equal(south_pole.angleBetween(eq2,eq3).rad, pi-1., decimal=12)
    np.testing.assert_almost_equal(south_pole.angleBetween(eq3,eq2).rad, 1.-pi, decimal=12)

    np.testing.assert_almost_equal(eq1.angleBetween(north_pole,eq2).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal(eq2.angleBetween(north_pole,eq1).rad, -pi/2., decimal=12)

    np.testing.assert_almost_equal(north_pole.area(eq1,eq2), 1., decimal=12)
    np.testing.assert_almost_equal(north_pole.area(eq2,eq1), 1., decimal=12)
    np.testing.assert_almost_equal(south_pole.area(eq1,eq2), 1., decimal=12)
    np.testing.assert_almost_equal(south_pole.area(eq2,eq1), 1., decimal=12)

    # For arbitrary points, we can check that the spherical triangle satisfies
    # the spherical trig laws.
    cA = coord.CelestialCoord(0.234 * radians, 0.342 * radians)
    cB = coord.CelestialCoord(-0.193 * radians, 0.882 * radians)
    cC = coord.CelestialCoord(0.721 * radians, -0.561 * radians)

    a = cB.distanceTo(cC)
    b = cC.distanceTo(cA)
    c = cA.distanceTo(cB)
    A = cA.angleBetween(cB,cC)
    B = cB.angleBetween(cC,cA)
    C = cC.angleBetween(cA,cB)
    E = abs(A.rad)+abs(B.rad)+abs(C.rad)-pi
    s = (a+b+c)/2.

    # Law of cosines:
    np.testing.assert_almost_equal(c.cos(), a.cos()*b.cos() + a.sin()*b.sin()*C.cos(), decimal=12)
    np.testing.assert_almost_equal(a.cos(), b.cos()*c.cos() + b.sin()*c.sin()*A.cos(), decimal=12)
    np.testing.assert_almost_equal(b.cos(), c.cos()*a.cos() + c.sin()*a.sin()*B.cos(), decimal=12)

    # Law of sines:
    np.testing.assert_almost_equal(A.sin() * b.sin(), B.sin() * a.sin(), decimal=12)
    np.testing.assert_almost_equal(B.sin() * c.sin(), C.sin() * b.sin(), decimal=12)
    np.testing.assert_almost_equal(C.sin() * a.sin(), A.sin() * c.sin(), decimal=12)

    # Alternate law of cosines:
    np.testing.assert_almost_equal(C.cos(), -A.cos()*B.cos() + A.sin()*B.sin()*c.cos(), decimal=12)
    np.testing.assert_almost_equal(A.cos(), -B.cos()*C.cos() + B.sin()*C.sin()*a.cos(), decimal=12)
    np.testing.assert_almost_equal(B.cos(), -C.cos()*A.cos() + C.sin()*A.sin()*b.cos(), decimal=12)

    # Spherical excess:
    np.testing.assert_almost_equal(cA.area(cB,cC), E, decimal=12)
    np.testing.assert_almost_equal(cA.area(cC,cB), E, decimal=12)
    np.testing.assert_almost_equal(cB.area(cA,cC), E, decimal=12)
    np.testing.assert_almost_equal(cB.area(cC,cA), E, decimal=12)
    np.testing.assert_almost_equal(cC.area(cB,cA), E, decimal=12)
    np.testing.assert_almost_equal(cC.area(cA,cB), E, decimal=12)

    # L'Huilier's formula for spherical excess:
    np.testing.assert_almost_equal(
            tan(E/4)**2,
            (s/2).tan() * ((s-a)/2).tan() * ((s-b)/2).tan() * ((s-c)/2).tan(), decimal=12)


@timer
def test_gnomonic_projection():
    """Test the gnomonic projection.
    """
    # Test that a small triangle has the correct properties.
    # Gnomonic projections turn great circles into straight lines.
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

    uA, vA = center.project(cA, projection='gnomonic')
    uB, vB = center.project(cB, projection='gnomonic')
    uC, vC = center.project(cC, projection='gnomonic')
    uA = uA / arcsec
    vA = vA / arcsec
    uB = uB / arcsec
    vB = vB / arcsec
    uC = uC / arcsec
    vC = vC / arcsec

    # Check that project_rad does the same thing
    uA2, vA2 = center.project_rad(cA.ra.rad, cA.dec.rad, projection='gnomonic')
    np.testing.assert_array_almost_equal(uA, uA2, decimal=12)
    np.testing.assert_array_almost_equal(vA, vA2, decimal=12)

    # Check the deprojection
    cA2 = center.deproject(uA*arcsec, vA*arcsec, projection='gnomonic')
    np.testing.assert_almost_equal(cA.ra.rad, cA2.ra.rad, decimal=12)
    np.testing.assert_almost_equal(cA.dec.rad, cA2.dec.rad, decimal=12)
    cA3 = center.deproject_rad(uA, vA, projection='gnomonic')
    np.testing.assert_array_almost_equal([cA.ra.rad, cA.dec.rad], cA3, decimal=12)

    # The angles are not preserved
    a = sqrt( (uB-uC)**2 + (vB-vC)**2 )
    b = sqrt( (uC-uA)**2 + (vC-vA)**2 )
    c = sqrt( (uA-uB)**2 + (vA-vB)**2 )
    cosA = ((uB-uA)*(uC-uA) + (vB-vA)*(vC-vA)) / (b*c)
    cosB = ((uC-uB)*(uA-uB) + (vC-vB)*(vA-vB)) / (c*a)
    cosC = ((uA-uC)*(uB-uC) + (vA-vC)*(vB-vC)) / (a*b)

    print('gnomonic cosA = ',cosA,cos(A))
    print('gnomonic cosB = ',cosB,cos(B))
    print('gnomonic cosC = ',cosC,cos(C))

    # The area is not preserved
    area = 0.5 * abs( (uB-uA)*(vC-vA) - (uC-uA)*(vB-vA) )
    area *= (arcsec / radians)**2
    print('gnomonic area = ',area,E)

    # The deproject jacobian should tell us how the area changes
    dudx, dudy, dvdx, dvdy = center.jac_deproject(uA*arcsec, vA*arcsec, 'gnomonic').ravel()
    jac_area = abs(dudx*dvdy - dudy*dvdx)
    np.testing.assert_almost_equal(jac_area, E/area, decimal=5)

    dudx, dudy, dvdx, dvdy = center.jac_deproject_arcsec(uA, vA, 'gnomonic').ravel()
    np.testing.assert_almost_equal(jac_area, abs(dudx*dvdy - dudy*dvdx))

    # center projects to 0,0 with unit area
    u, v = center.project(center, 'gnomonic')
    np.testing.assert_almost_equal(u.rad, 0., decimal=12)
    np.testing.assert_almost_equal(v.rad, 0., decimal=12)
    c2 = center.deproject(u,v, 'gnomonic')
    np.testing.assert_almost_equal(c2.ra.rad, center.ra.rad, decimal=12)
    np.testing.assert_almost_equal(c2.dec.rad, center.dec.rad, decimal=12)
    np.testing.assert_almost_equal(np.linalg.det(center.jac_deproject(u, v, 'gnomonic')), 1.)

    # Finally, check the claim that great circles turn into straight lines.
    # Meridians are great circles.  Also symmetric points across the equator.
    center = coord.CelestialCoord(0.234 * radians, 0.342 * radians)
    gc_list = [ [ (0.1, -0.5), (0.1, -0.3), (0.1, 0.1), (0.1, 0.4) ],
                [ (0.3, -0.5), (0.3, -0.3), (0.3, 0.1), (0.3, 0.4) ],
                [ (0.5, -0.5), (0.5, -0.3), (0.5, 0.1), (0.5, 0.4) ],
                [ (0.7, -0.5), (0.7, -0.3), (0.7, 0.1), (0.7, 0.4) ],
                [ (-0.1, -0.4), (0.1, 0.0), (0.3, 0.4) ],
                [ (0.3, -0.2), (0.4, 0.0), (0.5, 0.2) ],
                [ (0.3, -0.1), (0.0, 0.0), (-0.3, 0.1) ] ]
    for gc in gc_list:
        print('gc = ',gc)
        # The projected points for a single great circle
        proj = [center.project(coord.CelestialCoord(c[0] * radians, c[1] * radians)) for c in gc]
        # Calculate the slope from the first and last points
        u0, v0 = proj[0]
        un, vn = proj[-1]
        m0 = (vn.rad - v0.rad) / (un.rad - u0.rad)
        print('m0 = ',m0)
        for (u,v) in proj[1:-1]:
            m1 = (v.rad - v0.rad) / (u.rad - u0.rad)
            print('u,v = ',u,v,', m1 = ',m1)
            np.testing.assert_almost_equal(m1, m0, decimal=12)

    # gnomonic is the default.  Make sure that works right.
    uA3, vA3 = center.project(cA)
    uB3, vB3 = center.project(cB)
    uC3, vC3 = center.project(cC)
    np.testing.assert_array_almost_equal(uA, uA3 / arcsec, decimal=12)
    np.testing.assert_array_almost_equal(vA, vA3 / arcsec, decimal=12)
    np.testing.assert_array_almost_equal(uB, uB3 / arcsec, decimal=12)
    np.testing.assert_array_almost_equal(vB, vB3 / arcsec, decimal=12)
    np.testing.assert_array_almost_equal(uC, uC3 / arcsec, decimal=12)
    np.testing.assert_array_almost_equal(vC, vC3 / arcsec, decimal=12)

    # And make sure that invalid projection strings raise exceptions
    np.testing.assert_raises(ValueError, center.project, cA, 'TAN')
    np.testing.assert_raises(ValueError, center.project, cA, projection=3)
    np.testing.assert_raises(ValueError, center.project, 3, 4)
    np.testing.assert_raises(ValueError, center.project, u, v)
    np.testing.assert_raises(ValueError, center.project_rad, 3, 4, 'TAN')
    np.testing.assert_raises(ValueError, center.project_rad, 3, 4, projection=3)
    np.testing.assert_raises(TypeError, center.project_rad, cA)
    np.testing.assert_raises(ValueError, center.deproject, u, v, 'TAN')
    np.testing.assert_raises(ValueError, center.deproject, u, v, projection=3)
    np.testing.assert_raises(TypeError, center.deproject, 3, 4)
    np.testing.assert_raises(TypeError, center.deproject, cA)
    np.testing.assert_raises(ValueError, center.deproject_rad, 3, 4, 'TAN')
    np.testing.assert_raises(ValueError, center.deproject_rad, 3, 4, projection=3)
    np.testing.assert_raises(TypeError, center.deproject_rad, u, v)
    np.testing.assert_raises(TypeError, center.deproject_rad, cA)
    np.testing.assert_raises(ValueError, center.jac_deproject, u, v, 'TAN')
    np.testing.assert_raises(ValueError, center.jac_deproject, u, v, projection=3)
    np.testing.assert_raises(AttributeError, center.jac_deproject, 3, 4)
    np.testing.assert_raises(TypeError, center.jac_deproject, cA)
    np.testing.assert_raises(ValueError, center.jac_deproject_arcsec, 3, 4, 'TAN')
    np.testing.assert_raises(ValueError, center.jac_deproject_arcsec, 3, 4, projection=3)
    np.testing.assert_raises(TypeError, center.jac_deproject_arcsec, u, v)
    np.testing.assert_raises(TypeError, center.jac_deproject_arcsec, cA)


@timer
def test_stereographic_projection():
    """Test the stereographic projection.
    """
    # Test that a small triangle has the correct properties.
    # Stereographic projections preserve angles, but not area.
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

    uA, vA = center.project(cA, projection='stereographic')
    uB, vB = center.project(cB, projection='stereographic')
    uC, vC = center.project(cC, projection='stereographic')
    uA = uA / arcsec
    vA = vA / arcsec
    uB = uB / arcsec
    vB = vB / arcsec
    uC = uC / arcsec
    vC = vC / arcsec

    # The easiest way to compute the angles is from the dot products:
    # a.b = ab cos(C)
    a = sqrt( (uB-uC)**2 + (vB-vC)**2 )
    b = sqrt( (uC-uA)**2 + (vC-vA)**2 )
    c = sqrt( (uA-uB)**2 + (vA-vB)**2 )
    cosA = ((uB-uA)*(uC-uA) + (vB-vA)*(vC-vA)) / (b*c)
    cosB = ((uC-uB)*(uA-uB) + (vC-vB)*(vA-vB)) / (c*a)
    cosC = ((uA-uC)*(uB-uC) + (vA-vC)*(vB-vC)) / (a*b)

    print('stereographic cosA = ',cosA,cos(A))
    print('stereographic cosB = ',cosB,cos(B))
    print('stereographic cosC = ',cosC,cos(C))
    np.testing.assert_almost_equal(cosA, cos(A), decimal=5)
    np.testing.assert_almost_equal(cosB, cos(B), decimal=5)
    np.testing.assert_almost_equal(cosC, cos(C), decimal=5)

    # Check that project_rad does the same thing
    uA2, vA2 = center.project_rad(cA.ra.rad, cA.dec.rad, projection='stereographic')
    np.testing.assert_array_almost_equal(uA, uA2, decimal=10)
    np.testing.assert_array_almost_equal(vA, vA2, decimal=10)

    # Check the deprojection
    cA2 = center.deproject(uA*arcsec, vA*arcsec, projection='stereographic')
    np.testing.assert_almost_equal(cA.ra.rad, cA2.ra.rad, decimal=12)
    np.testing.assert_almost_equal(cA.dec.rad, cA2.dec.rad, decimal=12)
    cA3 = center.deproject_rad(uA, vA, projection='stereographic')
    np.testing.assert_array_almost_equal([cA.ra.rad, cA.dec.rad], cA3, decimal=12)

    # The area is not preserved
    area = 0.5 * abs( (uB-uA) * (vC-vA) - (uC-uA) * (vB-vA) )
    area *= (arcsec / radians)**2
    print('stereographic area = ',area,E)

    # The deproject jacobian should tell us how the area changes
    dudx, dudy, dvdx, dvdy = center.jac_deproject(uA*arcsec, vA*arcsec, 'stereographic').ravel()
    jac_area = abs(dudx*dvdy - dudy*dvdx)
    np.testing.assert_almost_equal(jac_area, E/area, decimal=5)

    dudx, dudy, dvdx, dvdy = center.jac_deproject_arcsec(uA, vA, 'stereographic').ravel()
    np.testing.assert_almost_equal(jac_area, abs(dudx*dvdy - dudy*dvdx))

    # center projects to 0,0 with unit area
    u, v = center.project(center, 'stereographic')
    np.testing.assert_almost_equal(u.rad, 0., decimal=12)
    np.testing.assert_almost_equal(v.rad, 0., decimal=12)
    c2 = center.deproject(u,v, 'stereographic')
    np.testing.assert_almost_equal(c2.ra.rad, center.ra.rad, decimal=12)
    np.testing.assert_almost_equal(c2.dec.rad, center.dec.rad, decimal=12)
    np.testing.assert_almost_equal(np.linalg.det(center.jac_deproject(u, v, 'stereographic')), 1.)


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
    np.testing.assert_allclose(area, E, rtol=1.e-8, err_msg="lambert didn't preserve area")

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

    print('lambert cosA = ',cosA,cos(A)) # 0.385002890444 0.371329625008
    print('lambert cosB = ',cosB,cos(B)) # 0.13803694316 0.102895820012
    print('lambert cosC = ',cosC,cos(C)) # 0.860935750473 0.885364487706
    assert cosA > cos(A) + 0.01
    assert cosB > cos(B) + 0.03
    assert cosC < cos(C) - 0.02

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
    # Test that a small triangle has the correct properties.
    # Postel projections preserve distance from the center.
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

    uA, vA = center.project(cA, projection='postel')
    uB, vB = center.project(cB, projection='postel')
    uC, vC = center.project(cC, projection='postel')
    uA = uA / arcsec
    vA = vA / arcsec
    uB = uB / arcsec
    vB = vB / arcsec
    uC = uC / arcsec
    vC = vC / arcsec

    dA = sqrt( uA**2 + vA**2 )
    dB = sqrt( uB**2 + vB**2 )
    dC = sqrt( uC**2 + vC**2 )
    print('postel dA = ',dA,center.distanceTo(cA))
    print('postel dB = ',dB,center.distanceTo(cB))
    print('postel dC = ',dC,center.distanceTo(cC))
    np.testing.assert_almost_equal(dA, center.distanceTo(cA) / arcsec, decimal=10)
    np.testing.assert_almost_equal(dB, center.distanceTo(cB) / arcsec, decimal=10)
    np.testing.assert_almost_equal(dC, center.distanceTo(cC) / arcsec, decimal=10)

    # Check that project_rad does the same thing
    uA2, vA2 = center.project_rad(cA.ra.rad, cA.dec.rad, projection='postel')
    np.testing.assert_array_almost_equal(uA, uA2, decimal=10)
    np.testing.assert_array_almost_equal(vA, vA2, decimal=10)

    # Check the deprojection
    cA2 = center.deproject(uA*arcsec, vA*arcsec, projection='postel')
    np.testing.assert_almost_equal(cA.ra.rad, cA2.ra.rad, decimal=12)
    np.testing.assert_almost_equal(cA.dec.rad, cA2.dec.rad, decimal=12)
    cA3 = center.deproject_rad(uA, vA, projection='postel')
    np.testing.assert_array_almost_equal([cA.ra.rad, cA.dec.rad], cA3, decimal=12)

    # The angles are not preserved
    a = sqrt( (uB-uC)**2 + (vB-vC)**2 )
    b = sqrt( (uC-uA)**2 + (vC-vA)**2 )
    c = sqrt( (uA-uB)**2 + (vA-vB)**2 )
    cosA = ((uB-uA)*(uC-uA) + (vB-vA)*(vC-vA)) / (b*c)
    cosB = ((uC-uB)*(uA-uB) + (vC-vB)*(vA-vB)) / (c*a)
    cosC = ((uA-uC)*(uB-uC) + (vA-vC)*(vB-vC)) / (a*b)

    print('postel cosA = ',cosA,cos(A))
    print('postel cosB = ',cosB,cos(B))
    print('postel cosC = ',cosC,cos(C))

    # The area is not preserved
    area = 0.5 * abs( (uB-uA)*(vC-vA) - (uC-uA)*(vB-vA) )
    area *= (arcsec / radians)**2
    print('postel area = ',area,E)

    # The deproject jacobian should tell us how the area changes
    dudx, dudy, dvdx, dvdy = center.jac_deproject(uA*arcsec, vA*arcsec, 'postel').ravel()
    jac_area = abs(dudx*dvdy - dudy*dvdx)
    np.testing.assert_almost_equal(jac_area, E/area, decimal=5)

    dudx, dudy, dvdx, dvdy = center.jac_deproject_arcsec(uA, vA, 'postel').ravel()
    np.testing.assert_almost_equal(jac_area, abs(dudx*dvdy - dudy*dvdx))

    # center projects to 0,0 with unit area
    u, v = center.project(center, 'postel')
    np.testing.assert_almost_equal(u.rad, 0., decimal=12)
    np.testing.assert_almost_equal(v.rad, 0., decimal=12)
    c2 = center.deproject(u,v, 'postel')
    np.testing.assert_almost_equal(c2.ra.rad, center.ra.rad, decimal=12)
    np.testing.assert_almost_equal(c2.dec.rad, center.dec.rad, decimal=12)
    np.testing.assert_almost_equal(np.linalg.det(center.jac_deproject(u, v, 'postel')), 1.)


@timer
def test_precess():
    """Test precession between epochs.
    """
    # I don't have much of a test here.  The formulae are what they are.
    # But it should at least be the case that a precession trip that ends up
    # back at the original epoch should leave the coord unchanged.
    orig = coord.CelestialCoord(0.234 * radians, 0.342 * radians)

    c1 = orig.precess(2000., 1950.)
    c2 = c1.precess(1950., 1900.)
    c3 = c2.precess(1900., 2000.)
    np.testing.assert_almost_equal(c3.ra.rad, orig.ra.rad, decimal=10)
    np.testing.assert_almost_equal(c3.dec.rad, orig.dec.rad, decimal=10)

    # The no op is exact.
    c4 = orig.precess(2000., 2000.)
    assert c4 == orig

    # I found a website that does precession calculations, so check that we are
    # consistent with them.
    # http://www.bbastrodesigns.com/coordErrors.html
    dra_1950 = -(2. + 39.07/60.)/60. * hours / radians
    ddec_1950 = -(16. + 16.3/60.)/60. * degrees / radians
    print('delta from website: ',dra_1950,ddec_1950)
    print('delta from precess: ',(c1.ra-orig.ra),(c1.dec-orig.dec))
    np.testing.assert_almost_equal(dra_1950, c1.ra.rad-orig.ra.rad, decimal=5)
    np.testing.assert_almost_equal(ddec_1950, c1.dec.rad-orig.dec.rad, decimal=5)

    dra_1900 = -(5. + 17.74/60.)/60. * hours / radians
    ddec_1900 = -(32. + 35.4/60.)/60. * degrees / radians
    print('delta from website: ',dra_1900,ddec_1900)
    print('delta from precess: ',(c2.ra-orig.ra),(c2.dec-orig.dec))
    np.testing.assert_almost_equal(dra_1900, c2.ra.rad-orig.ra.rad, decimal=5)
    np.testing.assert_almost_equal(ddec_1900, c2.dec.rad-orig.dec.rad, decimal=5)


@timer
def test_galactic():
    """Test the conversion from equatorial to galactic coordinates and back.
    """
    # According to wikipedia: http://en.wikipedia.org/wiki/Galactic_coordinate_system
    # the galactic center is located at 17h:45.6m, -28.94d
    # But I get more precise values from https://arxiv.org/pdf/1010.3773.pdf
    center = coord.CelestialCoord(
        coord.Angle.from_hms('17:45:37.1991'),
        coord.Angle.from_dms('-28:56:10.2207'))
    print('center.galactic = ',center.galactic())
    el,b = center.galactic()
    np.testing.assert_almost_equal(el.wrap().rad, 0., decimal=8)
    np.testing.assert_almost_equal(b.rad, 0., decimal=8)

    # Go back from galactic coords to CelestialCoord
    center2 = coord.CelestialCoord.from_galactic(el,b)
    np.testing.assert_almost_equal(center2.ra.rad, center.ra.rad, decimal=10)
    np.testing.assert_almost_equal(center2.dec.rad, center.dec.rad, decimal=10)

    # The north pole is at 12h:51.4m, 27.13d again with more precise values from the above paper.
    north = coord.CelestialCoord(
        coord.Angle.from_hms('12:51:26.27549'),
        coord.Angle.from_dms('27:07:41.7043'))
    print('north.galactic = ',north.galactic())
    el,b = north.galactic()
    np.testing.assert_almost_equal(b.rad, pi/2., decimal=8)
    north2 = coord.CelestialCoord.from_galactic(el,b)
    np.testing.assert_almost_equal(north2.ra.rad, north.ra.rad, decimal=10)
    np.testing.assert_almost_equal(north2.dec.rad, north.dec.rad, decimal=10)

    south = coord.CelestialCoord(
        coord.Angle.from_hms('00:51:26.27549'),
        coord.Angle.from_dms('-27:07:41.7043'))
    print('south.galactic = ',south.galactic())
    el,b = south.galactic()
    np.testing.assert_almost_equal(b.rad, -pi/2., decimal=8)
    south2 = coord.CelestialCoord.from_galactic(el,b)
    np.testing.assert_almost_equal(south2.ra.rad, south.ra.rad, decimal=10)
    np.testing.assert_almost_equal(south2.dec.rad, south.dec.rad, decimal=10)

    anticenter = coord.CelestialCoord(
        coord.Angle.from_hms('05:45:37.1991'),
        coord.Angle.from_dms('28:56:10.2207'))
    print('anticenter.galactic = ',anticenter.galactic())
    el,b = anticenter.galactic()
    np.testing.assert_almost_equal(el.rad, pi, decimal=8)
    np.testing.assert_almost_equal(b.rad, 0., decimal=8)
    anticenter2 = coord.CelestialCoord.from_galactic(el,b)
    np.testing.assert_almost_equal(anticenter2.ra.rad, anticenter.ra.rad, decimal=10)
    np.testing.assert_almost_equal(anticenter2.dec.rad, anticenter.dec.rad, decimal=10)


@timer
def test_ecliptic():
    """Test the conversion from equatorial to ecliptic coordinates and back.
    """
    # Use locations of ecliptic poles from http://en.wikipedia.org/wiki/Ecliptic_pole
    north_pole = coord.CelestialCoord(
        coord.Angle.from_hms('18:00:00.00'),
        coord.Angle.from_dms('66:33:38.55'))
    el, b = north_pole.ecliptic()
    # North pole should have b=90 degrees, with el being completely arbitrary.
    np.testing.assert_almost_equal(b.rad, pi/2, decimal=6)

    # Check round-trip
    north_pole2 = coord.CelestialCoord.from_ecliptic(el, b)
    np.testing.assert_almost_equal(north_pole2.ra.rad, north_pole.ra.rad, decimal=8)
    np.testing.assert_almost_equal(north_pole2.dec.rad, north_pole.dec.rad, decimal=8)

    south_pole = coord.CelestialCoord(
        coord.Angle.from_hms('06:00:00.00'),
        coord.Angle.from_dms('-66:33:38.55'))
    el, b = south_pole.ecliptic()
    np.testing.assert_almost_equal(b.rad, -pi/2, decimal=6)
    south_pole2 = coord.CelestialCoord.from_ecliptic(el, b)
    np.testing.assert_almost_equal(south_pole2.ra.rad, south_pole.ra.rad, decimal=8)
    np.testing.assert_almost_equal(south_pole2.dec.rad, south_pole.dec.rad, decimal=8)

    # Also confirm that positions that should be the same in equatorial and ecliptic coordinates are
    # actually the same:
    vernal = coord.CelestialCoord(0.*radians, 0.*radians)
    el, b = vernal.ecliptic()
    np.testing.assert_almost_equal(b.rad, 0., decimal=6)
    np.testing.assert_almost_equal(el.rad, 0., decimal=6)
    vernal2 = coord.CelestialCoord.from_ecliptic(el,b)
    np.testing.assert_almost_equal(vernal2.ra.rad, vernal.ra.rad, decimal=8)
    np.testing.assert_almost_equal(vernal2.dec.rad, vernal.dec.rad, decimal=8)

    autumnal = coord.CelestialCoord(pi*radians, 0.*radians)
    el, b = autumnal.ecliptic()
    np.testing.assert_almost_equal(el.wrap(pi*radians).rad, pi, decimal=6)
    np.testing.assert_almost_equal(b.rad, 0., decimal=6)
    autumnal2 = coord.CelestialCoord.from_ecliptic(el,b)
    np.testing.assert_almost_equal(autumnal2.ra.rad, autumnal.ra.rad, decimal=8)
    np.testing.assert_almost_equal(autumnal2.dec.rad, autumnal.dec.rad, decimal=8)


@timer
def test_ecliptic_date():
    """Test the date option of the ecliptic and from_ecliptic functions.
    """
    # Test the results of using a date to get ecliptic coordinates with respect to the sun,
    # instead of absolute ones. For this, use dates and times of vernal and autumnal equinox
    # in 2014 from
    #   http://wwp.greenwichmeantime.com/longest-day/
    # and the conversion to Julian dates from
    #   http://www.aavso.org/jd-calculator

    import datetime
    vernal_eq_date = datetime.datetime(2014,3,20,16,57,0)
    vernal = coord.CelestialCoord(0.*radians, 0.*radians)
    el, b = vernal.ecliptic(epoch=2014)
    el_rel, b_rel = vernal.ecliptic(epoch=2014, date=vernal_eq_date)
    # Vernal equinox: should have (el, b) = (el_rel, b_rel) = 0.0
    np.testing.assert_almost_equal(el.rad, el_rel.rad, decimal=4)
    np.testing.assert_almost_equal(b.rad, b_rel.rad, decimal=6)
    vernal2 = coord.CelestialCoord.from_ecliptic(el_rel, b_rel, date=vernal_eq_date)
    np.testing.assert_almost_equal(vernal2.ra.wrap().rad, vernal.ra.rad, decimal=8)
    np.testing.assert_almost_equal(vernal2.dec.rad, vernal.dec.rad, decimal=8)

    # Now do the autumnal equinox: should have (el, b) = (pi, 0) = (el_rel, b_rel) when we look at
    # the time of the vernal equinox.
    autumnal = coord.CelestialCoord(pi*radians, 0.*radians)
    el, b = autumnal.ecliptic(epoch=2014)
    el_rel, b_rel = autumnal.ecliptic(epoch=2014, date=vernal_eq_date)
    np.testing.assert_almost_equal(el_rel.wrap(pi*radians).rad, el.wrap(pi*radians).rad,
                                   decimal=4)
    np.testing.assert_almost_equal(b_rel.rad, b.rad, decimal=6)
    autumnal2 = coord.CelestialCoord.from_ecliptic(el_rel, b_rel, date=vernal_eq_date)
    np.testing.assert_almost_equal(autumnal2.ra.wrap(pi*radians).rad,
                                   autumnal.ra.wrap(pi*radians).rad, decimal=8)
    np.testing.assert_almost_equal(autumnal2.dec.rad, autumnal.dec.rad, decimal=8)

    # And check that if it's the date of the autumnal equinox (sun at (180, 0)) but we're looking at
    # the position of the vernal equinox (0, 0), then (el_rel, b_rel) = (-180, 0)
    autumnal_eq_date = datetime.datetime(2014,9,23,2,29,0)
    el_rel, b_rel = vernal.ecliptic(epoch=2014, date=autumnal_eq_date)
    np.testing.assert_almost_equal(el_rel.wrap(-pi*radians).rad, -pi, decimal=4)
    np.testing.assert_almost_equal(b_rel.rad, 0., decimal=6)

    # And check that if it's the date of the vernal equinox (sun at (0, 0)) but we're looking at
    # the position of the autumnal equinox (180, 0), then (el_rel, b_rel) = (180, 0)
    el_rel, b_rel = autumnal.ecliptic(epoch=2014, date=vernal_eq_date)
    np.testing.assert_almost_equal(el_rel.wrap(pi*radians).rad, pi, decimal=4)
    np.testing.assert_almost_equal(b_rel.rad, 0., decimal=6)

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
