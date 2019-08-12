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

# We'll use this a lot, so just import it.
from math import pi

@timer
def test_init():
    """Test basic construction and use of Angle classes
    """
    # Various ways to make 45 degrees
    theta1 = pi/4. * coord.radians
    theta2 = 45 * coord.degrees
    theta3 = coord.Angle(3, coord.hours)
    theta4 = coord.Angle(theta=45 * 60, unit=coord.arcmin)
    theta5 = coord.Angle(45. * 3600, unit=coord.arcsec)
    gradians = coord.AngleUnit(2. * pi / 400.)
    theta6 = 50 * gradians
    theta7 = coord.Angle(theta1)  # Copy constructor

    # Equivalent ways to access theta in radians
    assert theta1.rad == pi/4., "(pi/4 radians).rad != pi/4"
    assert theta1 / coord.radians == pi/4., "(pi/4 radians) / radians != pi/4"

    # Others involved a calculation so aren't required to be precisely equal.
    np.testing.assert_almost_equal(theta2.rad, pi/4., decimal=12, err_msg='45 degrees != pi/4 rad')
    np.testing.assert_almost_equal(theta3.rad, pi/4., decimal=12, err_msg='3 hours != pi/4 rad')
    np.testing.assert_almost_equal(theta4.rad, pi/4., decimal=12, err_msg='45*60 arcmin != pi/4 rad')
    np.testing.assert_almost_equal(theta5.rad, pi/4., decimal=12, err_msg='45*3600 arcsec != pi/4 rad')
    np.testing.assert_almost_equal(theta6.rad, pi/4., decimal=12, err_msg='50 grad != pi/4 rad')
    np.testing.assert_almost_equal(theta7.rad, pi/4., decimal=12, err_msg='copy != pi/4 rad')

    for theta in [theta1, theta2, theta3, theta4, theta5, theta6, theta7]:
        # Access angle in degrees
        np.testing.assert_almost_equal(theta.deg, 45., 12)
        np.testing.assert_almost_equal(theta / coord.degrees, 45., 12)

        # Access angle in other units
        np.testing.assert_almost_equal(theta / gradians, 50., 12)
        np.testing.assert_almost_equal(theta / coord.hours, 3., 12)
        np.testing.assert_almost_equal(theta / coord.arcmin, 45.*60, 11)
        np.testing.assert_almost_equal(theta / coord.arcsec, 45.*60*60, 10)


@timer
def test_invalid():
    """Test invalid ways to try to make an Angle
    """
    theta = 1.3 * coord.degrees

    # No unit
    np.testing.assert_raises(TypeError, coord.Angle, 3.4)

    # Invalid value type
    np.testing.assert_raises(TypeError, coord.Angle, theta, coord.degrees)
    np.testing.assert_raises(ValueError, coord.Angle, 'spam', coord.degrees)

    # Wrong order
    np.testing.assert_raises(TypeError, coord.Angle, coord.degrees, 1.3)

    # Invalid unit type
    np.testing.assert_raises(TypeError, coord.Angle, 1.3, 3)
    np.testing.assert_raises(TypeError, coord.Angle, 1.3, theta)

    # Invalid value type for copy constructor
    np.testing.assert_raises(TypeError, coord.Angle, 1.3)


@timer
def test_arith():
    """Test various arithmetic operations on Angles
    """
    theta1 = 45. * coord.degrees
    theta2 = 135. * coord.degrees
    np.testing.assert_almost_equal((theta1 + theta2).rad, pi, decimal=12)
    np.testing.assert_almost_equal((theta2 - theta1).rad, pi/2., decimal=12)
    np.testing.assert_almost_equal((4*theta1).rad, pi, decimal=12)
    np.testing.assert_almost_equal((theta2*2).rad, 1.5*pi, decimal=12)
    np.testing.assert_almost_equal((4*theta1 - theta2).rad, pi/4., decimal=12)
    np.testing.assert_almost_equal((theta1/2.).rad, pi/8., decimal=12)
    np.testing.assert_almost_equal((-theta1).rad, -pi/4., decimal=12)

    theta2 += theta1
    np.testing.assert_almost_equal(theta2.rad, pi, decimal=12)
    theta2 -= 2*theta1
    np.testing.assert_almost_equal(theta2.rad, pi/2., decimal=12)
    theta2 *= 4
    np.testing.assert_almost_equal(theta2.rad, 2.*pi, decimal=12)
    theta2 /= -2
    np.testing.assert_almost_equal(theta2.rad, -pi, decimal=12)

    # Check invalid arithmetic
    np.testing.assert_raises(TypeError, coord.Angle.__add__, theta1, 23)
    np.testing.assert_raises(TypeError, coord.Angle.__add__, theta1, '23')
    np.testing.assert_raises(TypeError, coord.Angle.__add__, theta1, coord.degrees)
    np.testing.assert_raises(TypeError, coord.Angle.__sub__, theta1, 23)
    np.testing.assert_raises(TypeError, coord.Angle.__sub__, theta1, '23')
    np.testing.assert_raises(TypeError, coord.Angle.__sub__, theta1, coord.degrees)
    np.testing.assert_raises(TypeError, coord.Angle.__mul__, theta1, '23')
    np.testing.assert_raises(TypeError, coord.Angle.__mul__, theta1, 23*coord.degrees)
    np.testing.assert_raises(TypeError, coord.Angle.__mul__, theta1, coord.degrees)
    np.testing.assert_raises(TypeError, coord.Angle.__div__, theta1, 23*coord.degrees)
    np.testing.assert_raises(TypeError, coord.Angle.__div__, theta1, '23')

@timer
def test_wrap():
    """Test the Angle.wrap function
    """
    # Two values that differ by a full cycle
    theta1 = 45 * coord.degrees
    theta2 = (45 + 360) * coord.degrees

    # Not wrapped intrinsically.  They differ by 2pi
    assert 6.28 < theta2.rad - theta1.rad < 6.29

    # Wrapping theta2 makes them equal
    print('theta1 = ',theta1,' => ',theta1.wrap())
    print('theta2 = ',theta2,' => ',theta2.wrap())
    np.testing.assert_almost_equal(theta2.wrap().rad, theta1.rad, decimal=12)

    # Any multiple of 360 degrees is equivalent
    theta3 = (45 - 360) * coord.degrees
    assert -6.29 < theta3.rad - theta1.rad < -6.28
    np.testing.assert_almost_equal(theta3.wrap().rad, theta1.rad, decimal=12)

    theta4 = (45 + 10 * 360) * coord.degrees
    assert 62.8 < theta4.rad - theta1.rad < 62.9
    np.testing.assert_almost_equal(theta4.wrap().rad, theta1.rad, decimal=12)

    # Can wrap around other angles besides 0.
    pir = pi * coord.radians
    np.testing.assert_almost_equal(theta2.wrap(pir).rad, theta1.rad, decimal=12)
    np.testing.assert_almost_equal(theta2.rad, theta1.wrap(2*pir).rad, decimal=12)
    np.testing.assert_almost_equal(theta2.rad, theta1.wrap(3*pir).rad, decimal=12)
    np.testing.assert_almost_equal(theta3.rad, theta1.wrap(-pir).rad, decimal=12)
    np.testing.assert_almost_equal(theta3.rad, theta1.wrap(-2*pir).rad, decimal=12)
    np.testing.assert_almost_equal(theta2.wrap(27 * coord.radians).rad,
                                   theta1.wrap(27 * coord.radians).rad, decimal=12)
    np.testing.assert_almost_equal(theta3.wrap(-127 * coord.radians).rad,
                                   theta1.wrap(-127 * coord.radians).rad, decimal=12)
    np.testing.assert_almost_equal(theta4.wrap(-927 * coord.radians).rad,
                                   theta1.wrap(-927 * coord.radians).rad, decimal=12)

@timer
def test_compare():
    """Test comparisons of Angles.
    """
    # Values in order
    theta1 = -45 * coord.degrees
    theta2 = 45 * coord.degrees
    theta2b = 45 * coord.degrees
    theta3 = 47 * coord.degrees
    theta4 = (45 + 360) * coord.degrees
    theta4b= (45 + 360) * coord.degrees
    theta5 = (47 + 360) * coord.degrees

    assert theta1 < theta2
    assert theta1 < theta3
    assert theta1 < theta4
    assert theta1 < theta5

    assert theta1 <= theta2
    assert theta1 <= theta3
    assert theta1 <= theta4
    assert theta1 <= theta5

    assert theta5 > theta1
    assert theta5 > theta2
    assert theta5 > theta3
    assert theta5 > theta4

    assert theta5 >= theta1
    assert theta5 >= theta2
    assert theta5 >= theta3
    assert theta5 >= theta4

    assert theta2 == theta2b
    assert theta2 <= theta2b
    assert theta2 >= theta2b

    assert theta2 == theta2
    assert theta2 <= theta2
    assert theta2 >= theta2

    assert theta4 == theta4b
    assert theta4 <= theta4b
    assert theta4 >= theta4b

    # Have to explicitly wrap if you want it.
    assert theta2.wrap() < theta5.wrap()
    assert theta5.wrap() > theta2.wrap()
    assert theta3.wrap() > theta4.wrap()
    assert theta4.wrap() < theta3.wrap()


@timer
def test_trig():
    """Test the various trig functions with Angles
    """
    theta1 = 45 * coord.degrees
    np.testing.assert_almost_equal(theta1.sin(), math.sqrt(0.5), decimal=12)
    np.testing.assert_almost_equal(theta1.cos(), math.sqrt(0.5), decimal=12)
    np.testing.assert_almost_equal(theta1.tan(), 1., decimal=12)
    np.testing.assert_almost_equal(theta1.sincos(), math.sqrt(0.5), decimal=12)

    theta2 = 30 * coord.degrees
    np.testing.assert_almost_equal(theta2.sin(), 0.5, decimal=12)
    np.testing.assert_almost_equal(theta2.cos(), math.sqrt(0.75), decimal=12)
    np.testing.assert_almost_equal(theta2.tan(), math.sqrt(1./3.), decimal=12)
    np.testing.assert_almost_equal(theta2.sincos(), (0.5, math.sqrt(0.75)), decimal=12)

    theta3 = 59 * coord.degrees
    theta4 = 31 * coord.degrees
    np.testing.assert_almost_equal(theta3.sin(), theta4.cos(), decimal=12)
    np.testing.assert_almost_equal(theta3.cos(), theta4.sin(), decimal=12)
    np.testing.assert_almost_equal(theta3.tan(), 1./theta4.tan(), decimal=12)
    np.testing.assert_almost_equal(theta3.sincos(), theta4.sincos()[::-1], decimal=12)

@timer
def test_pickle():
    """Test that Angle instances pickle correctly
    """
    theta1 = 45 * coord.degrees
    theta2 = 29.9 * coord.radians
    theta3 = 1.e-4 * coord.arcsec
    do_pickle(theta1)
    do_pickle(theta2)
    do_pickle(theta3)

@timer
def test_eq():
    """Check that equal angles are equal, but unequal ones are not.
    """
    theta1 = pi/4. * coord.radians
    theta2 = 45 * coord.degrees
    assert theta1 == theta2, "pi/4 rad != 45 deg"

    theta3 = coord.Angle(theta1)  # Copy constructor
    assert theta3 == theta1, "copy not == to original"

    theta4 = theta1.wrap()
    assert theta1 == theta4, "(pi/4 rad).wrap() != pi/4 rad"

    # These should all test as unequal.  Note some non-Angles in the list.
    diff_list = [ theta1,
                  14 * coord.degrees,
                  -theta1,
                  theta1 * 2.,
                  theta1 + 360. * coord.degrees,
                  theta1 - 360. * coord.degrees,
                  pi/4.,
                  coord.Angle,
                  None ]
    all_obj_diff(diff_list)

@timer
def test_hms():
    """Test Angle.hms and from_hms functions
    """
    theta1 = 45 * coord.degrees
    hms1 = theta1.hms()
    print('hms1 = ',hms1)
    assert hms1 == '03:00:00'

    theta2 = 0 * coord.degrees
    hms2 = theta2.hms(plus_sign=True, prec=8, sep='hm')
    print('hms2 = ',hms2)
    assert hms2 == '+00h00m00.00000000'

    # Check that 360 degres is 24 hours
    theta3 = 360 * coord.degrees
    hms3 = theta3.hms(sep='hms', prec=0)
    print('hms3 = ',hms3)
    assert hms3 == '24h00m00s'

    # Check a negative small value
    theta4 = -0.2 * coord.degrees
    hms4 = theta4.hms(pad=False, prec=2, sep=(' hours ',' minutes ',' seconds'))
    print('hms4 = ',hms4)
    assert hms4 == '-0 hours 0 minutes 48.00 seconds'

    # Some random value
    theta5 = 1.3 * coord.radians
    hms5 = theta5.hms(prec=12)
    print('hms5 = ',hms5)
    h,m,s = hms5.split(':')
    theta5b = float(h) + float(m)/60. + float(s)/3600.
    np.testing.assert_almost_equal(theta5 / coord.hours, theta5b, decimal=12)

    # Negative angles < 1 degree are tricky
    theta6 = -6.112479 * coord.arcmin
    hms6 = theta6.hms(prec=2)
    print('hms6 = ',hms6)
    assert hms6 == '-00:00:24.45'
    alt_theta6 = -6.1125 * coord.arcmin

    # Check round trip
    np.testing.assert_almost_equal(theta1.rad, coord.Angle.from_hms(hms1).rad, decimal=12)
    np.testing.assert_almost_equal(theta2.rad, coord.Angle.from_hms(hms2).rad, decimal=12)
    np.testing.assert_almost_equal(theta3.rad, coord.Angle.from_hms(hms3).rad, decimal=12)
    np.testing.assert_almost_equal(theta4.rad, coord.Angle.from_hms(hms4).rad, decimal=12)
    np.testing.assert_almost_equal(theta5.rad, coord.Angle.from_hms(hms5).rad, decimal=12)
    np.testing.assert_almost_equal(alt_theta6.rad, coord.Angle.from_hms(hms6).rad, decimal=12)

    # This one doesn't round trip, but it is valid -> hms
    theta7 = 1.2345 * coord.radians
    hms7 = theta7.hms(sep='', prec=0)
    print('hms7 = ',hms7)
    assert hms7 == '044256'

    # These aren't constructible from hms() but they parse correctly
    theta8 = 17.3 * coord.hours
    hms8a = '17h18'
    hms8b = '+17h18m'
    hms8c = '17.3h'
    hms8d = '+17.1h12m'  # Weird, but allowed.  Maybe it shouldn't be.
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_hms(hms8a).rad, decimal=12)
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_hms(hms8b).rad, decimal=12)
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_hms(hms8c).rad, decimal=12)
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_hms(hms8d).rad, decimal=12)

    # Check invalid calls
    np.testing.assert_raises(ValueError, theta1.hms, sep='hmsss')
    np.testing.assert_raises(TypeError, theta1.hms, sep=coord.hours)
    np.testing.assert_raises(ValueError, theta1.hms, prec=-1)
    np.testing.assert_raises(TypeError, theta1.hms, prec='lots')
    np.testing.assert_raises(ValueError, coord.Angle.from_hms, 'a')
    np.testing.assert_raises(ValueError, coord.Angle.from_hms, '-')
    np.testing.assert_raises(ValueError, coord.Angle.from_hms, '+')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, '-17')
    np.testing.assert_raises(ValueError, coord.Angle.from_hms, '')
    np.testing.assert_raises(ValueError, coord.Angle.from_hms, '01:21:31:15')



@timer
def test_dms():
    """Test Angle.dms and from_dms functions
    """
    theta1 = 45 * coord.degrees
    dms1 = theta1.dms()
    print('dms1 = ',dms1)
    assert dms1 == '45:00:00'

    theta2 = 0 * coord.degrees
    dms2 = theta2.dms(plus_sign=True, prec=8, sep='dm')
    print('dms2 = ',dms2)
    assert dms2 == '+00d00m00.00000000'

    theta3 = 360 * coord.degrees
    dms3 = theta3.dms(sep='dms', prec=0)
    print('dms3 = ',dms3)
    assert dms3 == '360d00m00s'

    theta4 = -1 * coord.degrees
    dms4 = theta4.dms()
    dms4 = theta4.dms(pad=False, prec=2, sep=(' degrees ',' minutes ',' seconds'))
    print('dms4 = ',dms4)
    assert dms4 == '-1 degrees 0 minutes 0.00 seconds'

    theta5 = 1.3 * coord.radians
    dms5 = theta5.dms(prec=12)
    print('dms5 = ',dms5)
    d,m,s = dms5.split(':')
    theta5b = float(d) + float(m)/60. + float(s)/3600.
    np.testing.assert_almost_equal(theta5 / coord.degrees, theta5b, decimal=12)

    # Negative angles < 1 degree are tricky
    theta6 = -6.112479 * coord.arcmin
    dms6 = theta6.dms(prec=2)
    print('dms6 = ',dms6)
    assert dms6 == '-00:06:06.75'
    alt_theta6 = -6.1125 * coord.arcmin

    # Check round trip
    np.testing.assert_almost_equal(theta1.rad, coord.Angle.from_dms(dms1).rad, decimal=12)
    np.testing.assert_almost_equal(theta2.rad, coord.Angle.from_dms(dms2).rad, decimal=12)
    np.testing.assert_almost_equal(theta3.rad, coord.Angle.from_dms(dms3).rad, decimal=12)
    np.testing.assert_almost_equal(theta4.rad, coord.Angle.from_dms(dms4).rad, decimal=12)
    np.testing.assert_almost_equal(theta5.rad, coord.Angle.from_dms(dms5).rad, decimal=12)
    np.testing.assert_almost_equal(alt_theta6.rad, coord.Angle.from_dms(dms6).rad, decimal=12)

    # This one doesn't round trip, but it is valid -> hms
    theta7 = 1.2345 * coord.radians
    dms7 = theta7.dms(sep='', prec=0)
    print('dms7 = ',dms7)
    assert dms7 == '704354'

    # These aren't constructible from dms() but they parse correctly
    theta8 = 17.3 * coord.degrees
    dms8a = '+17d18'
    dms8b = '17d18m'
    dms8c = '+17.3d'
    dms8d = '17.1d12m'  # Weird, but allowed.  Maybe it shouldn't be.
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_dms(dms8a).rad, decimal=12)
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_dms(dms8b).rad, decimal=12)
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_dms(dms8c).rad, decimal=12)
    np.testing.assert_almost_equal(theta8.rad, coord.Angle.from_dms(dms8d).rad, decimal=12)

    # Check invalid calls
    np.testing.assert_raises(ValueError, theta1.dms, sep='dmsss')
    np.testing.assert_raises(TypeError, theta1.dms, sep=coord.degrees)
    np.testing.assert_raises(ValueError, theta1.dms, prec=-1)
    np.testing.assert_raises(TypeError, theta1.dms, prec='lots')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, 'a')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, '-')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, '+')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, '+17')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, '')
    np.testing.assert_raises(ValueError, coord.Angle.from_dms, '01:21:31:15')


if __name__ == '__main__':
    test_init()
    test_invalid()
    test_arith()
    test_wrap()
    test_compare()
    test_trig()
    test_pickle()
    test_eq()
    test_hms()
    test_dms()
