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

    # Access angle in other units
    np.testing.assert_almost_equal(theta1 / gradians, 50., 12, 'pi/4 rad / gradians != 50')
    np.testing.assert_almost_equal(theta3 / gradians, 50., 12, '3 hours / gradians != 50')
    np.testing.assert_almost_equal(theta1 / coord.hours, 3., 12, 'pi/4 rad / hours != 3')
    np.testing.assert_almost_equal(theta5 / coord.hours, 3., 12, '45*3600 arcsec / hours != 3')

@timer
def test_invalid():
    """Test invalid ways to try to make an Angle
    """
    pass


@timer
def test_arith():
    """Test various arithmetic operations on Angles
    """
    pass


@timer
def test_wrap():
    """Test the Angle.wrap function
    """
    pass


@timer
def test_trig():
    """Test the various trig functions with Angles
    """
    pass


@timer
def test_pickle():
    """Test that Angle instances pickle correctly
    """
    pass


@timer
def test_eq():
    """Check that equal angles are equal, but unequal ones are not.
    """
    theta1 = pi/4. * coord.radians
    theta2 = 45 * coord.degrees
    assert theta1 == theta2, "pi/4 rad != 45 deg"

    theta3 = theta1.wrap()
    assert theta1 == theta3, "(pi/4 rad).wrap() != pi/4 rad"

    theta4 = coord.Angle(theta1)  # Copy constructor
    assert theta4 == theta1, "copy not == to original"

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
    pass


@timer
def test_dms():
    """Test Angle.dms and from_dms functions
    """
    pass


if __name__ == '__main__':
    test_init()
    test_invalid()
    test_arith()
    test_wrap()
    test_trig()
    test_pickle()
    test_eq()
    test_hms()
    test_dms()
