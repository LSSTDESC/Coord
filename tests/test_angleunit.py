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
    """Test the AngleUnit initializer
    """
    # This is from the AngleUnit doc string.  Make sure it works!
    gradians = coord.AngleUnit(2.*pi / 400.)
    assert gradians.value == 2.*pi / 400., 'gradians value not 2pi/400'

    # Can also use named keyword argument
    np.testing.assert_almost_equal(coord.AngleUnit(value=17).value, 17, decimal=12)

    # Other types are ok as the value argument so long as they are convertible to float.
    np.testing.assert_almost_equal(coord.AngleUnit(np.float64(0.17)).value, 0.17, decimal=12,
                                   err_msg='using np.float64 for value failed')
    np.testing.assert_almost_equal(coord.AngleUnit(np.float32(0.23)).value, 0.23, decimal=8,
                                   err_msg='using np.float32 for value failed')
    np.testing.assert_almost_equal(coord.AngleUnit('1.7').value, 1.7, decimal=12,
                                   err_msg='using string "1.7" for value failed')


@timer
def test_builtin_units():
    """Check the built-in AngleUnits
    """
    np.testing.assert_almost_equal(coord.radians.value, 1., decimal=12)
    np.testing.assert_almost_equal(coord.degrees.value, pi / 180., decimal=12)
    np.testing.assert_almost_equal(coord.hours.value, pi / 12., decimal=12)
    np.testing.assert_almost_equal(coord.arcmin.value, pi / 180. / 60., decimal=12)
    np.testing.assert_almost_equal(coord.arcsec.value, pi / 180. / 3600., decimal=12)


@timer
def test_invalid():
    """Check that invalid constructors raise the appropriate exceptions.
    """
    # Wrong type of value argument
    np.testing.assert_raises(TypeError, coord.AngleUnit, coord.degrees)
    # Also wrong type, but strings give a ValueError
    np.testing.assert_raises(ValueError, coord.AngleUnit, 'spam')
    # Wrong number of arguments
    np.testing.assert_raises(TypeError, coord.AngleUnit, 1, 3)
    np.testing.assert_raises(TypeError, coord.AngleUnit)
    # Wrong keyword argument
    np.testing.assert_raises(TypeError, coord.AngleUnit, the_value=0.2)


@timer
def test_div():
    """Test AngleUnit / AngleUnit
    """
    np.testing.assert_almost_equal(coord.degrees / coord.arcmin, 60., decimal=12)
    np.testing.assert_almost_equal(coord.degrees / coord.arcsec, 3600., decimal=12)
    np.testing.assert_almost_equal(coord.hours / coord.degrees, 15., decimal=12)
    np.testing.assert_almost_equal(coord.hours / coord.hours, 1., decimal=12)

    # Dividing by radians is just the value attribute.
    np.testing.assert_almost_equal(coord.degrees / coord.radians, coord.degrees.value, decimal=12)
    np.testing.assert_almost_equal(coord.arcsec / coord.radians, coord.arcsec.value, decimal=12)

    # We claim this is equivalent to (1 * AngleUnit) / AngleUnit.  Check.
    np.testing.assert_almost_equal((1. * coord.degrees) / coord.arcmin, 60., decimal=12)
    np.testing.assert_almost_equal((1. * coord.degrees) / coord.arcsec, 3600., decimal=12)

    # In general, it is associative.
    np.testing.assert_almost_equal((23.17 * coord.hours) / coord.radians,
                                   23.17 * (coord.hours / coord.radians), decimal=12)

    # Invalid to divide by something other than an AngleUnit
    np.testing.assert_raises(TypeError, coord.degrees.__div__, 12)
    np.testing.assert_raises(TypeError, coord.degrees.__div__, 12 * coord.arcsec)
    np.testing.assert_raises(TypeError, coord.degrees.__div__, 'arcmin')


@timer
def test_pickle():
    """Check that AngleUnits pickle correctly.

    Not a big deal here, but if we start including any C structs, they require special handling
    to make picklable, so useful to have this check in place.

    More usefully, the do_pickle function also checks that eval(repr(obj)) returns something
    equivalent to the original object along with a few other sanity checks such as obj != None.
    So it's a good idea to just run this on every kind of object we have.
    """
    do_pickle(coord.radians)
    do_pickle(coord.degrees)
    do_pickle(coord.hours)
    do_pickle(coord.arcmin)
    do_pickle(coord.arcsec)
    gradians = coord.AngleUnit(2.*pi / 400.)
    do_pickle(gradians)


@timer
def test_eq():
    """Check that equal units are equal, but unequal ones are not.
    """
    # The alt versions should be equal to the built-in versions.
    alt_radians = coord.AngleUnit(1.)
    alt_degrees = coord.AngleUnit(pi/180.)
    assert coord.radians == alt_radians
    assert coord.degrees == alt_degrees

    # These should all test as unequal.  Note some non-AngleUnits in the list.
    diff_list = [ coord.radians, alt_degrees, coord.hours, coord.arcmin, coord.arcsec,
                  1.0, coord.AngleUnit, None ]
    all_obj_diff(diff_list)


@timer
def test_from_name():
    """Test the AngleUnit from_name static method

    Note in particular that shorter or longer versions are allowed.  As are caps or lowercase.
    """
    assert coord.AngleUnit.from_name('radians') == coord.radians
    assert coord.AngleUnit.from_name('Radian') == coord.radians
    assert coord.AngleUnit.from_name('rad') == coord.radians
    assert coord.AngleUnit.from_name('RADians FTW') == coord.radians  # unlikely usage, but valid.
    assert coord.AngleUnit.from_name('degrees') == coord.degrees
    assert coord.AngleUnit.from_name('degree') == coord.degrees
    assert coord.AngleUnit.from_name('degr') == coord.degrees
    assert coord.AngleUnit.from_name('DEG') == coord.degrees
    assert coord.AngleUnit.from_name('hours') == coord.hours
    assert coord.AngleUnit.from_name('hour') == coord.hours
    assert coord.AngleUnit.from_name('hrs') == coord.hours
    assert coord.AngleUnit.from_name('hr') == coord.hours
    assert coord.AngleUnit.from_name('arcmin') == coord.arcmin
    assert coord.AngleUnit.from_name('arcminutes') == coord.arcmin
    assert coord.AngleUnit.from_name('arcsec') == coord.arcsec
    assert coord.AngleUnit.from_name('arcseconds') == coord.arcsec

    np.testing.assert_raises(ValueError, coord.AngleUnit.from_name, 'gradians')
    np.testing.assert_raises(ValueError, coord.AngleUnit.from_name, 'spam')

    # Make sure all the angles in AngleUnit.valid_names work
    for name in coord.AngleUnit.valid_names:
        assert isinstance(coord.AngleUnit.from_name(name), coord.AngleUnit)

if __name__ == '__main__':
    test_init()
    test_builtin_units()
    test_invalid()
    test_div()
    test_pickle()
    test_eq()
    test_from_name()
