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
    pass


@timer
def test_builtin_units():
    """Check the built-in AngleUnits
    """
    pass


@timer
def test_invalid():
    """Check that invalid constructors raise the appropriate exceptions.
    """
    pass


@timer
def test_div():
    """Test AngleUnit / AngleUnit
    """
    pass


@timer
def test_pickle():
    """Check that everything pickles correctly.

    Not a big deal here, but if we start including any C structs, they require special handling
    to make picklable, so useful to have this check in place.

    More usefully, the do_pickle function also checks that eval(repr(obj)) returns something
    equivalent to the original object along with a few other sanity checks such as obj != None.
    So it's a good idea to just run this on every kind of object we have.
    """
    pass


@timer
def test_from_name():
    """Test the AngleUnit from_name static method

    Note in particular that shorter or longer versions are allowed.  As are caps or lowercase.
    """
    pass


if __name__ == '__main__':
    test_init()
    test_builtin_units()
    test_invalid()
    test_div()
    test_pickle()
    test_from_name()
