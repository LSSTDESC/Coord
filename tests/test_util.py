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
import datetime

from helper_util import *

# Note: sun_position_ecliptic and ecliptic_obliquity in util.py are implicitly tested
# by the test_ecliptic_date tests.

@timer
def test_julian():
    """Test date_to_julian_day function
    """
    # Use Julian Day Number from Date Calculator here:
    #    http://keisan.casio.com/exec/system/1227779487
    # The Julian Day for noon on 1/1/2000 is also mentioned here:
    #    https://en.wikipedia.org/wiki/Julian_day
    dates = [ datetime.date(1582, 10, 14),
              datetime.date(2000, 1, 1),
              datetime.datetime(2000, 1, 1, 12, 0, 0),
              datetime.datetime(2000, 1, 1, 0, 0, 0),
              datetime.datetime(2000, 1, 1, 18, 0, 0)
            ]
    julian_days = [ 2299160, 2451545, 2451545, 2451544.5, 2451545.25 ]

    for date, jd in zip(dates, julian_days):
        jd1 = coord.util.date_to_julian_day(date)
        np.testing.assert_almost_equal(jd1, jd, 6, 'Julian day not as expected for %s'%date)

    np.testing.assert_raises(ValueError, coord.util.date_to_julian_day, 2000)
    np.testing.assert_raises(ValueError, coord.util.date_to_julian_day, '1875 AD')



if __name__ == '__main__':
    test_julian()
