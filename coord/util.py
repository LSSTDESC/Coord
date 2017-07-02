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

"""
Some utility functions used by the coord module.
"""

import numpy as np
import math
import datetime

from .angleunit import degrees
from .angle import Angle

def sun_position_ecliptic(date):
    """Helper routine to calculate the position of the sun in ecliptic coordinates given a
    python datetime object.

    It is most precise for dates between 1950-2050, and is based on

        http://en.wikipedia.org/wiki/Position_of_the_Sun#Ecliptic_coordinates

    :param date:    The date as either a datetime.datetime instance or a datetime.date instance.

    :returns the angular position of the sun along the ecliptic.
    """
    # We start by getting the number of days since Greenwich noon on 1 January 2000 (J2000).
    jd = date_to_julian_day(date)
    n = jd - 2451545.0
    L = 280.46*degrees + (0.9856474*degrees) * n
    g = 357.528*degrees + (0.9856003*degrees) * n
    lam = L + (1.915*degrees)*g.sin() + (0.020*degrees)*(2*g).sin()
    return lam

def date_to_julian_day(date):
    """Helper routine to return the Julian day for a given date.

    If `date` is a datetime.datetime instance, then it uses the full time info.
    If `date` is a datetime.date, then it does the calculation for noon of that day.

    :returns: the (possibly fractional) Julian day for the given date.
    """
    # From http://code-highlights.blogspot.com/2013/01/julian-date-in-python.html
    if not (isinstance(date, datetime.date) or isinstance(date, datetime.datetime)):
        raise ValueError("Date must be a python datetime object!")
    a = (14. - date.month)//12
    y = date.year + 4800 - a
    m = date.month + 12*a - 3
    retval = date.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045
    if isinstance(date, datetime.datetime):
        dayfrac = (date.hour + date.minute/60. + date.second/3600.)/24
        # The default is the value at noon, so we want to add 0 if dayfrac = 0.5
        dayfrac -= 0.5
        retval += dayfrac
    return retval

def ecliptic_obliquity(epoch):
    """Helper routine to return the obliquity of the ecliptic for a given date.

    :param epoch:   The epoch at which to calculate the obliquity.

    :returns the obliquity as an Angle instance
    """
    # We need to figure out the time in Julian centuries from J2000 for this epoch.
    t = (epoch - 2000.)/100.
    # Then we use the last (most recent) formula listed under
    # http://en.wikipedia.org/wiki/Ecliptic#Obliquity_of_the_ecliptic, from
    # JPL's 2010 calculations.
    ep = Angle.from_dms('23:26:21.406')
    ep -= Angle.from_dms('00:00:46.836769')*t
    ep -= Angle.from_dms('00:00:0.0001831')*(t**2)
    ep += Angle.from_dms('00:00:0.0020034')*(t**3)
    # There are even higher order terms, but they are probably not important for any reasonable
    # calculation someone would do with this package.
    return ep
