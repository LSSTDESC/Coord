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
import math

class AngleUnit(object):
    """
    A class for defining angular units used by Angle objects.

    **Initialization:**

        An AngleUnit takes a single argument for initialization, a float that specifies the size
        of the desired angular unit in radians.  For example:

            :meth:`coord.AngleUnit.__init__`
 
            >>> gradian = AngleUnit(2. * math.pi / 400.)
            >>> print(gradian)
            coord.AngleUnit(0.015707963267948967)

    **Built-in units:**

        There are five built-in AngleUnits which are always available for use:

            :coord.radians:   coord.AngleUnit(1.)
            :coord.degrees:   coord.AngleUnit(pi / 180.)
            :coord.hours:     coord.AngleUnit(pi / 12.)
            :coord.arcmin:    coord.AngleUnit(pi / 180. / 60.)
            :coord.arcsec:    coord.AngleUnit(pi / 180. / 3600.)

    **Attribute:**

        An AngleUnit as the following (read-only) attribute:

            :value:     The measure of the unit in radians.
    """
    valid_names = ['rad', 'deg', 'hr', 'hour', 'arcmin', 'arcsec']

    def __init__(self, value):
        """
        :param value:   The measure of the unit in radians.
        """
        self._value = float(value)

    @property
    def value(self):
        """A read-only attribute giving the measure of the AngleUnit in radians."""
        return self._value

    def __rmul__(self, theta):
        """float * AngleUnit returns an Angle"""
        from .angle import Angle  # Can't do this at the top, since circular reference.
        return Angle(theta, self)

    def __div__(self, unit):
        """AngleUnit / AngleUnit returns a float giving the relative scaling.

        Note: At least to within machine precision, it is the case that

            (x * angle_unit1) / angle_unit2 == x * (angle_unit1 / angle_unit2)
        """
        if not isinstance(unit, AngleUnit):
            raise TypeError("Cannot divide AngleUnit by %s"%unit)
        return self.value / unit.value

    __truediv__ = __div__

    @staticmethod
    def from_name(unit):
        """Convert a string into the corresponding AngleUnit.

        Only the start of the string is checked, so for instance 'radian' or 'radians' is
        equivalent to 'rad'.

        Valid options are:

            :rad:           AngleUnit(1.)
            :deg:           AngleUnit(pi / 180.)
            :hour or hr:    AngleUnit(pi / 12.)
            :arcmin:        AngleUnit(pi / 180. / 60.)
            :arcsec:        AngleUnit(pi / 180. / 3600.)

        Note: these valid names are listed in AngleUnit.valid_names.

        :param unit:    The string name of the unit to return

        :returns: an AngleUnit
        """
        unit = unit.strip().lower()
        if unit.startswith('rad') :
            return radians
        elif unit.startswith('deg') :
            return degrees
        elif unit.startswith('hour') :
            return hours
        elif unit.startswith('hr') :
            return hours
        elif unit.startswith('arcmin') :
            return arcmin
        elif unit.startswith('arcsec') :
            return arcsec
        else :
            raise ValueError("Unknown Angle unit: %s"%unit)

    def __repr__(self):
        if self == radians:
            return 'coord.radians'
        elif self == degrees:
            return 'coord.degrees'
        elif self == hours:
            return 'coord.hours'
        elif self == arcmin:
            return 'coord.arcmin'
        elif self == arcsec:
            return 'coord.arcsec'
        else:
            return 'coord.AngleUnit(%r)'%self.value

    def __eq__(self, other):
        return isinstance(other,AngleUnit) and self.value == other.value

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(('coord.AngleUnit', self.value))

# Convenient pre-set built-in units
# (These are typically the only ones we will use.)

radians = AngleUnit(1.)
hours = AngleUnit(math.pi / 12.)
degrees = AngleUnit(math.pi / 180.)
arcmin = AngleUnit(math.pi / 10800.)
arcsec = AngleUnit(math.pi / 648000.)
