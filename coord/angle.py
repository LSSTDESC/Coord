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

"""@file angle.py
A class for handling Angles with proper angular units
"""

import math
import numpy as np

from .angleunit import AngleUnit, radians, degrees, hours

class Angle(object):
    """A class representing an Angle.

    Initialization
    --------------

    Angles are a value with an AngleUnit.

    You typically create an Angle by multiplying a number by a coord.AngleUnit, for example::

        >>> pixel = 0.27 * coord.arcsec
        >>> ra = 13.4 * coord.hours
        >>> dec = -32 * coord.degrees
        >>> from math import pi
        >>> theta = pi/2. * coord.radians

    You can also initialize explicitly, taking a value and a unit::

        >>> phi = coord.Angle(90, coord.degrees)

    There are five built-in AngleUnits which are always available for use:

        :coord.radians:   # = coord.AngleUnit(1.)
        :coord.degrees:   # = coord.AngleUnit(pi / 180.)
        :coord.hours:     # = coord.AngleUnit(pi / 12.)
        :coord.arcmin:    # = coord.AngleUnit(pi / 180. / 60.)
        :coord.arcsec:    # = coord.AngleUnit(pi / 180. / 3600.)

    Radian access method
    --------------------

    Since extracting the value in radians is extremely common, we have an accessor method to do this
    quickly::

        >>> x = theta.rad()
        >>> print x
        1.57079632679

    It is equivalent to the more verbose::

        >>> x = theta / coord.radians

    but without actually requiring the floating point operation of dividing by 1.

    Operations
    ----------

    Allowed arithmetic with Angles include the following:
    (In the list below, `x` is a float, `unit` is a coord.AngleUnit, `theta` is a coord.Angle)::

        >>> theta = x * unit
        >>> x = theta / unit
        >>> theta3 = theta1 + theta2
        >>> theta3 = theta1 - theta2
        >>> theta2 = theta1 * x
        >>> theta2 = x * theta1
        >>> theta2 = theta1 / x
        >>> theta2 = -theta1
        >>> theta2 += theta1
        >>> theta2 -= theta1
        >>> theta *= x
        >>> theta /= x
        >>> x = unit1 / unit2   # equivalent to x = (1 * unit1) / unit2

    Operations on NumPy arrays containing Angles are permitted, provided that they are within the
    bounds of the allowed operations on Angles listed above (e.g., addition/subtraction of Angles,
    multiplication of an Angle by a float, but not multiplication of Angles together).

    There are convenience function for getting the sin, cos, and tan of an angle, along with
    one for getting sin and cos together, which should be more efficient than doing sin and
    cos separately:

        >>> sint = theta.sin()  # equivalent to sint = math.sin(theta.rad())
        >>> cost = theta.cos()  # equivalent to cost = math.cos(theta.rad())
        >>> tant = theta.tan()  # equivalent to tant = math.tan(theta.rad())
        >>> sint, cost = theta.sincos()

    Wrapping
    --------

    Depending on the context, theta = 2pi radians and theta = 0 radians are the same thing.
    If you want your angles to be wrapped to [-pi,pi) radians, you can do this by calling

        >>> theta = theta.wrap()

    This could be appropriate before testing for the equality of two angles for example, or
    calculating the difference between them.
    """
    def __init__(self, theta, unit=None):
        """
        :param theta:   The numerical value of the angle.
        :param unit:    The units theta is measured in.
        """
        # We also want to allow angle1 = Angle(angle2) as a copy, so check for that.
        if isinstance(theta,Angle):
            if unit is not None:
                raise TypeError("Cannot provide unit if theta is already an Angle instance")
            self._rad = theta._rad
        elif unit is None:
            raise TypeError("Must provide unit for Angle.__init__")
        elif not isinstance(unit, AngleUnit):
            raise TypeError("Invalid unit %s of type %s"%(unit,type(unit)))
        else:
            # Normal case
            self._rad = float(theta) * unit.value

    def rad(self):
        """Return the Angle in radians.

        Equivalent to angle / coord.radians
        """
        return self._rad

    def __neg__(self):
        return _Angle(-self._rad)

    def __add__(self, other):
        if not isinstance(other, Angle):
            raise TypeError("Cannot add %s of type %s to an Angle"%(other,type(other)))
        return _Angle(self._rad + other._rad)

    def __sub__(self, other):
        if not isinstance(other, Angle):
            raise TypeError("Cannot subtract %s of type %s from an Angle"%(other,type(other)))
        return _Angle(self._rad - other._rad)

    def __mul__(self, other):
        if other != float(other):
            raise TypeError("Cannot multiply Angle by %s of type %s"%(other,type(other)))
        return _Angle(self._rad * other)

    __rmul__ = __mul__

    def __div__(self, other):
        if isinstance(other, AngleUnit):
            return self._rad / other.value
        elif other == float(other):
            return _Angle(self._rad / other)
        else:
            raise TypeError("Cannot divide Angle by %s of type %s"%(other,type(other)))

    __truediv__ = __div__

    def wrap(self, center=None):
        """Wrap Angle to lie in the range [-pi, pi) radians.

        Depending on the context, theta = 2pi radians and theta = 0 radians are the same thing.
        If you want your angles to be wrapped to [-pi, pi) radians, you can do this by calling

            >>> theta = theta.wrap()

        This could be appropriate before testing for the equality of two angles for example, or
        calculating the difference between them.

        If you want to wrap to a different range than [-pi, pi), you can set the `center` argument
        to be the desired center of the the range.  e.g. for return values to fall in [0, 2pi),
        you could call

            >>> theta = theta.wrap(center=180. * coord.degrees)

        @param center   The center point of the wrapped range. [default: 0 radians]

        @returns the equivalent angle within the range [center-pi, center+pi)
        """
        if center is None: center = _Angle(0.)
        start = center._rad - math.pi
        offset = (self._rad - start) // (2.*math.pi)  # How many full cycles to subtract
        return _Angle(self._rad - offset * 2.*math.pi)

    def sin(self):
        """Return the sin of an Angle."""
        return math.sin(self._rad)

    def cos(self):
        """Return the cos of an Angle."""
        return math.cos(self._rad)

    def tan(self):
        """Return the tan of an Angle."""
        return math.tan(self._rad)

    def sincos(self):
        """Return both the sin and cos of an Angle as a numpy array [sint, cost].

        (On some systems, this may be slightly faster than doing each separately.)
        """
        # Note: This is admittedly pretty gratuitous.  We're writing Python after all, so the
        # difference between this and two trig calls is probably neglible given all the normal
        # Python overhead.  Mostly, I added this as an excuse to set up the C++ extension stuff
        # to make sure people could see how to add further C++-layer optimizations as needed.
        import coord
        sc = np.empty(2)
        coord._lib.coord_sincos(self._rad, coord._ffi.cast('double*', sc.ctypes.data))
        return sc

    def __str__(self):
        return str(self._rad) + ' radians'

    def __repr__(self):
        return 'coord.Angle(%r, coord.radians)'%self.rad()

    def __eq__(self, other):
        return isinstance(other,Angle) and self.rad() == other.rad()

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(('coord.Angle', self._rad))

    @staticmethod
    def _make_dms_string(decimal, sep):
        if decimal >= 0:
            sign = '+'
        else:
            sign = '-'
            decimal = -decimal
        # Round to nearest 1.e-8 seconds
        decimal = int(3600.e8 * decimal + 0.5)
        d = decimal // 360000000000
        decimal -= d * 360000000000
        m = decimal // 6000000000
        decimal -= m * 6000000000
        s = decimal // 100000000
        decimal -= s * 100000000
        return '%s%02d%s%02d%s%02d.%08d'%(sign,d,sep,m,sep,s,decimal)

    def hms(self, sep=":"):
        """Return an HMS representation of the angle as a string: hh:mm:ss.decimal.

        The returned representation will have 0 <= hh < 24.

        An optional `sep` parameter can change the : to something else (e.g. a space or
        nothing at all).

        Note: the reverse process is effected by HMS_Angle:

            >>> angle = -5.357 * coord.hours
            >>> hms = angle.hms()
            >>> print hms
            +18:38:34.80000000
            >>> angle2 = coord.HMS_Angle(hms)
            >>> print angle2 / coord.hours
            18.643
            >>> print angle2 / coord.hours - 24
            -5.357
            >>> print angle2 - angle - 24 * coord.hours
            0.0 radians

        @param sep      The token to put between the hh and mm, and beteen mm and ss. [default: ':']

        @returns a string of the HMS representation of the angle.
        """
        # HMS convention is usually to have the hours between 0 and 24, not -12 and 12
        h = self.wrap(12. * hours) / hours
        return self._make_dms_string(h,sep)

    def dms(self, sep=":"):
        """Return a DMS representation of the angle as a string: (+/-)ddmmss.decimal
        An optional `sep` parameter can change the : to something else (e.g. a space or
        nothing at all).

        Note: the reverse process is effected by DMS_Angle:

            >>> angle = -(5 * coord.degrees + 13 * coord.arcmin + 23 * coord.arcsec)
            >>> dms = angle.dms()
            >>> print dms
            -05:13:23.00000000
            >>> angle2 = coord.DMS_Angle(dms)
            >>> print angle2 / coord.degrees
            -5.22305555556
            >>> print angle2 - angle
            0.0 radians

        @param sep      The token to put between the dd and mm, and beteen mm and ss. [default: ':']

        @returns a string of the DMS representation of the angle.
        """
        d = self.wrap() / degrees
        return self._make_dms_string(d,sep)

    @staticmethod
    def parse_dms(dms):
        """Convert a string of the form dd:mm:ss.decimal into decimal degrees."""
        sign = 1
        if dms[0] == '-':
            sign = -1
            dms = dms[1:]

        d, m, s = dms.split(':')

        return sign * (int(d) + int(m)/60. + float(s)/3600.)


def _Angle(theta):
    """Equivalent to either `theta * coord.radians` or `Angle(theta, coord.radians)`, but without
    the normal overhead (which isn't much to be honest, but this is nonetheless slightly quicker).
    """
    ret = Angle.__new__(Angle)
    ret._rad = theta
    return ret


def HMS_Angle(str):
    """Convert a string of the form hh:mm:ss.decimal into an Angle.

    There may be an initial + or - (or neither), then two digits for the hours, two for the
    minutes, and two for the seconds.  Then there may be a decimal point followed by more
    digits.  There may be a colon separating hh, mm, and ss, or whitespace, or nothing at all.
    In fact, the code will ignore any non-digits between the hours, minutes, and seconds.

    Note: the reverse process is effected by Angle.hms():

        >>> angle = -5.357 * coord.hours
        >>> hms = angle.hms()
        >>> print hms
        +18:38:34.80000000
        >>> angle2 = coord.HMS_Angle(hms)
        >>> print angle2 / coord.hours
        18.643
        >>> print angle2 / coord.hours - 24
        -5.357
        >>> print angle2 - angle - 24 * coord.hours
        0.0 radians

    @returns the corresponding Angle instance
    """
    return Angle.parse_dms(str) * hours

def DMS_Angle(str):
    """Convert a string of the form dd:mm:ss.decimal into an Angle.

    There may be an initial + or - (or neither), then two digits for the degrees, two for the
    minutes, and two for the seconds.  Then there may be a decimal point followed by more
    digits.  There may be a colon separating dd, mm, and ss, or whitespace, or nothing at all.
    In fact, the code will ignore any non-digits between the degrees, minutes, and seconds.

    @returns the corresponding Angle instance
    """
    return Angle.parse_dms(str) * degrees

