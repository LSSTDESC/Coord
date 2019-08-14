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

import math
import numpy as np

from .angleunit import AngleUnit, radians, degrees, hours, arcmin, arcsec

class Angle(object):
    """A class representing an Angle.  Angles are a value with an AngleUnit.

    **Initialization:**

        You typically create an Angle by multiplying a number by a coord.AngleUnit, for example:

            ..

            >>> pixel = 0.27 * arcsec
            >>> ra = 13.4 * hours
            >>> dec = -32 * degrees
            >>> from math import pi
            >>> theta = pi/2. * radians

        You can also initialize explicitly, taking a value and a unit:

            :meth:`coord.Angle.__init__`

            >>> unit = AngleUnit(math.pi / 100)  # gradians
            >>> phi = Angle(90, unit)

    **Built-in units:**

        There are five built-in AngleUnits which are always available for use:

            :coord.radians:   coord.AngleUnit(1.)
            :coord.degrees:   coord.AngleUnit(pi / 180.)
            :coord.hours:     coord.AngleUnit(pi / 12.)
            :coord.arcmin:    coord.AngleUnit(pi / 180. / 60.)
            :coord.arcsec:    coord.AngleUnit(pi / 180. / 3600.)

    **Attribute:**

        Since extracting the value in radians is extremely common, we have a read-only attribute
        to do this quickly:

            :rad:       The measure of the unit in radians.

        For example:

            ..

            >>> theta = 90 * degrees
            >>> print(theta.rad)
            1.5707963267948966

        It is equivalent to the more verbose:

            ..

            >>> x = theta / radians
            >>> print(x)
            1.5707963267948966

        but without actually requiring the floating point operation of dividing by 1.

    **Arithmetic:**

        Allowed arithmetic with Angles include the following.
        In the list below,

            - ``x`` is an arbitrary ``float`` value
            - ``unit1`` and ``unit2`` are arbitrary `AngleUnit` instances
            - ``theta1`` and ``theta2`` are arbitrary `Angle` instances

            >>> x = 37.8
            >>> unit1 = arcmin
            >>> unit2 = degrees

            >>> theta1 = x * unit1
            >>> theta2 = x * unit2
            >>> x2 = theta1 / unit2
            >>> theta = theta1 + theta2
            >>> theta = theta1 - theta2
            >>> theta = theta1 * x
            >>> theta = x * theta1
            >>> theta = theta1 / x
            >>> theta = -theta1
            >>> theta += theta1
            >>> theta -= theta1
            >>> theta *= x
            >>> theta /= x
            >>> x = unit1 / unit2   # equivalent to x = (1 * unit1) / unit2

        The above operations on NumPy arrays containing Angles are permitted as well.

    **Trigonometry:**

        There are convenience function for getting the sin, cos, and tan of an angle, along with
        one for getting sin and cos together, which should be more efficient than doing sin and
        cos separately:

            | :meth:`coord.Angle.sin`
            | :meth:`coord.Angle.cos`
            | :meth:`coord.Angle.tan`
            | :meth:`coord.Angle.sincos`

            >>> sint = theta.sin()  # equivalent to sint = math.sin(theta.rad)
            >>> cost = theta.cos()  # equivalent to cost = math.cos(theta.rad)
            >>> tant = theta.tan()  # equivalent to tant = math.tan(theta.rad)
            >>> sint, cost = theta.sincos()

        These functions mean that numpy trig functions will work on Angles or arrays of Angles:

            ..

            >>> sint = np.sin(theta)
            >>> cost = np.cos(theta)
            >>> tant = np.tan(theta)

    **Wrapping:**

        Depending on the context, theta = 2pi radians and theta = 0 radians may mean the same thing.
        If you want your angles to be wrapped to [-pi,pi) radians, you can do this by calling

            :meth:`coord.Angle.wrap`

            >>> theta = theta.wrap()

        This could be appropriate before testing for the equality of two angles for example, or
        calculating the difference between them.

        There is also an option to wrap into a different 2 pi range if so desired by specifying
        the center of the range.
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

    @property
    def rad(self):
        """Return the Angle in radians.

        Equivalent to angle / coord.radians
        """
        return self._rad

    @property
    def deg(self):
        """Return the Angle in degrees.

        Equivalent to angle / coord.degrees
        """
        return self / degrees

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
        """Wrap Angle to lie in the range [-pi, pi) radians (or other range  of 2pi radians)

        Depending on the context, theta = 2pi radians and theta = 0 radians are the same thing.
        If you want your angles to be wrapped to [-pi, pi) radians, you can do this as follows:

            ..

            >>> theta = Angle(700 * degrees)
            >>> theta = theta.wrap()
            >>> print(theta.deg)
            -19.99999999999998

        This could be appropriate before testing for the equality of two angles for example, or
        calculating the difference between them.

        If you want to wrap to a different range than [-pi, pi), you can set the ``center`` argument
        to be the desired center of the the range.  e.g. for return values to fall in [0, 2pi),
        you could call

            ..

            >>> theta = theta.wrap(center=180. * degrees)
            >>> print(theta / degrees)
            340.0

        :param center:  The center point of the wrapped range. [default: 0 radians]

        :returns: the equivalent angle within the range [center-pi, center+pi)
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
        return 'coord.Angle(%r, coord.radians)'%self.rad

    def __eq__(self, other):
        return isinstance(other,Angle) and self.rad == other.rad

    def __ne__(self, other):
        return not self.__eq__(other)

    def __le__(self, other):
        if not isinstance(other, Angle):
            raise TypeError("Cannot compare %s of type %s to an Angle"%(other,type(other)))
        return self._rad <= other._rad

    def __lt__(self, other):
        if not isinstance(other, Angle):
            raise TypeError("Cannot compare %s of type %s to an Angle"%(other,type(other)))
        return self._rad < other._rad

    def __ge__(self, other):
        if not isinstance(other, Angle):
            raise TypeError("Cannot compare %s of type %s to an Angle"%(other,type(other)))
        return self._rad >= other._rad

    def __gt__(self, other):
        if not isinstance(other, Angle):
            raise TypeError("Cannot compare %s of type %s to an Angle"%(other,type(other)))
        return self._rad > other._rad

    def __hash__(self):
        return hash(('coord.Angle', self._rad))

    @staticmethod
    def _make_dms_string(decimal, sep, prec, pad, plus_sign):
        # Account for the sign properly
        if decimal < 0:
            sign = '-'
            decimal = -decimal
        elif plus_sign:
            sign = '+'
        else:
            sign = ''

        # Figure out the 3 sep tokens
        sep1 = sep2 = ''
        sep3 = None
        if len(sep) == 1:
            sep1 = sep2 = sep
        elif len(sep) == 2:
            sep1, sep2 = sep
        elif len(sep) == 3:
            sep1, sep2, sep3 = sep

        # Round to nearest 1.e-8 seconds (or 10**-prec if given)
        round_prec = 8 if prec is None else prec
        digits = 10**round_prec

        decimal = int(3600 * digits * decimal + 0.5)

        d = decimal // (3600 * digits)
        decimal -= d * (3600 * digits)
        m = decimal // (60 * digits)
        decimal -= m * (60 * digits)
        s = decimal // digits
        decimal -= s * digits

        # Make the string
        if pad:
            d_str = '%02d'%d
            m_str = '%02d'%m
            s_str = '%02d'%s
        else:
            d_str = '%d'%d
            m_str = '%d'%m
            s_str = '%d'%s
        string = '%s%s%s%s%s%s.%0*d'%(sign,d_str,sep1,m_str,sep2,s_str,round_prec,decimal)
        if not prec:
            string = string.rstrip('0')
            string = string.rstrip('.')
        if sep3:
            string = string + sep3
        return string

    def hms(self, sep=":", prec=None, pad=True, plus_sign=False):
        """Return an HMS representation of the angle as a string: +-hh:mm:ss.decimal.

        An optional ``sep`` parameter can change the : to something else (e.g. a space or
        nothing at all).

        Note: the reverse process is effected by :meth:`Angle.from_hms`:

            ..

            >>> angle = -5.357 * hours
            >>> hms = angle.hms()
            >>> print(hms)
            -05:21:25.2
            >>> angle2 = Angle.from_hms(hms)
            >>> print(angle2 / hours)
            -5.356999999999999

        :param sep:         The token to put between the hh and mm and beteen mm and ss.  This may
                            also be a string of 2 or 3 items, e.g. 'hm' or 'hms'.  Or even a
                            tuple of strings such as ('hours ', 'minutes ', 'seconds').
                            [default: ':']
        :param prec:        The number of digits of precision after the decimal point.
                            [default: None]
        :param pad:         Whether to pad with a leading 0 if necessary to make h,m,s 2 digits.
                            [default: True]
        :param plus_sign:   Whether to use a plus sign for positive angles. [default: False]

        :returns: a string of the HMS representation of the angle.
        """
        if not len(sep) <= 3:
            raise ValueError("sep must be a string or tuple of length <= 3")
        if prec is not None and not prec >= 0:
            raise ValueError("prec must be >= 0")
        return self._make_dms_string(self/hours, sep, prec, pad, plus_sign)

    def dms(self, sep=":", prec=None, pad=True, plus_sign=False):
        """Return a DMS representation of the angle as a string: +-dd:mm:ss.decimal
        An optional ``sep`` parameter can change the : to something else (e.g. a space or
        nothing at all).

        Note: the reverse process is effected by :meth:`Angle.from_dms`:

            ..

            >>> angle = -(5 * degrees + 21 * arcmin + 25.2 * arcsec)
            >>> dms = angle.dms()
            >>> print(dms)
            -05:21:25.2
            >>> angle2 = Angle.from_dms(dms)
            >>> print(angle2 / degrees)
            -5.356999999999999

        :param sep:         The token to put between the hh and mm and beteen mm and ss.  This may
                            also be a string of 2 or 3 items, e.g. 'dm' or 'dms'.  Or even a
                            tuple of strings such as ('degrees ', 'minutes ', 'seconds').
                            [default: ':']
        :param prec:        The number of digits of precision after the decimal point.
                            [default: None]
        :param pad:         Whether to pad with a leading 0 if necessary to make h 2 digits.
                            [default: True]
        :param plus_sign:   Whether to use a plus sign for positive angles. [default: False]

        :returns: a string of the DMS representation of the angle.
        """
        if not len(sep) <= 3:
            raise ValueError("sep must be a string or tuple of length <= 3")
        if prec is not None and not prec >= 0:
            raise ValueError("prec must be >= 0")
        return self._make_dms_string(self/degrees, sep, prec, pad, plus_sign)

    @staticmethod
    def from_hms(str):
        """Convert a string of the form hh:mm:ss.decimal into an Angle.

        There may be an initial + or - (or neither), then two digits for the hours, two for the
        minutes, and two for the seconds.  Then there may be a decimal point followed by more
        digits.  There may be a colon separating hh, mm, and ss, or whitespace, or nothing at all.
        In fact, the code will ignore any non-digits between the hours, minutes, and seconds.

        Note: the reverse process is effected by Angle.hms():

            ..

            >>> angle = -5.357 * hours
            >>> hms = angle.hms()
            >>> print(hms)
            -05:21:25.2
            >>> angle2 = Angle.from_hms(hms)
            >>> print(angle2 / hours)
            -5.356999999999999

        :param str:     The string to parse.

        :returns: the corresponding Angle instance
        """
        return Angle._parse_dms(str) * hours

    @staticmethod
    def from_dms(str):
        """Convert a string of the form dd:mm:ss.decimal into an Angle.

        There may be an initial + or - (or neither), then two digits for the degrees, two for the
        minutes, and two for the seconds.  Then there may be a decimal point followed by more
        digits.  There may be a colon separating dd, mm, and ss, or whitespace, or nothing at all.
        In fact, the code will ignore any non-digits between the degrees, minutes, and seconds.

        Note: the reverse process is effected by Angle.dms():

            ..

            >>> angle = -(5 * degrees + 21 * arcmin + 25.2 * arcsec)
            >>> dms = angle.dms()
            >>> print(dms)
            -05:21:25.2
            >>> angle2 = Angle.from_dms(dms)
            >>> print(angle2 / degrees)
            -5.356999999999999

        :param str:     The string to parse.

        :returns: the corresponding Angle instance
        """
        return Angle._parse_dms(str) * degrees

    @staticmethod
    def _parse_dms(dms):
        """Convert a string of the form dd:mm:ss.decimal into decimal degrees.
        """
        import re
        tokens = tuple(filter(None, re.split('([\.\d]+)', dms.strip())))
        if len(tokens) <= 1:
            raise ValueError("string is not of the expected format")
        sign = 1
        try:
            dd = float(tokens[0])
        except ValueError:
            if tokens[0].strip() == '-':
                sign = -1
            tokens = tokens[1:]
            dd = float(tokens[0])
            if len(tokens) <= 1:
                raise ValueError("string is not of the expected format")
        if len(tokens) <= 2:
            return sign * dd
        mm = float(tokens[2])
        if len(tokens) <= 4:
            return sign * (dd + mm/60)
        if len(tokens) >= 7:
            raise ValueError("string is not of the expected format")
        ss = float(tokens[4])
        return sign * (dd + mm/60. + ss/3600.)

def _Angle(theta):
    """Equivalent to ``Angle(theta, coord.radians)``, but without the normal overhead (which isn't
    much to be honest, but this is nonetheless slightly quicker).

    :param theta:   The numerical value of the angle in radians.
    """
    ret = Angle.__new__(Angle)
    ret._rad = theta
    return ret

