.. image:: https://github.com/LSSTDESC/Coord/actions/workflows/ci.yml/badge.svg
        :target: https://github.com/LSSTDESC/Coord/actions/workflows/ci.yml
.. image:: https://codecov.io/gh/LSSTDESC/Coord/branch/main/graph/badge.svg
        :target: https://codecov.io/gh/LSSTDESC/Coord

Coord is a Python module that provides basic functionality related to angles and
celestial coordinates.

It is similar in functionality to the astropy.coordinates module, but with more of an
emphasis on efficiency.  Some functions are more than 100 times faster than the corresponding
functionality in astropy.  On the other hand, the API is somewhat more restrictive than
the API used by astropy, so the appropriate module to use will depend on your needs.

Notable functionality:

* Spherical geometric calculations of the distance between two points, angles in spherical
  triangles, and areas of spherical triangles.
* Tangent-plane projection of a coordinate relative to a given center point of the projection.
* Deprojection back into spherical coordinates.
* Analytic calculation of the jacobian of the tangent projection.
* Precession calculations between different epochs.
* Conversion to galactic and ecliptic coordinate systems (and back).

The code was originally written by Mike Jarvis for TreeCorr and then also used in GalSim.
It has been republished here as a stand-alone module to enable others (within the LSST DESC
or not) to have easier access to these functions.

One missing feature (for which pull requests would be welcome) is conversion between FK5 and ICRS
systems (or even FK4).  If you care about the slight differences between these systems, then you
should probably stick to astropy, which does handle these distinctions.

Aside: The module was used as a pedagogical tool at the LSST DESC DE School on July 9, 2017.
A video of the `lesson <http://www.lsst-desc.org/DEschool#MikeJarvis>`_ can be viewed
`here <https://www.youtube.com/watch?v=u3x5OEXgtnU>`_.

Licence
=======

The code is licensed under the MIT License, which basically means you can use it in any way
you want, so long as you keep the copyright notice at the top of each source file and/or include
the accompanying LICENSE file with the distribution.

Installation
============

You can install Coord with pip::

    $ pip install LSSTDESC.Coord --user

or if you have root access, you might prefer::

    $ sudo pip install LSSTDESC.Coord

or possibly with neither ``sudo`` nor ``--user`` if your Python distro is in a writable directory.

If you use anaconda you can install from conda forge::

    $ conda install -c conda-forge lsstdesc.coord

If you prefer to download or clone the repo and install manually, you can install with
setup.py using one of the usual variants::

    $ python setup.py install --prefix={prefix}

or::

    $ sudo python setup.py install
