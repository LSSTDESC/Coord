.. image:: https://api.travis-ci.org/LSSTDESC/Coord.svg?branch=master
        :target: https://travis-ci.org/LSSTDESC/Coord
.. image:: https://codecov.io/gh/LSSTDESC/Coord/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/LSSTDESC/Coord

Coord is a Python module that provides basic functionality related to angles,
celestial coordinates, and spherical geometry calculations.

It is similar in functionality to the astropy.coordinates module, but with more of an
emphasis on efficiency.  Some functions are up to 10 times faster than the corresponding
functionality in astropy.  On the other hand, the API is somewhat more restrictive than
the API used by astropy.

The code was were originally written by Mike Jarvis for TreeCorr and then also used in GalSim.
It has been republished here as a stand-alone module to enable others (within the LSST DESC
or not) to have easier access to these functions.

In addition, the module will be used as a pedagogical tool at the LSST DESC DE School
on July 9, 2017.  (See below for more details about this.)

Licence
=======

The code is licensed under the MIT License, which basically means you can use it in any way
you want, so long as you keep the copyright notice at the top of each source file and/or include
the accompanying LICENSE file with the distribution.

Installation
============

This module is mostly pure Python.  However, it does have a small bit of C++ code which is
wrapped with cffi, which in turn depends on libffi.  Most modern Python installations will have
this already installed, so you probably don't have to do anything special to use it.  However,
if not, installing libffi is fairly straightforward.

Installing libffi on a Mac
--------------------------

You should only need to run the command::

    $ xcode-select --install

Installing libffi on Linux
--------------------------

If you have root access, then one of the following should work for you::

    $ apt-get install libffi-dev
    $ yum install libffi-devel

Installing libffi from source
-----------------------------

If you don't have root access (and don't want to bother your sysadmin), then installing from
source into a {prefix} directory where you have write access (e.g. your home directory) is also
not very hard::

    $ wget ftp://sourceware.org:/pub/libffi/libffi-3.2.1.tar.gz
    $ tar xfz libffi-3.2.1.tar.gz
    $ cd libffi-3.2.1
    $ ./configure --prefix={prefix}
    $ make
    $ make install
    $ cp */include/ffi*.h {prefix}/include
    $ cd ..

Installing Coord
----------------

Once you have done one of the above (or ignored that in the hope that you already have libffi
installed), you can install Coord simply with::

    $ python setup.py install --prefix={prefix}

or if you have root access, you might prefer::

    $ sudo python setup.py install

Or any of the other usual variants of this command as you prefer.


DESC DE School
==============

For the DE School lesson on July 10, 2017, we will be designing some tests for this module, so
I havve intentionally not ported over any of the existing tests that I have in the TreeCorr and
GalSim repositories related to that version of the code.  If you are planning to participate in
the DE School lesson, then you are on your honor not to look at those tests.

Prior to attending DE School, it would be a good idea to try to clone this repository and install
the code.  If you have trouble, please ask on the #desc-bnl-sb-2017 Slack channel.  It will be
helpful if you have the code already installed on your laptop, although you may certainly still
participate even if not.

After the end of the class I will add all the unit tests that I currently have for this module
along with (hopefully) some that the class will have written.

