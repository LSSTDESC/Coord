DESC DE School
==============

This module was used as a pedagogical tool for an LSST DESC DE School lesson on July 10, 2017.

A video of the `lesson <http://www.lsst-desc.org/DEschool#MikeJarvis>`_ can be viewed
`here <https://www.youtube.com/watch?v=u3x5OEXgtnU>`_.

Prior to the lesson, only a few of the unit tests had been copied over to this repository
from the TreeCorr and/or GalSim repositories, where earlier versions of the code lived.
The lesson culminated in students having an opportunity to practice designing unit tests by
filling in some of the missing tests.

Original instructions for DE School stutents:
---------------------------------------------

Prior to attending DE School, it would be a good idea to try to clone this repository and install
the code.  If you have trouble, please ask for help on the #desc-bnl-sb-2017 Slack channel.
It will be helpful if you have the code already installed on your laptop, although you may
certainly still participate even if not.

After the end of the class I will add all the unit tests that I currently have for this module
along with (hopefully) some that the class will have written.

1. Clone the repository onto your laptop.::

    $ git clone git@github.com:LSSTDESC/Coord.git

2. Follow the above instructions to install Coord.

3. Make sure you have nose and astropy installed.::

    $ pip install nose astropy

4. Try running the existing tests.::

    $ cd tests
    $ nosetests

5. Make sure that you are a member of the LSSTDESC GitHub team.  If not, post a request on the
   #desc-ci Slack channel that you want to be added.

6. Look through the documentation here:

    https://lsstdesc.github.io/Coord/

7. Make a branch for yourself.  Pick some appropriate name for your branch, like your last name,
   or your initials.  Something you think would likely be unique to you.::

    $ git pull
    $ git checkout -b your_branch_name
    $ git push -u origin your_branch_name

8. Look through the existing tests and pick a test that is currently unimplemented.
   There are lots of stubs of test functions that you can pick from that currently just say
   ``pass``, or feel free to do something different if you have an idea of what to test and
   don't seem something appropriate.

   Write the test.  Think about the various considerations I talked about in the lesson.

   Be as complete or not as you want.  You only have half an hour, so don't worry if you don't
   feel like you covered everything about the concept you are trying to test.

9. One the test is working, commit and push.::

    $ git add -p
    $ git commit -m "Add some tests of the .... function"
    $ git push

10. Make a pull request of your new test.  Start here:

     https://github.com/LSSTDESC/Coord/compare

    Select the "compare" tab, and find your branch name.

    Write a short title.  e.g. "Add test of stereographic projection"

    Write a slightly longer summary in the text box.  (Normally -- for this exercise there might
    not be much more to write than what you already wrote in the title.)

    Click "Create pull request"

    Wait for Travis to run your test.  If it fails, check why and edit your code appropriately.
    You can add more commits to your branch, which will be included in the pull request.
