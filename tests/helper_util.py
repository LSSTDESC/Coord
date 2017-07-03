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

try:
    import cPickle as pickle
except ImportError:
    import pickle
import copy
from collections import Hashable, Counter
import functools
import time
import coord

# This file has some helper functions that are used by tests from multiple files to help
# avoid code duplication.

# The functions in this file are based heavily on the corresponding functions in the
# galsim_test_helpers in the GalSim/tests directory.

def do_pickle(obj, func=(lambda x : x), copyable=True, reprable=True):
    """Check that the object is picklable.  Also that it has basic == and != functionality.

    :param obj:         The object to be checked
    :param func:        If desired, some function of the object to use for equality checking.
                        [default: (lambda x : x)]
    :param copyable:    Whether the object should be expected to be equal after copy(obj).
                        e.g. things with an internal random state may not be equal after copying.
                        [default: True]
    :param reprable:    Whether the object should be expected to survive a round trip through
                        obj = eval(repr(obj)).  This is often a desirable property of objects,
                        but it may not always be possible (or efficient). [default: True]
    """
    # Test pickle round-trip returns an equivalent (but not identical) object
    print('Try pickling ',obj)  # Note: this implicitly checks that str(obj) works.
    obj2 = pickle.loads(pickle.dumps(obj))
    assert obj2 is not obj, "Obj %s not identical after pickling"%obj
    f1 = func(obj)
    f2 = func(obj2)
    assert f2 == f1, "func(obj) = %r\nfunc(obj2) = %r"%(f1, f2)

    # Test the hash values are equal for two equivalent objects.
    if isinstance(obj, Hashable):
        assert hash(obj) == hash(obj2), "Objects %s and %s have different hashes"%(obj,obj2)

    # Test that copy makes an equivalent copy.
    obj3 = copy.copy(obj)
    assert obj3 is not obj, "Obj %s not identical after copy"%obj
    if copyable:
        f1 = func(obj)  # Might need to redo this in case function changes an internal state.
        f3 = func(obj3)
        assert f3 == f1, "func(obj) = %r\nfunc(obj3) = %r"%(f1, f3)

    # A deepcopy should always work.
    obj4 = copy.deepcopy(obj)
    assert obj4 is not obj, "Obj %s not identical after deep copy"%obj
    f1 = func(obj)
    f4 = func(obj4)
    assert f4 == f1, "func(obj) = %r\nfunc(obj4) = %r"%(f1, f4)

    # Also test that the repr is an accurate representation of the object.
    # The gold standard is that eval(repr(obj)) == obj.  So check that here as well.

    if reprable:
        obj5 = eval(repr(obj))
        f1 = func(obj)
        f5 = func(obj5)
        assert f5 == f1, "func(obj) = %r\nfunc(obj5) = %r"%(f1, f5)
    else:
        # Even if we're not actually doing the test, still make the repr to check for syntax errors.
        repr(obj)


def all_obj_diff(objs):
    """Helper function that verifies that each element in `objs` is unique and, if hashable,
    produces a unique hash.
    """
    # Check that all objects are unique.
    # Would like to use `assert len(objs) == len(set(objs))` here, but this requires that the
    # elements of objs are hashable (and that they have unique hashes!, which is what we're trying
    # to test!.  So instead, we just loop over all combinations.
    for i, obji in enumerate(objs):
        # Could probably start the next loop at `i+1`, but we start at 0 for completeness
        # (and to verify a != b implies b != a)
        for j, objj in enumerate(objs):
            if i == j:
                continue
            assert obji != objj, ("Found equivalent objects {0} == {1} at indices {2} and {3}"
                                  .format(obji, objj, i, j))

    # Now check that all hashes are unique (if the items are hashable).
    if not isinstance(objs[0], Hashable):
        return
    hashes = [hash(obj) for obj in objs]
    try:
        assert len(hashes) == len(set(hashes)), "At least two of the hashes are equal"
    except AssertionError as e:
        for k, v in Counter(hashes).items():
            if v <= 1:
                continue
            print("Found multiple equivalent object hashes:")
            for i, obj in enumerate(objs):
                if hash(obj) == k:
                    print(i, repr(obj))
        raise e


def timer(f):
    """A decorator to put before each test function to have it output how long it took to run the
    test.
    """
    @functools.wraps(f)
    def f2(*args, **kwargs):
        t0 = time.time()
        result = f(*args, **kwargs)
        t1 = time.time()
        fname = repr(f).split()[1]
        print('time for %s = %.2f' % (fname, t1-t0))
        return result
    return f2
