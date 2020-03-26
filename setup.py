from __future__ import print_function
import sys,os,glob,re

import setuptools
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from setuptools.command.easy_install import easy_install

print('Python version = ',sys.version)
py_version = "%d.%d"%sys.version_info[0:2]  # we check things based on the major.minor version.

with open('requirements.txt') as f:
    required = f.read().splitlines()

sources = glob.glob(os.path.join('src','*.cpp'))
print('sources = ',sources)
headers = glob.glob(os.path.join('include','*.h'))
print('headers = ',headers)

copt =  {
    'gcc' : ['-O3','-ffast-math'],
    'icc' : ['-O3'],
    'clang' : ['-O3','-ffast-math', '-stdlib=libc++'],
    'unknown' : [],
}
lopt =  {
    'gcc' : [],
    'icc' : [],
    'clang' : ['-stdlib=libc++'],
    'unknown' : [],
}

undef_macros = []

# If we build with debug, also undefine NDEBUG flag
if "--debug" in sys.argv:
    undef_macros+=['NDEBUG']
    # Usually already there, but make sure -g in included if we are debugging.
    for name in copt.keys():
        if name != 'unknown':
            copt[name].append('-g')
    debug = True

local_tmp = 'tmp'

def get_compiler_type(compiler, check_unknown=True, output=False):
    """Try to figure out which kind of compiler this really is.
    In particular, try to distinguish between clang and gcc, either of which may
    be called cc or gcc.
    """
    if debug: output=True
    cc = compiler.compiler_so[0]
    if cc == 'ccache':
        cc = compiler.compiler_so[1]
    cmd = [cc,'--version']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    lines = p.stdout.readlines()
    if output:
        print('compiler version information: ')
        for line in lines:
            print(line.decode().strip())
    # Python3 needs this decode bit.
    # Python2.7 doesn't need it, but it works fine.
    line = lines[0].decode(encoding='UTF-8')
    if line.startswith('Configured'):
        line = lines[1].decode(encoding='UTF-8')

    if 'clang' in line:
        return 'clang'
    elif 'gcc' in line:
        return 'gcc'
    elif 'GCC' in line:
        return 'gcc'
    elif 'clang' in cc:
        return 'clang'
    elif 'gcc' in cc or 'g++' in cc:
        return 'gcc'
    elif 'icc' in cc or 'icpc' in cc:
        return 'icc'
    elif check_unknown:
        if output:
            print('Unknown compiler.')
        for cc_type in ['gcc', 'clang']:
            if output:
                print('Check if the compiler works like ',cc_type)
            if try_cc(compiler, cc_type):
                return cc_type
        # I guess none of them worked.  Now we really do have to bail.
        if output:
            print("None of these compile options worked.  Not adding any optimization flags.")
        return 'unknown'
    else:
        return 'unknown'

def try_compile(cpp_code, compiler, cflags=[], lflags=[], prepend=None):
    """Check if compiling some code with the given compiler and flags works properly.
    """
    # Put the temporary files in a local tmp directory, so that they stick around after failures.
    if not os.path.exists(local_tmp): os.makedirs(local_tmp)

    # We delete these manually if successful.  Otherwise, we leave them in the tmp directory
    # so the user can troubleshoot the problem if they were expecting it to work.
    with tempfile.NamedTemporaryFile(delete=False, suffix='.cpp', dir=local_tmp) as cpp_file:
        cpp_file.write(cpp_code.encode())
        cpp_name = cpp_file.name

    # Just get a named temporary file to write to:
    with tempfile.NamedTemporaryFile(delete=False, suffix='.os', dir=local_tmp) as o_file:
        o_name = o_file.name

    # Another named temporary file for the executable
    with tempfile.NamedTemporaryFile(delete=False, suffix='.exe', dir=local_tmp) as exe_file:
        exe_name = exe_file.name

    # Try compiling with the given flags
    cc = [compiler.compiler_so[0]]
    if prepend:
        cc = [prepend] + cc
    cmd = cc + compiler.compiler_so[1:] + cflags + ['-c',cpp_name,'-o',o_name]
    if debug:
        print('cmd = ',' '.join(cmd))
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        lines = p.stdout.readlines()
        p.communicate()
        if debug and p.returncode != 0:
            print('Trying compile command:')
            print(' '.join(cmd))
            print('Output was:')
            print('   ',b'   '.join(lines).decode())
        returncode = p.returncode
    except (IOError,OSError) as e:
        if debug:
            print('Trying compile command:')
            print(cmd)
            print('Caught error: ',repr(e))
        returncode = 1
    if returncode != 0:
        # Don't delete files in case helpful for troubleshooting.
        return False

    # Link
    cc = compiler.linker_so[0]
    cmd = [cc] + compiler.linker_so[1:] + lflags + [o_name,'-o',exe_name]
    if debug:
        print('cmd = ',' '.join(cmd))
    try:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        lines = p.stdout.readlines()
        p.communicate()
        if debug and p.returncode != 0:
            print('Trying link command:')
            print(' '.join(cmd))
            print('Output was:')
            print('   ',b'   '.join(lines).decode())
        returncode = p.returncode
    except (IOError,OSError) as e:
        if debug:
            print('Trying link command:')
            print(' '.join(cmd))
            print('Caught error: ',repr(e))
        returncode = 1

    if returncode:
        # The linker needs to be a c++ linker, which isn't 'cc'.  However, I couldn't figure
        # out how to get setup.py to tell me the actual command to use for linking.  All the
        # executables available from build_ext.compiler.executables are 'cc', not 'c++'.
        # I think this must be related to the bugs about not handling c++ correctly.
        #    http://bugs.python.org/issue9031
        #    http://bugs.python.org/issue1222585
        # So just switch it manually and see if that works.
        if 'clang' in cc:
            cpp = cc.replace('clang', 'clang++')
        elif 'icc' in cc:
            cpp = cc.replace('icc', 'icpc')
        elif 'gcc' in cc:
            cpp = cc.replace('gcc', 'g++')
        elif ' cc' in cc:
            cpp = cc.replace(' cc', ' c++')
        elif cc == 'cc':
            cpp = 'c++'
        else:
            comp_type = get_compiler_type(compiler)
            if comp_type == 'gcc':
                cpp = 'g++'
            elif comp_type == 'clang':
                cpp = 'clang++'
            elif comp_type == 'icc':
                cpp = 'g++'
            else:
                cpp = 'c++'
        cmd = [cpp] + compiler.linker_so[1:] + lflags + [o_name,'-o',exe_name]
        if debug:
            print('cmd = ',' '.join(cmd))
        try:
            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            lines = p.stdout.readlines()
            p.communicate()
            if debug and p.returncode != 0:
                print('Trying link command:')
                print(' '.join(cmd))
                print('Output was:')
                print('   ',b'   '.join(lines).decode())
            returncode = p.returncode
        except (IOError,OSError) as e:
            if debug:
                print('Trying to link using command:')
                print(' '.join(cmd))
                print('Caught error: ',repr(e))
            returncode = 1

    # Remove the temp files
    if returncode != 0:
        # Don't delete files in case helpful for troubleshooting.
        return False
    else:
        os.remove(cpp_name)
        os.remove(o_name)
        if os.path.exists(exe_name):
            os.remove(exe_name)
        return True

def try_cpp(compiler, cflags=[], lflags=[], prepend=None):
    """Check if compiling a simple bit of c++ code with the given compiler works properly.
    """
    from textwrap import dedent
    cpp_code = dedent("""
    #include <iostream>
    #include <vector>
    #include <cmath>

    int main() {
        int n = 500;
        std::vector<double> x(n,0.);
        for (int i=0; i<n; ++i) x[i] = 2*i+1;
        double sum=0.;
        for (int i=0; i<n; ++i) sum += std::log(x[i]);
        return sum;
    }
    """)
    return try_compile(cpp_code, compiler, cflags, lflags, prepend=prepend)

def check_ffi_compile(compiler, cflags=[], lflags=[]):
    ffi_code = """
#include "ffi/ffi.h"
int main() {
    return 0;
}
"""
    print("Checking if you have ffi installed on your system...")
    if try_compile(ffi_code, compiler, cflags, lflags):
        print('Found "ffi/ffi.h"')
        return True
    else:
        print('Unable to compile file with #include "ffi/ffi.h"')
        print('Trying ffi.h instead...')
        ffi_code = ffi_code.replace('ffi/ffi.h', 'ffi.h')
        if try_compile(ffi_code, compiler, cflags, lflags):
            print('Found "ffi.h"')
            return True
        else:
            print("Unable to compile when including either ffi/ffi.h or just ffi.h")
            return False

# Based on recipe 577058: http://code.activestate.com/recipes/577058/
def query_yes_no(question, default="yes", timeout=30):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is one of "yes" or "no".
    """
    valid = {"yes":"yes",   "y":"yes",  "ye":"yes",
             "no":"no",     "n":"no"}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while 1:
        sys.stdout.write(question + prompt)
        sys.stdout.flush()
        i, _, _ = select.select( [sys.stdin], [], [], timeout )

        if i:
            choice = sys.stdin.readline().strip()
        else:
            sys.stdout.write("\nPrompt timed out after %s seconds.\n"%timeout)
            return default

        if default is not None and choice == '':
            return default
        elif choice in valid.keys():
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")


def check_ffi(compiler, cflags=[], lflags=[]):
    try:
        import cffi
    except ImportError:
        # Then cffi will need to be installed.
        # It requires libffi, so check if it is available.
        if check_ffi_compile(compiler, cflags, lflags):
            return
        # libffi needs to be installed.  Give a helpful message about how to do so.
        prefix = '/SOME/APPROPRIATE/PREFIX'
        prefix_param = [param for param in sys.argv if param.startswith('--prefix=')]
        if len(prefix_param) == 1:
            prefix = prefix_param[0].split('=')[1]
            prefix = os.path.expanduser(prefix)
        msg = """
WARNING: Coord uses cffi, which in turn requires libffi to be installed.
         As the latter is not a python package, pip cannot download and
         install it.  However, it is fairly straightforward to install.

On Linux, you can use one of the following:

    apt-get install libffi-dev
    yum install libffi-devel

On a Mac, it should be available after you do:

    xcode-select --install

If neither of those work for you, you can install it yourself with the
following commands:

    wget ftp://sourceware.org:/pub/libffi/libffi-3.2.1.tar.gz
    tar xfz libffi-3.2.1.tar.gz
    cd libffi-3.2.1
    ./configure --prefix={0}
    make
    make install
    cp */include/ffi*.h {0}/include
    cd ..

If you have already done this, then check the command (given above) that failed.  You may
need to add a directory to either C_INCLUDE_PATH, LIBRARY_PATH, or LD_LIBRARY_PATH to
make it succeed.
""".format(prefix)
        print(msg)
        q = "Stop the installation here to take care of this?"
        yn = query_yes_no(q, default='yes')
        if yn == 'yes':
            sys.exit(1)

def fix_compiler(compiler):
    # Remove any -Wstrict-prototypes in the compiler flags (since invalid for C++)
    try:
        compiler.compiler_so.remove("-Wstrict-prototypes")
    except (AttributeError, ValueError):
        pass

    # Figure out what compiler it will use
    comp_type = get_compiler_type(compiler, output=True)
    cc = compiler.compiler_so[0]
    already_have_ccache = False
    if cc == 'ccache':
        already_have_ccache = True
        cc = compiler.compiler_so[1]
    if cc == comp_type:
        print('Using compiler %s'%(cc))
    else:
        print('Using compiler %s, which is %s'%(cc,comp_type))

    extra_cflags = copt[comp_type]
    extra_lflags = lopt[comp_type]

    success = try_cpp(compiler, extra_cflags, extra_lflags)
    if not success:
        # In case libc++ doesn't work, try letting the system use the default stdlib
        try:
            extra_cflags.remove('-stdlib=libc++')
            extra_lflags.remove('-stdlib=libc++')
        except (AttributeError, ValueError):
            pass
        else:
            success = try_cpp(compiler, extra_cflags, extra_lflags)
    if not success:
        print("There seems to be something wrong with the compiler or cflags")
        print(str(compiler.compiler_so))
        raise OSError("Compiler does not work for compiling C++ code")

    check_ffi(compiler, extra_cflags, extra_lflags)

    # Check if we can use ccache to speed up repeated compilation.
    if not already_have_ccache and try_cpp(compiler, prepend='ccache'):
        print('Using ccache')
        compiler.set_executable('compiler_so', ['ccache'] + compiler.compiler_so)

    # Return the extra cflags, since those will be added to the build step in a different place.
    print('Using extra flags ',extra_cflags)
    return extra_cflags, extra_lflags


# Make a subclass of build_ext so we can do different things depending on which compiler we have.
# In particular, we want to use different compiler options for OpenMP in each case.
# cf. http://stackoverflow.com/questions/724664/python-distutils-how-to-get-a-compiler-that-is-going-to-be-used
class my_builder( build_ext ):
    def build_extensions(self):
        cflags, lflags = fix_compiler(self.compiler)

        # Add the appropriate extra flags for that compiler.
        for e in self.extensions:
            e.extra_compile_args = cflags
            for flag in lflags:
                e.extra_link_args.append(flag)

        # Now run the normal build function.
        build_ext.build_extensions(self)

# AFAICT, setuptools doesn't provide any easy access to the final installation location of the
# executable scripts.  This bit is just to save the value of script_dir so I can use it later.
# cf. http://stackoverflow.com/questions/12975540/correct-way-to-find-scripts-directory-from-setup-py-in-python-distutils/
class my_easy_install( easy_install ):

    # Match the call signature of the easy_install version.
    def write_script(self, script_name, contents, mode="t", *ignored):
        # Run the normal version
        easy_install.write_script(self, script_name, contents, mode, *ignored)
        # Save the script install directory in the distribution object.
        # This is the same thing that is returned by the setup function.
        self.distribution.script_install_dir = self.script_dir


ext=Extension("coord._coord",
              sources,
              depends=headers,
              undef_macros=undef_macros,
              include_dirs=['include'])

with open('README.rst') as file:
    long_description = file.read()

# Read in the coord version from coord/_version.py
# cf. http://stackoverflow.com/questions/458550/standard-way-to-embed-version-into-python-package
version_file=os.path.join('coord','_version.py')
verstrline = open(version_file, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    coord_version = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (version_file,))
print('Coord version is %s'%(coord_version))

dist = setup(name="LSSTDESC.Coord",
      version=coord_version,
      author="LSST DESC (contact: Mike Jarvis)",
      author_email="michael@jarvis.net",
      description="Python module for handling angles and celestial coordinates",
      long_description=long_description,
      license = "MIT License",
      url="https://github.com/LSSTDESC/Coord",
      download_url="https://github.com/LSSTDESC/Coord/releases/tag/v%s.zip"%coord_version,
      packages=['coord'],
      package_data={'coord' : headers },
      ext_modules=[ext],
      install_requires=required,
      cmdclass = {'easy_install': my_easy_install},
)

# setup.py doesn't put the .so file in the coord directory, so this bit makes it possible to
# import coord from the root directory.  Not really advisable, but everyone does it at some
# point, so might as well facilitate it.
build_lib = glob.glob(os.path.join('build','*','coord','_coord*.so'))
if len(build_lib) >= 1:
    lib = os.path.join('coord','_coord.so')
    if os.path.lexists(lib): os.unlink(lib)
    os.link(build_lib[0], lib)

