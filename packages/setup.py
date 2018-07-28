from setuptools import setup, Extension, find_packages, distutils
from setuptools.command.build_ext import build_ext
import sys, os

# if we want to use conda compilers (which pick up llvm-openmp and fftw3 automatically)
if sys.platform == "darwin":
    os.environ["CC"]="clang"
    os.environ["CXX"]="clang++"

class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

ext_modules = [
    Extension(
        'basic_cpp',
        ['pyms/basic/basic_cpp.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            'pyms/mlpy/'
        ],
        libraries=[],
        extra_compile_args=[],
        extra_link_args=[],
        language='c++'
    ),
    Extension(
        'bandpass_filter_cpp',
        ['pyms/basic/bandpass_filter_cpp.cpp','pyms/basic/bandpass_filter_kernel.cpp'],
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
            get_pybind_include(user=True),
            'pyms/mlpy/'
        ],
        libraries=[],
        extra_compile_args=[],# ignored? hardcode below ( https://github.com/pybind/python_example/issues/23 )
        extra_link_args=[],
        language='c++'
    )]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except distutils.errors.CompileError:
            return False
    return True

# Note C++11 is needed by pybind11, not by the project:
def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.
    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-fopenmp')
            opts.append('-lfftw3')
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))   # AHB: C++11 not used now, but David Stein says needed.
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)
        
pkgs = find_packages()
print('found these packages:', pkgs)

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="ml_pyms",
    version="0.0.1",
    author="Jeremy Magland",
    author_email="",
    description="python processors for legacy mountainsort 3",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/magland/mountainsort",
    packages=pkgs,
    package_data={
        '': ['*.mp'], # Include all processor files
    },
    ext_modules=ext_modules,
    cmdclass={'build_ext': BuildExt},
    install_requires=
    [
        'pybind11>=2.2',
        'numpy',
        'numpydoc',
        'scipy',
        'matplotlib',
        'sklearn'
    ],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ),
)
