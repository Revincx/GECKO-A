from setuptools import setup, Extension
import os

# Path to the C source files
lib_dir = os.path.join('..', 'LIB')

# C source files for the extension
sources = [
    'geckoa_module.c',
    os.path.join(lib_dir, 'keyparameter.c'),
    os.path.join(lib_dir, 'keyflag.c'),
    os.path.join(lib_dir, 'minidict.c'),
    os.path.join(lib_dir, 'tempoci.c'),
    os.path.join(lib_dir, 'references.c'),
    os.path.join(lib_dir, 'tempflag.c'),
    os.path.join(lib_dir, 'sortstring.c'),
    os.path.join(lib_dir, 'toolbox.c'),
]

# Include directories
include_dirs = [lib_dir]

# Define the extension module
geckoa_ext = Extension(
    'geckoa',
    sources=sources,
    include_dirs=include_dirs,
    extra_compile_args=['-std=c99', '-O2'],
    extra_link_args=['-lm'],
)

setup(
    name='geckoa',
    version='0.1.0',
    description='GECKO-A atmospheric chemistry mechanism generator',
    author='GECKO-A Team',
    ext_modules=[geckoa_ext],
    py_modules=['geckoa'],
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
    ],
)
