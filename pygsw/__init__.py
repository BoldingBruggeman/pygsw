
from __future__ import print_function
import sys
import os
import ctypes

try:
   import numpy
except ImportError:
   print('Unable to import NumPy. Please ensure it is installed.')
   sys.exit(1)

# Determine potential names of dynamic library.
if os.name=='nt':
   dllpaths = ('pygsw_.dll','libpygsw_.dll')
elif os.name == "posix" and sys.platform == 'darwin':
   dllpaths = ('libpygsw_.dylib',)
else:
   dllpaths = ('libpygsw_.so',)

def find_library(basedir):
    for dllpath in dllpaths:
        dllpath = os.path.join(basedir, dllpath)
        if os.path.isfile(dllpath):
            return dllpath

dllpath = find_library(os.path.dirname(os.path.abspath(__file__)))
if not dllpath:
    for basedir in sys.path:
        dllpath = find_library(basedir)
        if dllpath:
            break

if not dllpath:
   print('Unable to locate pygsw dynamic library %s.' % (' or '.join(dllpaths),))
   sys.exit(1)

# Load library.
pygsw_ = ctypes.CDLL(str(dllpath))

CONTIGUOUS = str('CONTIGUOUS')

# Access to model objects (variables, parameters, dependencies, couplings, model instances)
pygsw_.calculate_pt.argtypes = [
    ctypes.c_double,
    ctypes.c_double,
    ctypes.c_int,
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS)
]
pygsw_.calculate_pt.restype = None

def calculate_pt(long, lat, z, t, sp):
    z = numpy.array(z, copy=False, dtype=float, order='C')
    t = numpy.array(t, copy=False, dtype=float, order='C')
    sp = numpy.array(sp, copy=False, dtype=float, order='C')
    assert z.shape == t.shape
    assert z.shape == sp.shape
    pt = numpy.empty_like(z)
    pygsw_.calculate_pt(long, lat, numpy.size(z), z, t, sp, pt)
    return pt
