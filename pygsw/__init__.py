
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
    ctypes.c_int,
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS),
    numpy.ctypeslib.ndpointer(dtype=ctypes.c_double, flags=CONTIGUOUS)
]
pygsw_.calculate_pt.restype = None

def calculate_pt(long, lat, z, t, sp):
    t = numpy.array(t, copy=False, dtype=float, order='C')
    sp = numpy.array(sp, copy=False, dtype=float, order='C')
    assert t.shape == sp.shape
    long = numpy.array(numpy.broadcast_to(long % 360., t.shape), copy=False, dtype=float, order='C')
    lat = numpy.array(numpy.broadcast_to(lat, t.shape), copy=False, dtype=float, order='C')
    z = numpy.array(numpy.broadcast_to(z, t.shape), copy=False, dtype=float, order='C')
    pt = numpy.empty_like(t)
    pygsw_.calculate_pt(t.size, long, lat, z, t, sp, pt)
    return pt
