# cython: language_level=3

cimport cython

cimport numpy
import numpy

cdef extern void cgsw_sa_from_sp(int n, double* sp, double* p, double* long, double* lat, double* sa)
cdef extern void cgsw_pt0_from_t(int n, double* sa, double* t, double* p, double* pt)
cdef extern void cgsw_ct_from_pt(int n, double* sa, double* pt, double* ct)
cdef extern void cgsw_nsquared_3d(int nx, int ny, int nz, double* h, double* sa, double* ct, double* p, double* lat, double* n2)

def calculate_pt(lon, lat, z, t, sp):
    cdef numpy.ndarray lon_, lat_, p, t_, sp_, pt0, ct
    cdef double[::1] SA
    lon_, lat_, p, t_, sp_ = [numpy.array(a, copy=False, dtype=float, order='C') for a in numpy.broadcast_arrays(lon % 360., lat, -z, t, sp)]
    SA = sa_from_sp(lon_.reshape(-1), lat_.reshape(-1), p.reshape(-1), sp_.reshape(-1))
    pt0 = numpy.empty_like(t_)
    pt0_from_t(SA, t_.reshape(-1), p.reshape(-1), pt0.reshape(-1))
    return pt0

def sa_from_sp(const double[::1] lon not None, const double[::1] lat not None, const double[::1] p not None, const double[::1] sp not None, double[::1] out=None):
    if out is None:
        out = numpy.empty_like(sp)
    cgsw_sa_from_sp(out.shape[0], &sp[0], &p[0], &lon[0], &lat[0], &out[0])
    return out

def pt0_from_t(const double[::1] SA not None, const double[::1] t not None, const double[::1] p not None, double[::1] out=None):
    if out is None:
        out = numpy.empty_like(t)
    cgsw_pt0_from_t(out.shape[0], &SA[0], &t[0], &p[0], &out[0])
    return out

def ct_from_pt(const double[::1] SA not None, const double[::1] pt0 not None, double[::1] out=None):
    if out is None:
        out = numpy.empty_like(pt0)
    cgsw_ct_from_pt(out.shape[0], &SA[0], &pt0[0], &out[0])
    return out

def nsquared(const double[:,:,::1] h not None, const double[:,:,::1] SA not None, const double[:,:,::1] ct not None, const double[:,:,::1] p not None, const double[:,::1] lat not None, double[:,:,::1] out=None):
    if out is None:
        out = numpy.empty((h.shape[0], h.shape[1], h.shape[2]), dtype=h.dtype)
    cgsw_nsquared_3d(h.shape[2], h.shape[1], h.shape[0], &h[0,0,0], &SA[0,0,0], &ct[0,0,0], &p[0,0,0], &lat[0,0], &out[0,0,0])