# cython: language_level=3

cimport cython

cimport numpy
import numpy

cdef extern void cgsw_sa_from_sp(int n, const double* sp, const double* p, const double* long, const double* lat, double* sa)
cdef extern void cgsw_pt0_from_t(int n, const double* sa, const double* t, const double* p, double* pt)
cdef extern void cgsw_ct_from_pt(int n, const double* sa, const double* pt, double* ct)
cdef extern void cgsw_pt_from_ct(int n, const double* sa, const double* ct, double* pt)
cdef extern void cgsw_rho(int n, const double* sa, const double* ct, const double* p, double* rho)
cdef extern void cgsw_nsquared_3d(int nx, int ny, int nz, const int* mask, const double* h, const double* sa, const double* ct, const double* p, const double* lat, double* n2)

def calculate_pt(lon, lat, z, t, sp):
    cdef numpy.ndarray lon_, lat_, p, t_, sp_, pt0, ct
    cdef double[::1] SA
    lon_, lat_, p, t_, sp_ = [numpy.array(a, copy=False, dtype=float, order='C') for a in numpy.broadcast_arrays(lon % 360., lat, -z, t, sp)]
    SA = sa_from_sp(lon_.reshape(-1), lat_.reshape(-1), p.reshape(-1), sp_.reshape(-1))
    pt0 = numpy.empty_like(t_)
    pt0_from_t(SA, t_.reshape(-1), p.reshape(-1), pt0.reshape(-1))
    return pt0

def sa_from_sp(const double[::1] lon not None, const double[::1] lat not None, const double[::1] p not None, const double[::1] sp not None, double[::1] out=None):
    assert lon.shape[0] == lat.shape[0], 'Sizes of longitude (%i) and latitude (%i) do not match.' % (lon.shape[0], lat.shape[0])
    assert lon.shape[0] == p.shape[0], 'Sizes of longitude (%i) and pressure (%i) do not match.' % (lon.shape[0], p.shape[0])
    assert lon.shape[0] == sp.shape[0], 'Sizes of longitude (%i) and practical salinity (%i) do not match.' % (lon.shape[0], sp.shape[0])
    if out is None:
        out = numpy.empty_like(sp)
    assert sp.shape[0] == out.shape[0], 'Sizes of practical salinity (%i) and output (%i) do not match.' % (sp.shape[0], out.shape[0])
    cgsw_sa_from_sp(<int>out.shape[0], &sp[0], &p[0], &lon[0], &lat[0], &out[0])
    return out

def pt0_from_t(const double[::1] SA not None, const double[::1] t not None, const double[::1] p not None, double[::1] out=None):
    assert SA.shape[0] == t.shape[0], 'Sizes of absolute salinity (%i) and in-situ temperature (%i) do not match.' % (SA.shape[0], t.shape[0])
    assert SA.shape[0] == p.shape[0], 'Sizes of absolute salinity (%i) and pressure (%i) do not match.' % (SA.shape[0], p.shape[0])
    if out is None:
        out = numpy.empty_like(t)
    assert SA.shape[0] == out.shape[0], 'Sizes of absolute salinity (%i) and output (%i) do not match.' % (SA.shape[0], out.shape[0])
    cgsw_pt0_from_t(<int>out.shape[0], &SA[0], &t[0], &p[0], &out[0])
    return out

def ct_from_pt(const double[::1] SA not None, const double[::1] pt0 not None, double[::1] out=None):
    assert SA.shape[0] == pt0.shape[0], 'Sizes of absolute salinity (%i) and potential temperature (%i) do not match.' % (SA.shape[0], pt0.shape[0])
    if out is None:
        out = numpy.empty_like(pt0)
    assert SA.shape[0] == out.shape[0], 'Sizes of absolute salinity (%i) and output (%i) do not match.' % (SA.shape[0], out.shape[0])
    cgsw_ct_from_pt(<int>out.shape[0], &SA[0], &pt0[0], &out[0])
    return out

def pt_from_ct(const double[::1] SA not None, const double[::1] ct not None, double[::1] out=None):
    assert SA.shape[0] == ct.shape[0], 'Sizes of absolute salinity (%i) and conservative temperature (%i) do not match.' % (SA.shape[0], ct.shape[0])
    if out is None:
        out = numpy.empty_like(ct)
    assert SA.shape[0] == out.shape[0], 'Sizes of absolute salinity (%i) and output (%i) do not match.' % (SA.shape[0], out.shape[0])
    cgsw_pt_from_ct(<int>out.shape[0], &SA[0], &ct[0], &out[0])
    return out

def rho(const double[::1] SA not None, const double[::1] ct not None, const double[::1] p not None, double[::1] out=None):
    assert SA.shape[0] == ct.shape[0], 'Sizes of absolute salinity (%i) and conservative temperature (%i) do not match.' % (SA.shape[0], ct.shape[0])
    assert SA.shape[0] == p.shape[0], 'Sizes of absolute salinity (%i) and pressure (%i) do not match.' % (SA.shape[0], p.shape[0])
    if out is None:
        out = numpy.empty_like(SA)
    assert SA.shape[0] == out.shape[0], 'Sizes of absolute salinity (%i) and output (%i) do not match.' % (SA.shape[0], out.shape[0])
    cgsw_rho(<int>out.shape[0], &SA[0], &ct[0], &p[0], &out[0])
    return out

def nsquared(const int[:,::1] mask not None, const double[:,:,::1] h not None, const double[:,:,::1] SA not None, const double[:,:,::1] ct not None, const double[:,:,::1] p not None, const double[:,::1] lat not None, double[:,:,::1] out=None):
    assert h.shape[1] == mask.shape[0] and h.shape[2] == mask.shape[1]
    assert h.shape[0] == p.shape[0] and h.shape[1] == p.shape[1] and h.shape[2] == p.shape[2]
    assert h.shape[0] == SA.shape[0] and h.shape[1] == SA.shape[1] and h.shape[2] == SA.shape[2]
    assert h.shape[0] == ct.shape[0] and h.shape[1] == ct.shape[1] and h.shape[2] == ct.shape[2]
    assert h.shape[1] == lat.shape[0] and h.shape[2] == lat.shape[1]
    if out is None:
        out = numpy.empty((h.shape[0] - 1, h.shape[1], h.shape[2]), dtype=h.dtype)
    assert out.shape[0] == p.shape[0] - 1 and out.shape[1] == p.shape[1] and out.shape[2] == p.shape[2]
    if h.shape[0] == 1:
        return
    cgsw_nsquared_3d(<int>h.shape[2], <int>h.shape[1], <int>h.shape[0], &mask[0,0], &h[0,0,0], &SA[0,0,0], &ct[0,0,0], &p[0,0,0], &lat[0,0], &out[0,0,0])