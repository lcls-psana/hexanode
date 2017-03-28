
# Example:
# /reg/g/psdm/sw/releases/ana-current/ConfigSvc/pyext/_ConfigSvc.pyx

# library header:
# /reg/common/package/hexanodelib/0.0.1/x86_64-centos7-gcc485/resort64c.h
# ~/lib/hexanode-lib/sort_non-LMF_from_1_detector/resort64c.h

# passing numpy arrays:
# http://stackoverflow.com/questions/17855032/passing-and-returning-numpy-arrays-to-c-methods-via-cython

# issue with outdated numpy
# https://docs.scipy.org/doc/numpy/reference/c-api.deprecations.html#deprecation-mechanism-npy-no-deprecated-api

# issue with resort64c.h:
# - seen from hexanode/src and app, but 
# - not seen in compilation stage from hexanode/includes or hexanode/pyext/hexanode_ext.pyx
# Hack to fix:
# export CPATH=/reg/common/package/hexanodelib/0.0.1/x86_64-centos7-gcc485


from cpython.bool cimport *
from cpython.float cimport *
from cpython.int cimport *
from cpython.list cimport *
from cpython.object cimport *
from cpython.ref cimport *
from cpython.string cimport *

#from libcpp.string cimport string
#from libcpp.vector cimport vector


import numpy as np
cimport numpy as np

#np.import_array()

cdef extern from "<stdint.h>" nogil:
    ctypedef   signed char  int8_t
    ctypedef   signed short int16_t
    ctypedef   signed int   int32_t
    ctypedef   signed long  int64_t
    ctypedef unsigned char  uint8_t
    ctypedef unsigned short uint16_t
    ctypedef unsigned int   uint32_t


ctypedef fused dtypes2d :
    np.ndarray[np.double_t, ndim=2]
    np.ndarray[np.int16_t,  ndim=2]
    np.ndarray[np.uint16_t, ndim=2]

#------------------------------
#------------------------------
#------------------------------
# TESTS
#------------------------------
#------------------------------
#------------------------------

def met1():
    """cython method for test call from python"""
    cdef int nelem = 1
    print "HELL of cython is here"


cdef extern from "hexanode/cfib.h":
    double cfib(int n)


def fib(n):
    """Returns the nth Fibonacci number"""
    return cfib(n)

#------------------------------

cdef extern from "hexanode/ctest_nda.h":
    #void ctest_nda[T](T *arr, int r, int c) except +
    void ctest_nda_f8(double   *arr, int r, int c) except +
    void ctest_nda_i2(int16_t  *arr, int r, int c) except +
    void ctest_nda_u2(uint16_t *arr, int r, int c) except +

#------------------------------
# Most compact working case 

#def test_nda_v1(np.ndarray nda):
#    print 'nda.dtype =', str(nda.dtype)
#    if   nda.dtype == np.float64 : ctest_nda(<double*>   nda.data, nda.shape[0], nda.shape[1])
#    elif nda.dtype == np.int16   : ctest_nda(<int16_t*>  nda.data, nda.shape[0], nda.shape[1])
#    elif nda.dtype == np.uint16  : ctest_nda(<uint16_t*> nda.data, nda.shape[0], nda.shape[1])
#    else: raise ValueError("Array data type is unknown")

#------------------------------
# Specialized methods for each type
#------------------------------

def test_nda_f8(np.ndarray[np.double_t, ndim=2] nda): ctest_nda_f8(<double*>  nda.data, nda.shape[0], nda.shape[1])
def test_nda_i2(np.ndarray[np.int16_t, ndim=2] nda):  ctest_nda_i2(<int16_t*> nda.data, nda.shape[0], nda.shape[1])
def test_nda_u2(np.ndarray[np.uint16_t, ndim=2] nda): ctest_nda_u2(<uint16_t*>nda.data, nda.shape[0], nda.shape[1])

#def test_nda_f8(np.ndarray[np.double_t, ndim=2] nda): ctest_nda(&nda[0,0], nda.shape[0], nda.shape[1])
#def test_nda_i2(np.ndarray[np.int16_t, ndim=2] nda):  ctest_nda(&nda[0,0], nda.shape[0], nda.shape[1])
#def test_nda_u2(np.ndarray[np.uint16_t, ndim=2] nda): ctest_nda(&nda[0,0], nda.shape[0], nda.shape[1])


#def test_nda(np.ndarray nda):
#def test_nda(dtypes2d nda):
#    print 'nda.dtype =', str(nda.dtype)
#    if   nda.dtype == np.float64 : test_nda_f8(nda)
#    elif nda.dtype == np.int16   : test_nda_i2(nda)
#    elif nda.dtype == np.uint16  : test_nda_u2(nda)
#    else: raise ValueError("Array data type is unknown")


#------------------------------
#------------------------------
#------------------------------

#def test_nda(np.ndarray nda):
#def test_nda(dtypes2d nda):
#    #print 'nda.dtype =', str(nda.dtype)
#    if   nda.dtype == np.float64 : ctest_nda_f8(<double*>   nda.data, nda.shape[0], nda.shape[1])
#    elif nda.dtype == np.int16   : ctest_nda_i2(<int16_t*>  nda.data, nda.shape[0], nda.shape[1])
#    elif nda.dtype == np.uint16  : ctest_nda_u2(<uint16_t*> nda.data, nda.shape[0], nda.shape[1])
#    else: raise ValueError("Array data type is unknown")

#------------------------------
#------------------------------
#------------------------------

## DOES NOT WORK
#def test_nda_xxx(np.ndarray nda):
#    ctest_nda(nda.data, nda.shape[0], nda.shape[1])

#------------------------------

# Very compact, BUT GENERATES TOO MANY WARNINGS AT scons

#def test_nda_xxx(dtypes2d nda):
#    print 'nda.dtype =', str(nda.dtype)
#    ctest_nda(&nda[0,0], nda.shape[0], nda.shape[1])

#------------------------------
#------------------------------
#------------------------------
#------------------------------

cdef extern from "hexanode/wrap_resort64c.h":
    void test_resort()

def ctest_resort():
    print 'In hexanode_ext.ctest_resort - call cpp method which uses hexanode_proxy/resort64c.h'
    test_resort()

#------------------------------

#cdef extern from "resort64c.h":
#    cdef cppclass sort_class:
#        sort_class()
#        int32_t sort() except +
#        int32_t run_without_sorting() except +


#cdef class c_sort_class:
#    """ Python wrapper for C++ class sort_class from resort64c.h. 
#    """
#    cdef sort_class* thisptr  # holds a C++ instance which we're wrapping

#    def __cinit__(self):
#        print "In c_sort_class.__cinit__"
#        self.thisptr = new sort_class();

#    def __dealloc__(self):
#        print "In c_sort_class.__dealloc__"
#        del self.thisptr

#    def c_sort(self):
#        print "In c_sort_class.c_sort"
#        return self.thisptr.sort()

#    def c_run_without_sorting(self):
#        print "In c_sort_class.c_run_without_sorting"
#        return self.thisptr.run_without_sorting()


#        int32_t run_without_sorting()
#        bool create_scalefactors_calibrator(bool _BuildStack, double _runtime_u, double _runtime_v, double _runtime_w, double _runtime_inner_fraction, double _fu,double _fv,double _fw)
#        int32_t init_after_setting_parameters()
#        bool do_calibration()
#        bool feed_calibration_data(bool build_map, double w_offset, int32_t number_of_correction_points)
#        bool feed_calibration_data(bool build_map, double w_offset)

#        void shift_layer_w(int32_t add_or_sub_sign,double w_offset)
#        void shift_sums(int32_t add_or_sub_sign,double sumu_offset,double sumv_offset)
#        void shift_sums(int32_t add_or_sub_sign,double sumu_offset,double sumv_offset,double sumw_offset)
#        void shift_position_origin(int32_t add_or_sub_sign,double pos_x_mm_offset,double pos_y_mm_offset)

#        bool are_all_single()
#        bool clone(sort_class *clone)
#        #static version_number_class get_version_number()


#------------------------------

cdef extern from "hexanode_proxy/resort64c.h":
    cdef cppclass signal_corrector_class:
        signal_corrector_class()


cdef class c_signal_corrector_class:
    """ Python wrapper for C++ class sort_class from resort64c.h. 
    """

    cdef signal_corrector_class* cobj  # holds a C++ instance which we're wrapping

    def __cinit__(self):
        print "XXX In c_signal_corrector_class.__cinit__",\
              " - direct test of methods from hexanode_proxy/resort64c.h in hexanode_ext.class c_signal_corrector_class"
        self.cobj = new signal_corrector_class();

    def __dealloc__(self):
        print "In c_signal_corrector_class.__dealloc__"
        del self.cobj

#------------------------------

