#--------------------------------------------------------------------------
# File and Version Information:
#  $Id: SConscript 13182 2017-02-22 20:25:58Z davidsch@SLAC.STANFORD.EDU $
#
# Description:
#  SConscript file for package hexanode
#------------------------------------------------------------------------

# Do not delete following line, it must be present in 
# SConscript file for any SIT project
Import('*')

#
# For the standard SIT packages which build libraries, applications,
# and Python modules it is usually sufficient to call
# standardSConscript() function which defines rules for all
# above targets. Many standard packages do not need any special options, 
# but those which need can modify standardSConscript() behavior using
# a number of arguments, here is a complete list:
#
#    LIBS - list of additional libraries needed by this package
#    LIBPATH - list of directories for additional libraries
#    BINS - dictionary of executables and their corresponding source files
#    TESTS - dictionary of test applications and their corresponding source files
#    SCRIPTS - list of scripts in app/ directory
#    UTESTS - names of the unit tests to run, if not given then all tests are unit tests
#    PYEXTMOD - name of the Python extension module, package name used by default
#    CCFLAGS - additional flags passed to C/C++ compilers
#    NEED_QT - set to True to enable Qt support
#
#

# for our RPM installation, we will assume hexanode_proxy is installed
# and provides the lib Resort64c_x64. No other package should try to use it
# directly.
standardSConscript(PYEXTMOD="hexanode_ext", LIBS="Resort64c_x64")
#standardSConscript(PYEXTMOD="hexanode_ext",\
#  CCFLAGS="-I/reg/common/package/hexanodelib/0.0.1/x86_64-centos7-gcc485/",\
#  LIBPATH="/reg/common/package/hexanodelib/0.0.1/x86_64-centos7-gcc485",\
#  LIBS="Resort64c_x64"
#)
