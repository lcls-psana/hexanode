#!/usr/bin/env python

print 'Use command: mpirun -n 2 python hexanode/examples/ex-11-MPIDS-save-h5.py'

"""
ipython
import psana
ds = psana.MPIDataSource('exp=xpptut15:run=390:smd')
events = ds.events()
evt = events.next()
"""

from time import time
import psana
from expmon.HexDataIO import HexDataIO
from pyimgalgos.GlobalUtils import print_ndarr

#------------------------------

DS_NAME = 'exp=xpptut15:run=390:smd' # sys.argv[1] ?
EVSKIP = 0
EVENTS = EVSKIP + 200000
OFNAME = './xpptut15-r390-e%06d-smd.h5'% EVENTS

#------------------------------

def proc_data():

    DIO = HexDataIO(dic_src_channels={'AmoETOF.0:Acqiris.0':(6,7,8,9,10,11),'AmoITOF.0:Acqiris.0':(0,)})
    DIO.open_input_dataset(DS_NAME, pbits=0, do_mpids=True)
    #ds = psana.MPIDataSource('exp=xpptut15:run=390:smd')

    print 'DIO starttime: %s' % DIO.start_time()
    print 'DIO tdc_resolution [ns]: %.3f' % DIO.tdc_resolution()

    smldata = DIO.ds.small_data(OFNAME, keys_to_save=[], gather_interval=100)

    t0_sec = time()
    t1_sec = time()
 
    for nev,evt in enumerate(DIO.events()):

        if nev < EVSKIP: continue
        if nev > EVENTS: break

	if DIO.ds.rank == 0 : 
	    if nev<5\
	    or (nev<50 and (not nev%10))\
	    or (nev<500 and (not nev%100))\
	    or not nev%1000 :
                t1 = time()
                print 'Rank: %d event: %06d, dt(sec): %.3f' % (DIO.ds.rank, nev, t1 - t1_sec)
                t1_sec = t1

        if evt is None :
            print '  WARNING: evt is None, rank: %d' % DIO.ds.rank
            continue

        DIO.proc_waveforms_for_evt(evt)

        #print_ndarr(DIO._number_of_hits, '    number of hits', first=0, last=7)
        #print_ndarr(DIO._tdc_ns, '    TDC[ns]', first=0, last=5)
 
        # save per-event data
        smldata.event(nhits=DIO._number_of_hits, tdcns=DIO._tdc_ns)
 
    # get "summary" data
    #run_sum = smldata.sum(partial_run_sum)
    # save HDF5 file, including summary data
    print "consumed time (sec) = %.6f\n" % (time() - t0_sec)
    smldata.save(nevents=nev)

#------------------------------

if __name__ == "__main__" :
    import sys; global sys
    tname = sys.argv[1] if len(sys.argv) > 1 else '1'
    print 50*'_', '\nTest %s:' % tname

    proc_data()

    #if   tname == '1' : proc_data()
    #elif tname == '2' : test02()
    #elif tname == '3' : test03()
    #else : sys.exit('Test %s is not implemented' % tname)

    sys.exit('End of test %s' % tname)

#------------------------------
