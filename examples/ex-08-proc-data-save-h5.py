#!/usr/bin/env python
#------------------------------
"""Reads hexanode-acqiris events in loop
   applys CFD to 7 waveforms and find number of hits/channel and array of hit times
   saves number of hits/channel and array of hit times in hdf5 file
"""
#------------------------------
import os
import sys
#import numpy as np
from time import time
from expmon.HexDataIO import HexDataIO
#from pyimgalgos.GlobalUtils import print_ndarr

#------------------------------

DS_NAME = 'exp=xpptut15:run=390' # sys.argv[1]
EVSKIP = 0
EVENTS = EVSKIP + 200000
OFNAME = './xpptut15-r390-e200k-nhits-tdcns.h5' #'./test.h5'

#------------------------------

def proc_data():

    DIO = HexDataIO(dic_src_channels={'AmoETOF.0:Acqiris.0':(6,7,8,9,10,11),'AmoITOF.0:Acqiris.0':(0,)})
    DIO.open_input_dataset(DS_NAME, pbits=0) # pbits=1022

    print 'DIO starttime: : %s' % DIO.start_time()
    print 'DIO stoptime   : %s' % DIO.stop_time()
    tdc_res_ns = DIO.tdc_resolution()
    print 'DIO tdc_resolution : %.3f' % tdc_res_ns

    DIO.open_output_h5file(OFNAME)

    event_number = 0
    t0_sec = time()
    t1_sec = time()
    while DIO.read_next_event() :

        #number_of_channels = DIO.get_number_of_channels()
        event_number = DIO.get_event_number()

        if event_number < EVSKIP : continue
        if event_number >= EVENTS : break

        #print 'Event number: %06d' % event_number
	if event_number<5\
	or (event_number<50 and (not event_number%10))\
	or (event_number<500 and (not event_number%100))\
	or not event_number%1000 :
            t1 = time()
            print 'Event: %06d, dt(sec): %.3f' % (event_number, t1 - t1_sec)
            t1_sec = t1

        DIO.add_event_to_h5file()

    print "consumed time (sec) = %.6f\n" % (time() - t0_sec)
    DIO.close_output_h5file(pbits=1)

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
