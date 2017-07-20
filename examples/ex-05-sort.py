#!/usr/bin/env python
#------------------------------

import sys
import hexanode
import numpy as np
from time import time

#------------------------------

NUM_CHANNELS=32
NUM_IONS=16
FNAME_CALIBRATION_TABLE = "calibration_table.txt"
DO_PRINT = False # True

# #include "math.h"
# #include "hexanode_proxy/resort64c.h"
# #include <stdio.h>
# #include "stdlib.h"
# #include "memory.h"
# #include <iostream> // for cout, puts etc.
# #include <time.h>

# #include "hexanode/LMF_IO.h"
# #include "hexanode/SortUtils.h"
# 
# double	offset_sum_u, offset_sum_v, offset_sum_w;
# double	w_offset, pos_offset_x, pos_offset_y;
# int	command;
# sort_class * sorter;

#------------------------------

if __name__ == "__main__" :
    tname = sys.argv[1] if len(sys.argv) > 1 else '1'
    print 50*'_', '\nTest %s:' % tname
    print "syntax: sort_LMF filename\n"\
          "        This file will be sorted and\n"\
          "        a new file will be written.\n\n"
 
    if len(sys.argv) < 2 :
        print "Please provide a filename.\n"
        sys.exit(0)

    if len(sys.argv) > 2 :
        print "too many arguments\n"
        sys.exit(0)
             
    tdc_ns = np.zeros((NUM_CHANNELS, NUM_IONS), dtype=np.float64)
    number_of_hits = np.zeros((NUM_CHANNELS,), dtype=np.int32)

    command = -1;
 
#   // The "command"-value is set in the first line of "sorter.txt"
#   // 0 = only convert to new file format
#   // 1 = sort and write new file 
#   // 2 = calibrate fv, fw, w_offset
#   // 3 = create calibration table files

#   // create the sorter:
    sorter = None
    sorter = hexanode.py_sort_class()
    fname_cfg = "sorter.txt"
    status, command, offset_sum_u, offset_sum_v, offset_sum_w, w_offset, pos_offset_x, pos_offset_y=\
        hexanode.py_read_config_file(fname_cfg, sorter)
    print 'read_config_file status, command, offset_sum_u, offset_sum_v, offset_sum_w, w_offset, pos_offset_x, pos_offset_y=',\
                            status, command, offset_sum_u, offset_sum_v, offset_sum_w, w_offset, pos_offset_x, pos_offset_y

    if not status :
        print "WARNING: can't read config file %s" % fname_cfg
        del sorter
        sys.exit(0)

    print 'use_sum_correction', sorter.use_sum_correction
    print 'use_pos_correction', sorter.use_pos_correction
    if sorter is not None :
        if sorter.use_sum_correction or sorter.use_pos_correction :
            status = hexanode.py_read_calibration_tables(FNAME_CALIBRATION_TABLE, sorter)

    if command == -1 :
   	print "no config file was read. Nothing to do."
        if sorter is not None : del sorter
        sys.exit(0)

    print "Numeration of channels:"
    print 'Cu1', sorter.cu1 
    print 'Cu2', sorter.cu2 
    print 'Cv1', sorter.cv1 
    print 'Cv2', sorter.cv2 
    print 'Cw1', sorter.cw1 
    print 'Cw2', sorter.cw2 
    print 'Cmcp',sorter.cmcp 

#   char error_text[512];

    LMF_Filename = sys.argv[1]
    LMF = hexanode.lmf_io(NUM_CHANNELS, NUM_IONS)
    if not LMF.open_input_lmf(LMF_Filename) :
        print "Can't open file: %s" % LMF_Filename
        sys.exit(0)

    print 'LMF starttime: : %s' % LMF.start_time()
    print 'LMF stoptime   : %s' % LMF.stop_time()

#   // initialization of the sorter:
    print "init sorter... "

    sorter.set_tdc_resolution_ns(0.025)
    sorter.set_tdc_array_row_length(NUM_IONS)
    sorter.set_count(number_of_hits)
    sorter.set_tdc_pointer(tdc_ns)

    if command >= 2 :
        sorter.create_scalefactors_calibrator(True,\
                                              sorter.runtime_u,\
                                              sorter.runtime_v,\
                                              sorter.runtime_w, 0.78,\
                                              sorter.fu, sorter.fv, sorter.fw)

    error_code = sorter.init_after_setting_parameters()
    if error_code :
   	print "sorter could not be initialized\n"
        error_text = sorter.get_error_text(error_code, 512)
        print 'Error %d: %s' % (error_code, error_text)
        sys.exit(0)

    print "ok for sorter initialization\n"

    print "LMF.tdcresolution %f\n" % LMF.tdc_resolution

#   while (my_kbhit()); // empty keyboard buffer

    event_counter = 0

    t_sec = time()

    print "reading event data... \n"

    while LMF.read_next_event() :

        #number_of_channels = LMF.get_number_of_channels()
        event_number = LMF.get_event_number()
    	event_counter+=1
 
	if not event_number%10000 :
            print 'Event number: %06d' % event_number

#   	//if (event_counter%10000 == 0) {if (my_kbhit()) break;}

#       //==================================
#       // TODO by end user:
#   	// Here you must read in a data block from your data file
#   	// and fill the array tdc_ns[][] and number_of_hits[]

        #nhits = np.zeros((NUMBER_OF_CHANNELS,), dtype=np.int32)
        LMF.get_number_of_hits_array(number_of_hits)
        if LMF.error_flag :
            error_text = LMF.get_error_text(LMF.error_flag)
            print "LMF Error %d: %s" % (LMF.error_flag, error_text)
            sys.exit(0)
        if DO_PRINT : print '   number_of_hits_array', number_of_hits[:8]


        #dtdc = np.zeros((NUMBER_OF_CHANNELS, NUMBER_OF_HITS), dtype=np.float64)
        LMF.get_tdc_data_array(tdc_ns)
        if LMF.error_flag :
            error_text = LMF.get_error_text(LMF.error_flag)
            print "LMF Error %d: %s" % (LMF.error_flag, error_text)
            sys.exit(0)
        if DO_PRINT : print '   TDC data:\n', tdc_ns[0:8,0:5]


#   	// apply conversion to ns
        if True :
            tdc_ns *= LMF.tdc_resolution

#       //==================================
#   	// TODO by end user...

        if sorter.use_hex :        
  	    # shift the time sums to zero:
   	    sorter.shift_sums(+1, offset_sum_u, offset_sum_v, offset_sum_w)
   	    #shift layer w so that the middle lines of all layers intersect in one point:
   	    sorter.shift_layer_w(+1, w_offset)
        else :
            # shift the time sums to zero:
            sorter.shift_sums(+1, offset_sum_u, offset_sum_v)

   	# shift all signals from the anode so that the center of the detector is at x=y=0:
   	sorter.shift_position_origin(+1, pos_offset_x, pos_offset_y)
 
   	sorter.feed_calibration_data(True, w_offset) # for calibration of fv, fw, w_offset and correction tables

#   	if (sorter->scalefactors_calibrator) {
#   		if (sorter->scalefactors_calibrator->map_is_full_enough()) break;
#   	}
#   	
# 
#   	int number_of_particles = 0;
#   	if (command == 1) {  // sort the TDC-Data and reconstruct missing signals
#   		// sorts/reconstructs the detector signals and apply the sum- and NL-correction.
#   		number_of_particles = sorter->sort();
#   		// "number_of_particles" is the number of reconstructed particles
#   	} else {
#   		number_of_particles = sorter->run_without_sorting();
#   	}
# 
#   	double u = tdc_ns[Cu1][0] + tdc_ns[Cu2][0] - 2*tdc_ns[Cmcp][0];
#   	double v = tdc_ns[Cv1][0] + tdc_ns[Cv2][0] - 2*tdc_ns[Cmcp][0];
#   	double w = tdc_ns[Cw1][0] + tdc_ns[Cw2][0] - 2*tdc_ns[Cmcp][0];
# 
#   	if (false) {
#   	  printf("  Event %5i  number_of_particles: %i", event_number, number_of_particles);
#   	  for(int i=0; i<number_of_particles; i++) {
#   	    printf("\n    p:%1i x:%.3f y:%.3f t:%.3f", i, sorter->output_hit_array[i]->x, 
#                                                               sorter->output_hit_array[i]->y,
#   	  	                                        sorter->output_hit_array[i]->time);
#   	  }
#   	  printf("\n    u:%.3f v:%.3f w:%.3f\n", u, v, w);
#   	} 
# 
# 
#   	//printf("error	Seaqrch for TODO by end user...");
#   	//#error	TODO by end user:
#   	// write the results into a new data file.
#   	// the variable "number_of_particles" contains the number of
#   	// reconstructed particles.
#   	// the x and y  (in mm) and TOF (in ns) is stored in the array sorter->output_hit_array:
#   	// 
#   	// for the first particle:
#   	// sorter->output_hit_array[0]->x;
#   	// sorter->output_hit_array[0]->y;
#   	// sorter->output_hit_array[0]->time;
#   	// 
#   	// for the 2nd particle:
#   	// sorter->output_hit_array[1]->x;
#   	// sorter->output_hit_array[1]->y;
#   	// sorter->output_hit_array[1]->time;
#   	//
#   	// for each particle you can also retrieve the information about how the particle
#   	// was reconstructed (tog et some measure of the confidence):
#   	// sorter->output_hit_array[0]->method;
# 
# 
#   } // end of the while loop
# 
# 
    if command == 2 :
        print "TBD !!!!!!!!!! calibrating detector... "
        sorter.do_calibration()
        print "ok  after do_calibration \n"

#   	if (sorter->scalefactors_calibrator) {
#   	       printf("Good calibration factors are:\nf_U =%lg\nf_V =%lg\nf_W =%lg\nOffset on layer W=%lg\n",
#                           2.*sorter->fu, 
#                           2.*sorter->scalefactors_calibrator->best_fv,
#                           2.*sorter->scalefactors_calibrator->best_fw, 
#                           sorter->scalefactors_calibrator->best_w_offset);
#   	}


    if (command == 3) : # generate and print correction tables for sum- and position-correction
        print "TBD !!!!!!!!!!! creating calibration tables...\n"
#       create_calibration_tables(FNAME_CALIBRATION_TABLE, sorter);
        print "\nfinished creating calibration tables: %s\n" % FNAME_CALIBRATION_TABLE


    print "consumed time (sec) = %.6f\n" % (time() - t_sec)

    if sorter is not None : del sorter
    sys.exit(0)
