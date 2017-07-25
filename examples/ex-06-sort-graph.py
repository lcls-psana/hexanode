#!/usr/bin/env python
#------------------------------

import sys
import hexanode
import numpy as np
from time import time
from math import sqrt

from pyimgalgos.GlobalUtils import print_ndarr
from pyimgalgos.HBins import HBins
#------------------------------

NUM_CHANNELS=32
NUM_IONS=16
FNAME_CALIBRATION_TABLE = "calibration_table.txt"
DO_PRINT = False # True

#PLOT_PREFIX = '2017-07-25-figs-hexanode/plot'
PLOT_PREFIX = 'plot'

PLOT_TIME_SUMS     = True  # False # True #
PLOT_TIME_SUMS_CORR= True  # False # True #
PLOT_RAW_UVW       = True  # False # True #
PLOT_UVW           = True  # False # True #
PLOT_CORRELATIONS  = True  # False # True #
PLOT_XY_COMPONENTS = True  # False # True #
PLOT_XY_2D         = True  # False # True #
PLOT_MISC          = True  # False # True #

#------------------------------

class Store :
    """Store of shared parameters.
    """
    def __init__(self) :

         self.lst_raw_u = []
         self.lst_raw_v = []
         self.lst_raw_w = []

         self.lst_u = []
         self.lst_v = []
         self.lst_w = []

         self.lst_time_sum_u = []
         self.lst_time_sum_v = []
         self.lst_time_sum_w = []

         self.lst_time_sum_u_corr = []
         self.lst_time_sum_v_corr = []
         self.lst_time_sum_w_corr = []

         self.lst_Xuv = []
         self.lst_Xuw = []
         self.lst_Xvw = []

         self.lst_Yuv = []
         self.lst_Yuw = []
         self.lst_Yvw = []

         self.lst_Deviation = []
         self.lst_consist_indicator = []
         self.lst_rec_method = []

         # images 
         nbins = 100
         self.img_x_bins = HBins((-1., 4.), nbins, vtype=np.float32)
         self.img_y_bins = HBins(( 1., 6.), nbins, vtype=np.float32)
         self.img_xy_uv = np.zeros((nbins, nbins), dtype=np.float32)
         self.img_xy_uw = np.zeros((nbins, nbins), dtype=np.float32)
         self.img_xy_vw = np.zeros((nbins, nbins), dtype=np.float32)
         self.img_xy_1  = np.zeros((nbins, nbins), dtype=np.float32)
         self.img_xy_2  = np.zeros((nbins, nbins), dtype=np.float32)

#------------------------------

sp = Store() 

#------------------------------
from pyimgalgos.GlobalGraphics import hist1d, show, move_fig, save_fig, move, save, plotImageLarge, plotGraph

def plot_image(img, figsize=(11,10), axwin=(0.10, 0.08, 0.88, 0.88), cmap='inferno',\
               title='x-y image', xlabel='x', ylabel='y', titwin=None, fnm='img.png' ) : #'gray_r'
    """
    """
    img_range = (sp.img_x_bins.vmin(), sp.img_x_bins.vmax(), sp.img_y_bins.vmax(), sp.img_y_bins.vmin()) 
    imgnb = img[1:-2,1:-2]
    amp_range = (0, imgnb.mean() + 4*imgnb.std())
    #amp_range = (0, 0.2*img.max())
    axim = plotImageLarge(img, img_range=img_range, amp_range=amp_range, figsize=figsize,\
                          title=title, origin='upper', window=axwin, cmap=cmap) # 'Greys') #'gray_r'
    axim.set_xlabel(xlabel, fontsize=18)
    axim.set_ylabel(ylabel, fontsize=18)
    axim.set_title(title,   fontsize=12)

    move(sp.hwin_x0y0[0], sp.hwin_x0y0[1])
    save('%s-%s' % (sp.prefix, fnm), sp.do_save)
    #show()

#------------------------------

def h1d(hlst, bins=None, amp_range=None, weights=None, color=None, show_stat=True, log=False,\
        figsize=(6,5), axwin=(0.15, 0.12, 0.78, 0.80), title='Title', xlabel='x', ylabel='y', titwin=None, fnm='hist.png') :
    """Wrapper for hist1d, move, and save methods, using common store parameters
    """
    fig, axhi, hi = hist1d(np.array(hlst), bins, amp_range, weights, color, show_stat,\
                           log, figsize, axwin, title, xlabel, ylabel, titwin)

    move(sp.hwin_x0y0[0], sp.hwin_x0y0[1])
    save('%s-%s' % (sp.prefix, fnm), sp.do_save)
    return fig, axhi, hi


def plot_graph(x, y, figsize=(7,6), pfmt='r-', lw=2,\
               title='py vs. px', xlabel='px', ylabel='py', fnm='graph.png') :
    fig, ax = plotGraph(x, y, figsize=figsize, pfmt=pfmt, lw=lw)
    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_title(title,   fontsize=12)

    move(sp.hwin_x0y0[0], sp.hwin_x0y0[1])
    save('%s-%s' % (sp.prefix, fnm), sp.do_save)

#------------------------------

def plot_histograms(prefix='plot', do_save=True, hwin_x0y0=(0,400)) :
    """Plots/saves histograms
    """
    sp.prefix    = prefix
    sp.do_save   = do_save
    sp.hwin_x0y0 = hwin_x0y0
    #---------

    #---------
    if PLOT_TIME_SUMS :
    #---------
        nbins = 160
        limits = (80,120)
    #---------
        #print_ndarr(sp.lst_time_sum_u, 'U')
        h1d(np.array(sp.lst_time_sum_u), bins=nbins, amp_range=limits,\
            title ='Time sum U', xlabel='Time sum U (ns)', ylabel='Events',\
            fnm='time_sum_u_ns.png')
    #---------
        #print_ndarr(sp.lst_time_sum_v, 'V')
        h1d(np.array(sp.lst_time_sum_v), bins=nbins, amp_range=limits,\
            title ='Time sum V', xlabel='Time sum V (ns)', ylabel='Events',\
            fnm='time_sum_v_ns.png')
    #---------
        #print_ndarr(sp.lst_time_sum_w, 'W')
        h1d(np.array(sp.lst_time_sum_w), bins=nbins, amp_range=limits,\
            title ='Time sum W', xlabel='Time sum W (ns)', ylabel='Events',\
            fnm='time_sum_w_ns.png')
    #---------

    #---------
    if PLOT_TIME_SUMS_CORR :
    #---------
        nbins = 160
        limits = (-20,20)
    #---------
        h1d(np.array(sp.lst_time_sum_u_corr), bins=nbins, amp_range=limits,\
            title ='Time sum U corrected', xlabel='Time sum U (ns) corrected', ylabel='Events',\
            fnm='time_sum_u_ns_corr.png')
    #---------
        h1d(np.array(sp.lst_time_sum_v_corr), bins=nbins, amp_range=limits,\
            title ='Time sum V corrected', xlabel='Time sum V (ns) corrected', ylabel='Events',\
            fnm='time_sum_v_ns_corr.png')
    #---------
        h1d(np.array(sp.lst_time_sum_w_corr), bins=nbins, amp_range=limits,\
            title ='Time sum W corrected', xlabel='Time sum W (ns) corrected', ylabel='Events',\
            fnm='time_sum_w_ns_corr.png')
    #---------

    #---------
    if PLOT_UVW :
    #---------
        nbins = 200
        limits = (-50,50)
    #---------
        h1d(np.array(sp.lst_u), bins=nbins, amp_range=limits,\
            title ='U', xlabel='U (mm)', ylabel='Events',\
            fnm='u_mm.png')
    #---------
        h1d(np.array(sp.lst_v), bins=nbins, amp_range=limits,\
            title ='V', xlabel='V (mm)', ylabel='Events',\
            fnm='v_mm.png')
    #---------
        h1d(np.array(sp.lst_w), bins=nbins, amp_range=limits,\
            title ='W', xlabel='W (mm)', ylabel='Events',\
            fnm='w_mm.png')
    #---------

    #---------
    if PLOT_RAW_UVW :
    #---------
        nbins = 200
        limits = (-100,100)
    #---------
        h1d(np.array(sp.lst_raw_u), bins=nbins, amp_range=limits,\
            title ='Raw U', xlabel='U (ns)', ylabel='Events',\
            fnm='raw_u_ns.png')
    #---------
        h1d(np.array(sp.lst_raw_v), bins=nbins, amp_range=limits,\
            title ='Raw V', xlabel='V (ns)', ylabel='Events',\
            fnm='raw_v_ns.png')
    #---------
        h1d(np.array(sp.lst_raw_w), bins=nbins, amp_range=limits,\
            title ='Raw W', xlabel='W (ns)', ylabel='Events',\
            fnm='raw_w_ns.png')
    #---------

    #---------
    if PLOT_CORRELATIONS :
    #---------
         #print_ndarr(sp.lst_time_sum_u, 'time_sum_u')
         #print_ndarr(sp.lst_raw_u,      'lst_raw_u ')

         plot_graph(sp.lst_raw_u, sp.lst_time_sum_u, figsize=(8,7), pfmt='b-', lw=2,\
            title='t sum vs. U', xlabel='U (ns)', ylabel='t sum U (ns)',\
            fnm='t_sum_vs_raw_u.png')

         plot_graph(sp.lst_raw_v, sp.lst_time_sum_v, figsize=(8,7), pfmt='b-', lw=2,\
            title='t sum vs. V', xlabel='V (ns)', ylabel='t sum V (ns)',\
            fnm='t_sum_vs_raw_v.png')

         plot_graph(sp.lst_raw_w, sp.lst_time_sum_w, figsize=(8,7), pfmt='b-', lw=2,\
            title='t sum vs. W', xlabel='W (ns)', ylabel='t sum W (ns)',\
            fnm='t_sum_vs_raw_w.png')

         #---------

         plot_graph(sp.lst_raw_u, sp.lst_time_sum_u_corr, figsize=(8,7), pfmt='b-', lw=2,\
            title='t sum corrected vs. U', xlabel='U (ns)', ylabel='t sum corrected U (ns)',\
            fnm='t_sum_corr_vs_raw_u.png')

         plot_graph(sp.lst_raw_v, sp.lst_time_sum_v_corr, figsize=(8,7), pfmt='b-', lw=2,\
            title='t sum_corrected vs. V', xlabel='V (ns)', ylabel='t sum corrected V (ns)',\
            fnm='t_sum_corr_vs_raw_v.png')

         plot_graph(sp.lst_raw_w, sp.lst_time_sum_w_corr, figsize=(8,7), pfmt='b-', lw=2,\
            title='t sum_corrected vs. W', xlabel='W(ns)', ylabel='t sum corrected W (ns)',\
            fnm='t_sum_corr_vs_raw_w.png')

    #---------
    if PLOT_XY_COMPONENTS :
    #---------
        nbins = 200
        limits = (-25,25)
    #---------
        h1d(np.array(sp.lst_Xuv), bins=nbins, amp_range=limits,\
            title ='Xuv', xlabel='Xuv (mm)', ylabel='Events',\
            fnm='Xuv_mm.png')
    #---------
        h1d(np.array(sp.lst_Xuw), bins=nbins, amp_range=limits,\
            title ='Xuw', xlabel='Xuw (mm)', ylabel='Events',\
            fnm='Xuw_mm.png')
    #---------
        h1d(np.array(sp.lst_Xvw), bins=nbins, amp_range=limits,\
            title ='Xvw', xlabel='Xvw (mm)', ylabel='Events',\
            fnm='Xvw_mm.png')
    #---------
        h1d(np.array(sp.lst_Yuv), bins=nbins, amp_range=limits,\
            title ='Yuv', xlabel='Yuv (mm)', ylabel='Events',\
            fnm='Yuv_mm.png')
    #---------
        h1d(np.array(sp.lst_Yuw), bins=nbins, amp_range=limits,\
            title ='Yuw', xlabel='Yuw (mm)', ylabel='Events',\
            fnm='Yuw_mm.png')
    #---------
        h1d(np.array(sp.lst_Yvw), bins=nbins, amp_range=limits,\
            title ='Yvw', xlabel='Yvw (mm)', ylabel='Events',\
            fnm='Yvw_mm.png')
    #---------

    #---------
    if PLOT_MISC :
    #---------
        h1d(np.array(sp.lst_Deviation), bins=100, amp_range=(0,5),\
            title ='Deviation', xlabel='Deviation (mm)', ylabel='Events',\
            fnm='deviation_mm.png')
    #---------
        h1d(np.array(sp.lst_consist_indicator), bins=64, amp_range=(0,64),\
            title ='Consistence indicator', xlabel='Consistence indicator (bit)', ylabel='Events',\
            fnm='consistence_indicator.png')
    #---------

    #---------
    if PLOT_XY_2D :
    #---------
        plot_image(sp.img_xy_uv, fnm='xy_uv.png', title='XY_uv image',   xlabel='x', ylabel='y', titwin='XY_uv image')
        plot_image(sp.img_xy_uw, fnm='xy_uw.png', title='XY_uw image',   xlabel='x', ylabel='y', titwin='XY_uw image')
        plot_image(sp.img_xy_vw, fnm='xy_vw.png', title='XY_vw image',   xlabel='x', ylabel='y', titwin='XY_vw image')
        plot_image(sp.img_xy_1,  fnm='xy_1.png',  title='XY image hit1', xlabel='x', ylabel='y', titwin='XY image hit1')
        plot_image(sp.img_xy_2,  fnm='xy_2.png',  title='XY image hit2', xlabel='x', ylabel='y', titwin='XY image hit2')
    #---------
        h1d(np.array(sp.lst_rec_method), bins=64, amp_range=(0,32),\
            title ='Reconstruction method', xlabel='Method id (bit)', ylabel='Events',\
            fnm='reconstruction_method.png')
    #---------

#------------------------------

def py_sort() :
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

    Cu1  = sorter.cu1 
    Cu2  = sorter.cu2 
    Cv1  = sorter.cv1 
    Cv2  = sorter.cv2 
    Cw1  = sorter.cw1 
    Cw2  = sorter.cw2 
    Cmcp = sorter.cmcp
    print "Numeration of channels - u1:%i  u2:%i  v1:%i  v2:%i  w1:%i  w2:%i  mcp:%i"%\
          (Cu1, Cu2, Cv1, Cv2, Cw1, Cw2, Cmcp)

    inds_of_channels    = (Cu1, Cu2, Cv1, Cv2, Cw1, Cw2)
    incr_of_consistence = (  1,   2,   4,   8,  16,  32)
    inds_incr = zip(inds_of_channels, incr_of_consistence)
    
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

    #sorter.set_use_reflection_filter_on_u1(True)
    #sorter.set_use_reflection_filter_on_u2(True)

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
    osqrt3 = 1./sqrt(3.)

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

        time_sum_u = tdc_ns[Cu1,0] + tdc_ns[Cu2,0] - 2*tdc_ns[Cmcp,0]
        time_sum_v = tdc_ns[Cv1,0] + tdc_ns[Cv2,0] - 2*tdc_ns[Cmcp,0]
        time_sum_w = tdc_ns[Cw1,0] + tdc_ns[Cw2,0] - 2*tdc_ns[Cmcp,0]

        raw_u = tdc_ns[Cu1,0] - tdc_ns[Cu2,0]
        raw_v = tdc_ns[Cv1,0] - tdc_ns[Cv2,0]
        raw_w = tdc_ns[Cw1,0] - tdc_ns[Cw2,0]

        u = raw_u * sorter.fu
        v = raw_v * sorter.fv
        w = (raw_w + w_offset) * sorter.fw

        Xuv = u
        Xuw = u
        Xvw = v + w
        Yuv = (u - 2*v)*osqrt3
        Yuw = (2*w - u)*osqrt3
        Yvw = (w - v)*osqrt3

        dX = Xuv - Xvw
        dY = Yuv - Yvw
        Deviation = sqrt(dX*dX + dY*dY)

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

        #print 'map_is_full_enough', hexanode.py_sorter_scalefactors_calibration_map_is_full_enough(sorter)
        sfco = hexanode.py_scalefactors_calibration_class(sorter)


        #LMF.get_tdc_data_array(tdc_ns)
        time_sum_u_corr = tdc_ns[Cu1,0] + tdc_ns[Cu2,0] - 2*tdc_ns[Cmcp,0]
        time_sum_v_corr = tdc_ns[Cv1,0] + tdc_ns[Cv2,0] - 2*tdc_ns[Cmcp,0]
        time_sum_w_corr = tdc_ns[Cw1,0] + tdc_ns[Cw2,0] - 2*tdc_ns[Cmcp,0]

        # break loop if statistics is enough
        if sfco :
            if sfco.map_is_full_enough() : 
                 print 'sfo.map_is_full_enough(): %s  event number: %06d' % (sfco.map_is_full_enough(), event_number)
                 break


        # Sort the TDC-Data and reconstruct missing signals and apply the sum- and NL-correction.
        # number_of_particles is the number of reconstructed particles
   	number_of_particles = sorter.sort() if command == 1 else\
                              sorter.run_without_sorting()

   	if False :
   	    print "  Event %5i  number_of_particles: %i" % (event_number, number_of_particles)
   	    for i in range(number_of_particles) :
                hco= hexanode.py_hit_class(sorter, i)
   	        print "    p:%1i x:%.3f y:%.3f t:%.3f met:%d" % (i, hco.x, hco.y, hco.time, hco.method)

   	    print "    part1 u:%.3f v:%.3f w:%.3f" % (u, v, w)

#       // TODO by end user..."

#        reflection_u1 = sorter.u1_reflection_time_position
#        reflection_u2 = sorter.u2_reflection_time_position
#        reflection_v1 = sorter.u1_reflection_time_position
#        reflection_v2 = sorter.u2_reflection_time_position
#        reflection_w1 = sorter.u1_reflection_time_position
#        reflection_w2 = sorter.u2_reflection_time_position

#        print 'reflections:', reflection_u1,  reflection_u2,\
#                              reflection_v1,  reflection_v2,\
#                              reflection_w1,  reflection_w2

        if number_of_particles<1 : continue

        if PLOT_RAW_UVW or PLOT_CORRELATIONS :
            sp.lst_raw_u.append(raw_u)
            sp.lst_raw_v.append(raw_v)
            sp.lst_raw_w.append(raw_w)

        if PLOT_TIME_SUMS or PLOT_CORRELATIONS :
            sp.lst_time_sum_u.append(time_sum_u)
            sp.lst_time_sum_v.append(time_sum_v)
            sp.lst_time_sum_w.append(time_sum_w)

            sp.lst_time_sum_u_corr.append(time_sum_u_corr)
            sp.lst_time_sum_v_corr.append(time_sum_v_corr)
            sp.lst_time_sum_w_corr.append(time_sum_w_corr)

        if PLOT_UVW :
            sp.lst_u.append(u)
            sp.lst_v.append(v)
            sp.lst_w.append(w)

        if PLOT_XY_COMPONENTS :
            sp.lst_Xuv.append(Xuv)
            sp.lst_Xuw.append(Xuw)
            sp.lst_Xvw.append(Xvw)
                             
            sp.lst_Yuv.append(Yuv)
            sp.lst_Yuw.append(Yuw)
            sp.lst_Yvw.append(Yvw)

        if PLOT_MISC :
            sp.lst_Deviation.append(Deviation)
            
            # fill Consistence Indicator
            consistenceIndicator = 0
            for (ind, incr) in inds_incr :
              if number_of_hits[ind]>0 : consistenceIndicator += incr
            sp.lst_consist_indicator.append(consistenceIndicator)

        if PLOT_XY_2D :
            # fill 2-d images
            hco= hexanode.py_hit_class(sorter, 0)
            x1, y1, method = hco.x, hco.y, hco.method

            x2, y2 = (-10,-10) 
            if number_of_particles > 1 :
                hco2 = hexanode.py_hit_class(sorter, 1)
                x2, y2 = hco2.x, hco2.y

            sp.lst_rec_method.append(method)
            #print 'reconstruction method %d' % method

            ix1, ix2, ixuv, ixuw, ixvw = sp.img_x_bins.bin_indexes((x1, x2, Xuv, Xuw, Xvw))
            iy1, iy2, iyuv, iyuw, iyvw = sp.img_y_bins.bin_indexes((y1, y2, Yuv, Yuw, Yvw))

            sp.img_xy_1 [iy1,  ix1]  += 1
            sp.img_xy_2 [iy2,  ix2]  += 1
            sp.img_xy_uv[iyuv, ixuv] += 1
            sp.img_xy_uw[iyuw, ixuw] += 1 
            sp.img_xy_vw[iyvw, ixvw] += 1 

#   	// write the results into a new data file.
#   	// the variable "number_of_particles" contains the number of reconstructed particles.
#   	// the x and y (in mm) and TOF (in ns) is stored in the array sorter->output_hit_array:

#   	// for the i-th particle (i starts from 0):
#       // hco= hexanode.py_hit_class(sorter, i)
#       // hco.x, hco.y, hco.time

#   	// for each particle you can also retrieve the information about how the particle
#   	// was reconstructed (tog et some measure of the confidence):
#   	// hco.method

#   end of the while loop

    if command == 2 :
        print "calibrating detector... "
        sorter.do_calibration()
        print "ok - after do_calibration"
        sfco = hexanode.py_scalefactors_calibration_class(sorter)
        if sfco :
            print "Good calibration factors are:\n  f_U =%f\n  f_V =%f\n  f_W =%f\n  Offset on layer W=%f\n"%\
                  (2*sorter.fu, 2*sfco.best_fv, 2*sfco.best_fw, sfco.best_w_offset)


    if command == 3 : # generate and print correction tables for sum- and position-correction
        print "creating calibration tables..."
        status = hexanode.py_create_calibration_tables(FNAME_CALIBRATION_TABLE, sorter)
        print "finished creating calibration tables: %s status %s" % (FNAME_CALIBRATION_TABLE, status)


    print "consumed time (sec) = %.6f\n" % (time() - t_sec)

    if sorter is not None : del sorter

#------------------------------

if __name__ == "__main__" :
    tname = sys.argv[1] if len(sys.argv) > 1 else '1'
    print 50*'_', '\nTest %s:' % tname

    py_sort()
    plot_histograms(prefix=PLOT_PREFIX, do_save=True, hwin_x0y0=(0,0))
    show()
    sys.exit(0)

#------------------------------
