

#import numpy as np

import hexanode

def test1():
    print 'hexanode.test1'

def test2():
    print 'hexanode.test2 call hexanode_ext.met1()'
    hexanode.met1()  # cython

