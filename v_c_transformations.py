# -*- coding: utf-8 -*-
import numpy as np

#transforms a vector in geographical coordinates to cartesian coordinates
def gc_to_cc_vec(u,v,gc1,gc2):
    u1 = -u*np.sin(gc1) -v*np.sin(gc2)*np.cos(gc1)
    u2 = u*np.cos(gc1) - v*np.sin(gc2)*np.sin(gc1)
    u3 = v*np.cos(gc2)
    return [u1,u2,u3]

#transforms a vector in cartesian coordinates to geographical coordinates
def cc_to_gc_vec(gc1,gc2,u1,u2,u3):
    a1 = -u1*np.sin(gc1) + u2*np.cos(gc1)
    b1 = -u1*np.sin(gc2)*np.cos(gc1) -u2*np.sin(gc2)*np.sin(gc1) + u3*np.cos(gc2)
    return [a1,b1]

#Transforms coordinates from geographical to cartesian
def gc_to_cc(gc1,gc2):
    return [np.cos(gc2)*np.cos(gc1), np.cos(gc2)*np.sin(gc1),np.sin(gc2)]
            