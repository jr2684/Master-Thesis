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

def cc_to_gc(x,y,z):
    u1 = np.arctan2(y,x)
    u2 = np.arcsin(z/np.sqrt(x**2 + y**2 + z**2))
    return [u1,u2]

#measures the arc length between 2 points
def arc_length(x1,x2):
    lenx = np.sqrt(np.dot(x1,x1))
    leny = np.sqrt(np.dot(x2,x2))
    cc = np.dot(x1,x2)/(lenx*leny)
    if cc > 1:
        cc = 1
    if cc < -1:
        cc = -1
    return np.arccos(cc)            