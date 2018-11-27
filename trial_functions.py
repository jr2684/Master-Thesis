import numpy as np


def scalar_case(x,y):
    return np.cos(y)*np.sin(x) + np.sin(x)*np.cos(y)


def vector_case_0(x,y,component,k):
    if component == 1:
        #u component
        return np.cos(k*np.pi*(x-.25))*np.sin(k*np.pi*(y-.25))
    if component == 2:
        #v component
        return np.sin(k*np.pi*(x-.25))*np.cos(k*np.pi*(y-.25))    
    
    
'''parameters for the vector case
    k is the 
    R is the wave number
    a is the radius of the earth

'''

def vector_case_1(u1,u2,component):
    k = 7.484e-6
    a = 6.37122e6
    w = k
    R = 4
    if component == 1:
        #represents u component
        return a*w*(np.cos(u1))+a*k*(((np.cos(u1)))**(R-1))*((R*(np.sin(u1))**2)-np.cos(u1)**2)*np.cos(R*float(u2))
    elif component == 2:
        #represents v component
        return -a*k*R*(np.cos(u1)**(R-1))*np.sin(u1)*np.sin(R*u2)