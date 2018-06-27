# -*- coding: utf-8 -*-
import random
import numpy as np

"""In this cell we creat the points that are going to be used to test the RBF method against the
analytical function. In the instance below, we take a random sample of 1000 points that are in the
grid (-1,1)x(-1,1). Depending upon the geometry of the problem, or the conversion to a different coordinate 
system, we will use a different set of points
    Input: - start value: int or float. a in [a,b]
           - end value: int or float. b in [a,b]
           - num_of_points: int. Total number of points.
    Output: a list of random cartesian coordinates
"""

def random_point_generator(start_value,end_value,num_of_points):
    return np.array([[round(random.uniform(start_value,end_value),3),round(random.uniform(start_value,end_value),3)] for i in range(num_of_points)] )