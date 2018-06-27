# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial import distance
import RBF

'''
    The following method calculated the A matrix. In particular, we nest 2 different processes together. The first is the 
    method matrix_distances. This method takes calculates the euclidean distance between each point. The loop that follows
    then evaluates the RBF at that distance. It is this matrix of values that is then returned. An important side note is that
    the RBF function is hard coded into this method. One must then change the method for various parameter and model values
    
        input: grid - multidim. array
               function- str. Gives the type of function that will be used for the calculations
        output: grid - multidim. array
'''
def coefficient_matrix(grid,function):
    def matrix_distances(grid):
        matrix_distances_refined = []
        for j in range(len(grid)):
            temp = np.array([])
            temp = np.append(temp,distance.euclidean(grid[j],grid[0]))
            for i in range(len(grid)-1):
                temp = np.append(temp,distance.euclidean(grid[j], grid[i+1]))
            matrix_distances_refined.append(temp)
        return matrix_distances_refined
    matrix_evaluated_refined = []
    for i in matrix_distances(grid):
        matrix_evaluated_refined.append([RBF.calculate(function,1,point) for point in i])
    return matrix_evaluated_refined
    
    
"""
Here the test function is defined, along with a loop that produces the solution values f. We use these values
to test the accuracy of our RBF method. In general we won't be using data that has been created with a 
test function. Instead, it will come from some exterior source or from
another function call. 
    Input: - start_value: int. This is the minimum value that our grid will take.
           - end_value : int. This is the maximum value that our grid will take. 
           - step size : int. This is a step size we will use in our grid
    Output:
            -vector_values_refined: array of floats. This is an array of the analytic solutions. 

"""
def analytic_values(start_value,end_value,step_size):
    def test_function(x,y):
        return np.cos(y)*np.sin(x) + np.sin(x)*np.cos(y)
    xx = np.linspace(start_value,end_value,2*step_size + 1)
    yy = np.linspace(start_value,end_value,2*step_size + 1)
    XX,YY = np.meshgrid(xx,yy)
    zz = test_function(XX,YY)
    vector_values_refined = []
    for i in zz:
        for j in i:
            vector_values_refined.append(j)
    return vector_values_refined