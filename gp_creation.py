# -*- coding: utf-8 -*-
import numpy as np
from sklearn.neighbors import KDTree


"""
grid_refinement(left_origin,right_origin,refinement_level):

The function outputs an array consisting of all the gridpoints of a grid of dimensions left_origin x right_origin.
We can iteratively change the size of the grid via the variable refinement_level.
        Inputs:
            left_origin:int
            right_origin:int
            refinement_level:int
        Output:
            An array with entries of the form (i,j)
            
        Notes:
            
    
    """


def grid_refinement(left_origin,right_origin,refinement_level):
    #test_a = -1 #left_origin. These two values are what were used as a test 
    #test_b = 1 #right_origin
    refinement_level_tested = kosher_grid_values(refinement_level)
    test_n = 2*refinement_level_tested +1
    test_step = 1/refinement_level_tested
    test_x = np.linspace(left_origin,right_origin,test_n)
    test_y = np.linspace(left_origin,right_origin,test_n)
    test_array = []
    for i in np.arange(right_origin,(left_origin - test_step),-test_step):
        for j in np.arange(left_origin,right_origin + test_step,test_step):
            test_array.append((j,i))
    return(test_array)



"""
kosher_grid_values(integer):

    
Due to the nature of repeating decimals, we should only consider partitions that consist of multiples of 2
or 5 exclusively. In order to preserve equal graph spacings, we consider only a selected list of refinement 
values. If a value is entered that does not preserve the spacing, then the next largest value is taken in 
its place.Otherwise we introduce an error into our calculation. This isn't a fully tested feature and 
can easily be broken.
        Input:
            integer: int
            
        Output:
            an integer

    
"""


def kosher_grid_values(integer):
    kosher_refinement_values = [1,2,4,5,8,10,16,20,25,32,40,64,80,125,128,260,256,320,512,625,640,1280,2560,3125]
    if(integer not in kosher_refinement_values):
        i = 0
        while(integer > kosher_refinement_values[i]):
            i +=1
        integer = kosher_refinement_values[i]
    return integer
    
    
"""This method will find the n nearest neigbors of a point on a grid.
    Input: 
        point - (x,y) coordinate
        grid - numerical mesh that is being used for calculations
        num_of_neighbors - desired number of neighbors to be returned
        ToDo:(Write statement so that
        this value cannot exceed the dimensions of the grid)
    Output:
        tuple of length num_of_neighbors
"""
def nearest_neighbors(point,grid,num_of_neighbors):
    tree = KDTree(grid, leaf_size = 20)
    dist, ind = tree.query([point],k = num_of_neighbors)
    indexes = [grid[j] for j in ind[0]]
    return indexes