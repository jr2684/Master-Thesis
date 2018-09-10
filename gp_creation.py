# -*- coding: utf-8 -*-
import numpy as np
import random
from sklearn.neighbors import KDTree
import itertools
from operator import itemgetter
from v_c_transformations import gc_to_cc
from functions_for_coefficients_calculations import clean_machine_ep


"""
grid_refinement(left_origin,right_origin,refinement_level):

The function outputs a square array consisting of all the gridpoints of a grid of dimensions left_origin x right_origin.
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
grid_refinement_rectangle(left,right, top, bottom,refinement_level):

The function outputs a rectangular array consisting of all the gridpoints of a grid of dimensions
left x right x top x bottom
We can iteratively change the size of the grid via the variable refinement_level.
        Inputs:
            left:int
            right:int
            top:int
            bottom:int
            refinement_level:int
        Output:
            An array with entries of the form (i,j)
            
        Notes:
        Currently this is only fully functioning for the spherical case(this is because of the unique boundry condition 
of [-pi,pi). In order to make the interval open and exclude pi, some magic had to be performed.)
        
    
    """
def grid_refinement_rectangle_refactored(left,right,top,bottom,refinement_level):
    refinement_level_tested = kosher_grid_values(refinement_level)
    test_n = 2*refinement_level_tested +1
    x = sorted(np.linspace(right,left,test_n,endpoint = False),reverse = False)
    y = np.linspace(bottom,top,test_n)
    grid = sorted(list(itertools.product(x,y)),key = itemgetter(1),reverse = True)
    return grid    
    
"""
Creates a square mesh. Can be seen as a special case of grid_refinement_
rectangle_refactored
"""
def grid_refinement_square_refactored(left_origin,right_origin,refinement_level):
    refinement_level_tested = kosher_grid_values(refinement_level)
    test_n = 2*refinement_level_tested +1
    x = np.linspace(left_origin,right_origin,test_n)
    grid = sorted(list(itertools.product(x,x)),key = itemgetter(1),reverse = True)
    return grid



    
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

"""
The following two methods are used to create the grids for the vector case.These
differ from the scalar case in that the cartesian transformation produces 
nonsingularities in the coefficient matrix. This is avoided by removing the 
extrenuous values before performing the grid transformations
"""

#values_temp creates a rectangular grid that ranges on the x axis from [-pi,pi] and on the y axis from [-pi/2,pi/2].
#The three lines following are used to remove any values involving the poles(they all are the same when transformed)
#The final outputted array is analytic_positions, which consists of the unique grid values that will be used for the rest
#of the calculations
def grid_values(refinement_level):
    values_temp = grid_refinement_rectangle_refactored(-np.pi,np.pi,np.pi/2,-np.pi/2,refinement_level)
    vector_poles = [values_temp[0],values_temp[-1]]
    temp_list = values_temp.copy()
    modified = [edge for edge in temp_list if edge[1] not in [-np.pi/2,np.pi/2]]
    analytic_positions = [vector_poles[0]] + modified + [vector_poles[1]]
    return analytic_positions

def gc_to_cc_grid_values(refinement_level):
    #the geographical values are then transformed into cartesian values.
    transformed_values =[gc_to_cc(i[0],i[1]) for i in grid_values(refinement_level)]
    transformed_clean = clean_machine_ep(transformed_values)
    return transformed_clean
    
"""The following three methods are used to create different types of graphs. So far, the only type of graph that
    has been considered is a square(or rectangular) graph. The methods center_grid and midpoint grid will take a
    normal rectangular(or square) grid as its input and output a grid that consists of either a grid of the midpoints 
    or a grid of thhe centers.The method midpoints is a subroutine used inside of the other two methods. In particular
    these 2 methods will be used as intermediary grids for the vector interpolation. Since a coupler will be used to 
    interpolate between grids of different size, we will simulate these other grids via these two methods.
"""
    
def midpoints(grid):
    distance_horizontal = abs(round((sorted(grid,key = itemgetter(1),reverse = True)[1][0] - sorted(grid,key = itemgetter(1),reverse = True)[0][0]),2))
    distance_vertical = abs(round((sorted(grid,key = itemgetter(0),reverse = True)[1][1] - sorted(grid,key = itemgetter(0),reverse = True)[0][1]),2))
    midpoint_horiz = distance_horizontal/2
    midpoint_vert = distance_vertical/2
    return [midpoint_horiz,midpoint_vert]

def center_grid(grid,refinement,left,down):
    midpoint_distance_horizontal,midpoint_distance_vertical = midpoints(grid)
    test_n = 2*refinement +1
    x1 = [left + (2*i+1)*(midpoint_distance_horizontal) for i in range(test_n+1)]
    y1 = [down + (2*i+1)*(midpoint_distance_vertical) for i in range(test_n+1)]
    data = list(itertools.product(x1,y1))
    return data
    
def midpoint_grid(grid,refinement,left,down):
    midpoint_distance_horizontal,midpoint_distance_vertical = midpoints(grid)
    test_n = 2*refinement +1
    x2 = [left + (2*i+1)*(midpoint_distance_horizontal) for i in range(test_n+1)]
    y2 = [down + (2*i)*(midpoint_distance_vertical) for i in range(test_n+2)]
    data1 = list(itertools.product(x2,y2))
    
    x3 = [left + (2*i)*(midpoint_distance_horizontal) for i in range(test_n+2)]
    y3 = [down + (2*i+1)*(midpoint_distance_vertical) for i in range(test_n+1)]
    data2 = list(itertools.product(x3,y3))
    data = data1 + data2
    return data
    
"""produces a nonuniform grid of num_of_points vertices. These random grids should be the least
   accurate of all the the interpolation grids
"""

def random_grid_normal(num_of_points,left,right,down,up):          
    x = [round(random.uniform(left,right),3) for i in range(num_of_points)]
    y = [round(random.uniform(down,up),3) for i in range(num_of_points)]
    return [[x[i],y[i]] for i in range(num_of_points)]

def random_grid_beta(num_of_points,left,right,down,up,alpha,beta):          
    x = [right*round(random.betavariate(alpha,beta),3) for i in range(int(num_of_points/2))] +[left*round(random.betavariate(alpha,beta),3) for i in range(int(num_of_points/2))]
    y = [up*round(random.betavariate(alpha,beta),3) for i in range(int(num_of_points/2))]+[down*round(random.betavariate(alpha,beta),3) for i in range(int(num_of_points/2))]
    random.shuffle(x)
    random.shuffle(y)
    return [[x[i],y[i]] for i in range(num_of_points)]