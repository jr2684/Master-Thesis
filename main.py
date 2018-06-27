# -*- coding: utf-8 -*-

from scipy.spatial import distance
import RBF
import functions_for_coefficients_calculations as coe
import interpolant_points
import gp_creation as gp
import time
import numpy as np
import random

refinement_parameter = gp.kosher_grid_values(3)
#RBF_type = 'IQ'
epsilon = 1
possible_RBF = ['IMQ','gau','IQ','MQ'] 
possible_epsilon = [.977,.989,1,1.012,1.025]
grid_start_value = -1
grid_end_value = 1

def test_function(x,y):
    return np.cos(y)*np.sin(x) + np.sin(x)*np.cos(y)
    
for RBF_type in possible_RBF:
    start_time = time.time()
    for epsilon in possible_epsilon:
        start_time = time.time()
        #now we want to solve for the coefficients of the interpolation function
        coefficients_refined = np.linalg.solve(coe.coefficient_matrix(gp.grid_refinement(grid_start_value,grid_end_value,refinement_parameter),RBF_type),coe.analytic_values(grid_start_value,grid_end_value,refinement_parameter))   

        #In this cell we creat the points that are going to be used to test the RBF method against the
        #analytical function. In the instance below, we take a random sample of 1000 points that are in the
        #grid (-1,1)x(-1,1). Depending upon the geometry of the problem, or the conversion to a different coordinate 
        # system, we will use a different set of points
        num_of_data_points = 1000
        arbitrary_values_x = [round(random.uniform(-1,1),3) for i in range(num_of_data_points)]
        arbitrary_values_y = [round(random.uniform(-1,1),3) for i in range(num_of_data_points)]
        for_graph = [[arbitrary_values_x[i],arbitrary_values_y[i]] for i in range(num_of_data_points)]


        #here we will test the interpolation function with the coefficient values calculated above
        error_values = []
        for i in range(len(arbitrary_values_x)):
            value = 0
            for j in range(len(gp.grid_refinement(-1,1,refinement_parameter))):
                value = value + coefficients_refined[j]*RBF.calculate(RBF_type,epsilon,distance.euclidean((arbitrary_values_x[i],arbitrary_values_y[i]), gp.grid_refinement(-1,1,refinement_parameter)[j]))
            error_values.append(abs(test_function(arbitrary_values_x[i],arbitrary_values_y[i])-value))

        print(sum(error_values))
        print("--- %s seconds ---" % (time.time() - start_time))
        #RBF.Produce_error_plot(RBF_type,epsilon,arbitrary_values_x,arbitrary_values_y,error_values)

