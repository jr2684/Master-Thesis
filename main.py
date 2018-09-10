# -*- coding: utf-8 -*-

from scipy.spatial import distance
from RBF import rbf
import functions_for_coefficients_calculations as coe
import gp_creation as gp
import trial_functions as tf
import v_c_transformations as vc
import time
import numpy as np
import random
import functools
random.seed(9001)

case = 2
produce_plots = 0
refinement_parameter = gp.kosher_grid_values(3)
possible_RBF = ['IMQ','gau','IQ','MQ'] 
possible_epsilon = np.linspace(0,1,1000)
grid_start_value = -1
grid_end_value = 1

#Case 0 is the global case. Used to take a preliminary look at the global error
# and see if the method was working as expected   
if case == 0:
    for RBF_type in possible_RBF:
        start_time = time.time()
        for epsilon in possible_epsilon:
            a = rbf(RBF_type,epsilon)
            start_time = time.time()
            #now we want to solve for the coefficients of the interpolation function
            A = coe.coefficient_matrix(gp.grid_refinement(grid_start_value,grid_end_value,refinement_parameter),RBF_type)
            y = coe.analytic_values(grid_start_value,grid_end_value,refinement_parameter)
            L,D = coe.ldl_decomp(A)
            forward = np.linalg.solve(L,y)
            coefficients_refined = np.linalg.solve(D*L.H,forward)
            #coefficients_refined = np.linalg.solve(coe.coefficient_matrix(gp.grid_refinement(grid_start_value,grid_end_value,refinement_parameter),RBF_type),coe.analytic_values(grid_start_value,grid_end_value,refinement_parameter))   
    
            #In this cell we creat the points that are going to be used to test the RBF method against the
            #analytical function. In the instance below, we take a random sample of 1000 points that are in the
            #grid (-1,1)x(-1,1). Depending upon the geometry of the problem, or the conversion to a different coordinate 
            # system, we will use a different set of points
            num_of_data_points = 1000
            arbitrary_values_x = [round(random.uniform(grid_start_value,grid_end_value),3) for i in range(num_of_data_points)]
            arbitrary_values_y = [round(random.uniform(grid_start_value,grid_end_value),3) for i in range(num_of_data_points)]
            for_graph = [[arbitrary_values_x[i],arbitrary_values_y[i]] for i in range(num_of_data_points)]
    
    
            #here we will test the interpolation function with the coefficient values calculated above
            error_values = []
            for i in range(len(arbitrary_values_x)):
                value = 0
                for j in range(len(gp.grid_refinement(grid_start_value,grid_end_value,refinement_parameter))):
                    value = value + coefficients_refined[j]*a.calculate(distance.euclidean((arbitrary_values_x[i],arbitrary_values_y[i]), gp.grid_refinement(-1,1,refinement_parameter)[j]))
                error_values.append(abs(tf.scalar_case(arbitrary_values_x[i],arbitrary_values_y[i])-value))
    
            print(sum(error_values))
            print("--- %s seconds ---" % (time.time() - start_time))
            if(produce_plots == 1):
                a.Produce_error_plot(arbitrary_values_x,arbitrary_values_y,error_values,for_graph)
#case 1 is for the local case, where a reduced number of neighbors are considered for the interpolation
elif case == 1:
    a = rbf('gau',1)
    for grid_refinement_level in [10,20,40,80,260,320]:
        num_of_data_points = (2*grid_refinement_level)**2
        interpolation_grid = gp.random_grid_normal(num_of_data_points,-np.pi,np.pi,-np.pi/2,np.pi/2)
        X = gp.grid_values(grid_refinement_level)  #desired grid
        cartesian_X = [vc.gc_to_cc(i[0],i[1]) for i in X]
        cartesian_inter_grid = [vc.gc_to_cc(i[0],i[1]) for i in interpolation_grid]
        neighbors = 9
        #for neighbors in [3,9,15]:
        start_time = time.time()
        error_values = []
        for coordinate in cartesian_inter_grid:
            points = gp.nearest_neighbors(coordinate,cartesian_X,neighbors)
            solutions = [tf.vector_case_1(i[0],i[1],1) for i in points]
                           
            A = coe.coefficient_matrix(points,a.function,a.epsilon)
            L,D = coe.ldl_decomp(A)
            forward = np.linalg.solve(L,solutions)
            coefficients_refined = np.linalg.solve(D*L.H,forward)
            sum_elements = [coefficients_refined[j]*a.calculate(distance.euclidean(coordinate,points[j])) for j in range(len(coefficients_refined))]
            interpolation_value = functools.reduce(lambda a,x:a +x,sum_elements,0)
            #coefficients_refined = np.linalg.solve(coe.coefficient_matrix(points,a.function,a.epsilon),solutions)
            error_values.append(abs(tf.vector_case_1(coordinate[0],coordinate[1],1)-interpolation_value))
            #print(str(neighbors) + ' nearest neighbors: ' + str(abs(tf.scalar_case(coordinate[0],coordinate[1])-value)))
            #print(sum(sorted(error_values,reverse = True)[1:10]),statistics.median(error_values))
        print("--- %s seconds ---" % (time.time() - start_time),sum(error_values),grid_refinement_level,neighbors)
# the following code is working for a 2D vector. 
elif case == 2:
    error_values = []
    for function_type in ['gau','MQ','IMQ','IQ']:
        print(function_type)
        a = rbf(function_type,1)
        for grid_refinement_level in [2,4,8,16,32,64]:
            num_of_data_points = (2*grid_refinement_level)**2
            interpolation_grid = gp.random_grid_normal(num_of_data_points,-1,1,-1,1)
            X = gp.grid_refinement_square_refactored(-1,1,grid_refinement_level)  #desired grid
            #interpolation_grid = midpoint_grid(X,grid_refinement_level,-1,-1)
            #neighbors = 15
            for neighbors in [3,9,15]:
                start_time = time.time()
        
        #print(grid_refinement_level)
                for coordinate in interpolation_grid:
                    temp = []
                    points = gp.nearest_neighbors(coordinate,X,neighbors) #find the nearest neighbors
            solutions_u = [tf.vector_case_0(i[0],i[1],1,1) for i in points] #determine the rhs of the linear system
            solutions_v = [tf.vector_case_0(i[0],i[1],2,1) for i in points]
            A = coe.coefficient_matrix_optimized(points,a.function,a.epsilon) #compute the interpolation matrix
                #P, L, U = scipy.linalg.lu(A) #matrix decomp
    
                #forward_u = np.linalg.solve(L,solutions_u)
                #forward_v = np.linalg.solve(L,solutions_v)
    
            coefficients_refined_u = np.linalg.solve(A,solutions_u)
            coefficients_refined_v = np.linalg.solve(A,solutions_v)
    
            sum_elements_u = [coefficients_refined_u[j]*a.calculate(distance.euclidean(coordinate,points[j])) for j in range(len(coefficients_refined_u))]
            sum_elements_v = [coefficients_refined_v[j]*a.calculate(distance.euclidean(coordinate,points[j])) for j in range(len(coefficients_refined_v))]
            interpolation_value_u = functools.reduce(lambda a,x:a +x,sum_elements_u,0)
            interpolation_value_v = functools.reduce(lambda a,x:a +x,sum_elements_v,0)
    
                        #coefficients_refined = np.linalg.solve(coe.coefficient_matrix(points,a.function,a.epsilon),solutions)
                        #must transform back to (u,v)
    
            error_u,error_v = abs(interpolation_value_u -tf.vector_case_0(coordinate[0],coordinate[1],1,1)),abs(interpolation_value_v -tf.vector_case_0(coordinate[0],coordinate[1],2,1))
            temp.append([error_u,error_v])
            print(str(time.time() - start_time) + ' seconds')
            error_values.append(temp)
        print()
