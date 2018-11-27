# -*- coding: utf-8 -*-

from scipy.spatial import distance
from RBF import rbf
from numpy.linalg import LinAlgError
import functions_for_coefficients_calculations as coe
import gp_creation as gp
import trial_functions as tf
import v_c_transformations as vc
import time
import numpy as np
import random
import functools
import meshzoo

random.seed(9001)

case = 3
produce_plots = 0
refinement_parameter = gp.kosher_grid_values(3)
possible_RBF = ['gau'] 
num_neighbors = [9]
epsilonlist = [1]
grid_start_value = -1
grid_end_value = 1

scale_factor = []
scale_factor.extend(np.linspace(0.001,0.01,10))
scale_factor.extend(np.linspace(0.01,.1,10))
scale_factor.extend(np.linspace(.1,1,10))
scale_factor.extend(np.linspace(1,1.1,10))
scale_factor_list = sorted(list(set(scale_factor)))
    

"""
Case 1 is the global case. Used to take a preliminary look at the global error
and see if the method was working as expected  
""" 
if case == 1:
    for RBF_type in possible_RBF:
        start_time = time.time()
        epsilon = 1
        a = rbf(RBF_type,epsilon)
        A = coe.coefficient_matrix(gp.grid_refinement_square(grid_start_value,grid_end_value,refinement_parameter),RBF_type,epsilon)
        y = coe.analytic_values(grid_start_value,grid_end_value,refinement_parameter)
        #L,D = coe.ldl_decomp(A)
        #forward = np.linalg.solve(L,y)
        #coefficients_refined = np.linalg.solve(D*L.H,forward)
        coefficients_refined = np.linalg.solve(A,y)   
    
        num_of_data_points = 1000
        for_graph = gp.random_grid_normal(num_of_data_points,grid_start_value,grid_end_value,grid_start_value,grid_end_value,3)
    
        'here we will test the interpolation function with the coefficient values calculated above'
        error_values = []
        for i in range(len(for_graph)):
            value = 0
            for j in range(len(gp.grid_refinement_square(grid_start_value,grid_end_value,refinement_parameter))):
                value = value + coefficients_refined[j]*a.calculate(distance.euclidean((for_graph[i][0],for_graph[i][1]), gp.grid_refinement_square(-1,1,refinement_parameter)[j]))
            error_values.append(abs(tf.scalar_case(for_graph[i][0],for_graph[i][1])-value))
    
        print(sum(error_values))
        print("--- %s seconds ---" % (time.time() - start_time))
        if(produce_plots == 1):
            x_values = [i[0] for i in for_graph]
            y_values = [i[1] for i in for_graph]
            a.Produce_error_plot(x_values,y_values,error_values,for_graph)
    """
    the following code is working for a 2D vector in a cartesian plane.
    """
elif case == 2:
    error_values = []
    for function_type in possible_RBF:
        print(function_type)
        a = rbf(function_type,1)
        for grid_refinement_level in [2,4,8,16,32]: 
            print(grid_refinement_level)
            num_of_data_points = (2*grid_refinement_level)**2
            X = gp.grid_refinement_square(-1,1,grid_refinement_level)  #desired grid
            fill_distance = min(2*gp.midpoints(X)[0],2*gp.midpoints(X)[1])
            interpolation_grid = gp.random_grid_normal(num_of_data_points,-1,1,-1,1,len(str(fill_distance))-2)
            for neighbors in num_neighbors:
                for epsilon in scale_factor_list:
                    temp = []
                    for coordinate in interpolation_grid:
                        points = gp.nearest_neighbors(coordinate,X,neighbors) #find the nearest neighbors
                        solutions_u = [tf.vector_case_0(i[0],i[1],1,1) for i in points] #determine the rhs of the linear system
                        solutions_v = [tf.vector_case_0(i[0],i[1],2,1) for i in points]
                        A = coe.coefficient_matrix(points,a.function,a.epsilon) #compute the interpolation matrix
                        #P, L, U = scipy.linalg.lu(A) #matrix decomp
        
                        #forward_u = np.linalg.solve(L,solutions_u)
                        #forward_v = np.linalg.solve(L,solutions_v)
                        try: 
                            A_norm = np.linalg.norm(A)
                            A_inverse = np.linalg.inv(A)
                            A_inverse_norm = np.linalg.norm(A_inverse)
                            condition_A = A_norm * A_inverse_norm
                            coefficients_refined_u = np.linalg.solve(A,solutions_u)
                            coefficients_refined_v = np.linalg.solve(A,solutions_v)
        
                            sum_elements_u = [coefficients_refined_u[j]*a.calculate((distance.euclidean(coordinate,points[j]))) for j in range(len(coefficients_refined_u))]
                            sum_elements_v = [coefficients_refined_v[j]*a.calculate((distance.euclidean(coordinate,points[j]))) for j in range(len(coefficients_refined_v))]
                            interpolation_value_u = functools.reduce(lambda a,x:a +x,sum_elements_u,0)
                            interpolation_value_v = functools.reduce(lambda a,x:a +x,sum_elements_v,0)
                        except LinAlgError:
                            print("Danger Will Robinson!")
                        else:
                            error_u,error_v = abs(interpolation_value_u -tf.vector_case_0(coordinate[0],coordinate[1],1,1)),abs(interpolation_value_v -tf.vector_case_0(coordinate[0],coordinate[1],2,1))
                            error_u_2,error_v_2 = (abs(interpolation_value_u -tf.vector_case_0(coordinate[0],coordinate[1],1,1))**2),(abs(interpolation_value_v -tf.vector_case_0(coordinate[0],coordinate[1],2,1))**2)
                            temp.append([function_type,neighbors,epsilon,error_u,error_v,error_u_2,error_v_2,condition_A])
                    error_values.append(temp)
            print('scaling parameters: ',coe.error_calculations(error_values)[0])
            print()
            print('rmse(u): ',coe.error_calculations(error_values)[3])
            print()
            print('linfty(u): ',coe.error_calculations(error_values)[5])
            print()
            print('condition number(max): ',coe.error_calculations(error_values)[7])
            print()
                    
    

    """
    This code works for taking a 2D geographical system,converting it to cartesian for interpolation,
    and the back to a geographical system
    """      
elif case == 3:         
    error_values = []
    for grid_refinement_level in [4]:
        grid_points, cells = meshzoo.iso_sphere(grid_refinement_level)
        X = [vc.cc_to_gc(i[0],i[1],i[2]) for i in grid_points]
        interpolation_grid = gp.random_grid_normal(len(grid_points),-np.pi,np.pi,(-np.pi/2) +.3,(np.pi/2) -.3,3)

        for neighbors in num_neighbors:
            for epsilon in scale_factor_list:
                for function in possible_RBF:
                    a = rbf(function,epsilon)
                    temp = []
                    print(grid_refinement_level)
                    for coordinate in interpolation_grid:
                        cartesian_coordinate = vc.gc_to_cc(coordinate[0],coordinate[1])
                        points = gp.nearest_neighbors(coordinate,X,neighbors)
                        cartesian_points = [vc.gc_to_cc(i[0],i[1]) for i in points]
                        solutions_u = [tf.vector_case_1(i[0],i[1],1) for i in points]
                        solutions_v = [tf.vector_case_1(i[0],i[1],2) for i in points]
                        cartesian_solutions = [vc.gc_to_cc_vec(solutions_u[i],solutions_v[i],points[i][0],points[i][1]) for i in range(len(points))]
                        cartesian_solutions_x = [i[0] for i in cartesian_solutions]
                        cartesian_solutions_y = [i[1] for i in cartesian_solutions]
                        cartesian_solutions_z = [i[2] for i in cartesian_solutions]
        
                        A = coe.coefficient_matrix_3d(cartesian_points,a.function,a.epsilon)
                        #P, L, U = scipy.linalg.lu(A)
        
                        #forward_x = np.linalg.solve(L,cartesian_solutions_x)
                        #forward_y = np.linalg.solve(L,cartesian_solutions_y)
                        #forward_z = np.linalg.solve(L,cartesian_solutions_z)
                        try:
                            A_norm = np.linalg.norm(A)
                            A_inverse = np.linalg.inv(A)
                            A_inverse_norm = np.linalg.norm(A_inverse)
                            condition_A = A_norm * A_inverse_norm
                            coefficients_refined_x = np.linalg.solve(A,cartesian_solutions_x)
                            coefficients_refined_y = np.linalg.solve(A,cartesian_solutions_y)
                            coefficients_refined_z = np.linalg.solve(A,cartesian_solutions_z)
                        except LinAlgError:
                            print("Danger!")
                            temp.append([function,0,0,0,0,0])
                        else:
                            sum_elements_x = [coefficients_refined_x[j]*a.calculate(vc.arc_length(cartesian_coordinate,cartesian_points[j])) for j in range(len(coefficients_refined_x))]
                            sum_elements_y = [coefficients_refined_y[j]*a.calculate(vc.arc_length(cartesian_coordinate,cartesian_points[j])) for j in range(len(coefficients_refined_y))]
                            sum_elements_z = [coefficients_refined_z[j]*a.calculate(vc.arc_length(cartesian_coordinate,cartesian_points[j])) for j in range(len(coefficients_refined_z))]
                            interpolation_value_x = functools.reduce(lambda a,x:a +x,sum_elements_x,0)
                            interpolation_value_y = functools.reduce(lambda a,x:a +x,sum_elements_y,0)
                            interpolation_value_z = functools.reduce(lambda a,x:a +x,sum_elements_z,0)
    
                            #must transform back to (u,v)
                            reconstructed = vc.cc_to_gc_vec(coordinate[0],coordinate[1],interpolation_value_x,interpolation_value_y,interpolation_value_z)
                            error_u,error_v = abs(reconstructed[0] -tf.vector_case_1(coordinate[0],coordinate[1],1)),abs(reconstructed[1] -tf.vector_case_1(coordinate[0],coordinate[1],2))
                            error_u_2,error_v_2 = (abs(reconstructed[0] -tf.vector_case_1(coordinate[0],coordinate[1],1))**2),(abs(reconstructed[1] -tf.vector_case_1(coordinate[0],coordinate[1],2))**2)
                            temp.append([function,neighbors,epsilon,error_u,error_v,error_u_2,error_v_2,condition_A])
                    error_values.append(temp)
        
                print('scaling parameters: ',coe.error_calculations(error_values)[0])
                print()
                print('rmse(u): ',coe.error_calculations(error_values)[3])
                print()
                print('linfty(u): ',coe.error_calculations(error_values)[5])
                print()
                print('condition number(max): ',coe.error_calculations(error_values)[7])
                print()
                    
          