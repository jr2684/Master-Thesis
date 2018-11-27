# -*- coding: utf-8 -*-
import numpy as np
from scipy.spatial import distance
from RBF import rbf

'''
    The following method calculated the A matrix. In particular, we nest 2 different processes together. The first is the 
    method matrix_distances. This method takes calculates the euclidean distance between each point. The loop that follows
    then evaluates the RBF at that distance. It is this matrix of values that is then returned. An important side note is that
    the RBF function is hard coded into this method. One must then change the method for various parameter and model values
    
        input: grid - multidim. array
               function- str. Gives the type of function that will be used for the calculations
        output: grid - multidim. array
'''   
    
def coefficient_matrix(grid,function,epsilon):
    def matrix_distances_hermitian(grid):
        matrix_distances_refined = []
        for j in range(len(grid)):
            temp = np.array([])
            temp = np.append(temp,0)
            for i in range(len(grid)-1):
                if i<j:
                    temp = np.append(temp,0)
                else:
                    temp = np.append(temp,(distance.euclidean(grid[j], grid[i+1])))
            matrix_distances_refined.append(temp)
        return np.array(np.matrix(matrix_distances_refined) +np.matrix(matrix_distances_refined).T)
    matrix_evaluated_refined = []
    a = rbf(function,epsilon)
    for i in matrix_distances_hermitian(grid):
        matrix_evaluated_refined.append([a.calculate(point) for point in i])
    return matrix_evaluated_refined    
    
def coefficient_matrix_3d(grid,function,epsilon):
    def matrix_distances_hermitian(grid):
        matrix_distances_refined = []
        for j in range(len(grid)):
            temp = np.array([])
            temp = np.append(temp,0)
            for i in range(len(grid)-1):
                if i<j:
                    temp = np.append(temp,0)
                else:
                    temp = np.append(temp,(arc_length(grid[j], grid[i+1])))
            matrix_distances_refined.append(temp)
        return np.array(np.matrix(matrix_distances_refined) +np.matrix(matrix_distances_refined).T)
    matrix_evaluated_refined = []
    a = rbf(function,epsilon)
    for i in matrix_distances_hermitian(grid):
        matrix_evaluated_refined.append([a.calculate(point) for point in i])
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

    
"""
Checks if a matrix is square. 
    Input: Matrix
    Output: Boolean value
"""
def isSquare (m):
    return all (len (row) == len (m) for row in m)
    
"""
Since the standard library does not have a modified cholesky decomposition, the
following method is used for its calculation(A = LDL*, where L is a lower unit
triangular matrix, L* is its conjugate, and D is a diagonal matrix
    Input: square positive definite hermitian matrix
    output: L,D
"""

def ldl_decomp(A):
    if isSquare(A) is False:
        print("A must be square!")
        return None,None
    A = np.matrix(A)
    if not (A.H == A).all():
        print("A must be Hermitian!")
        return None, None
    else:
        S = np.diag(np.diag(A))
        Sinv = np.diag(1/np.diag(A))
        D = np.matrix(S.dot(S))
        Lch = np.linalg.cholesky(A)
        L = np.matrix(Lch.dot(Sinv))
        return L, D
        
"""
The following is an intermediary function that will be used to clean up the
coordinate values between transformations. Due to the machine epsilon, there 
are miniscule values (on the order of -16 to -33) that are remenants of the
functions used(This means converting between numerical systems). To clean
things up, every coordinate value is checked for these small values. If the values
present are smaller than the machine epsilon, then the value will be set to
0. In the following method, *args is used, thus allowing for an arbitrary number
of input values. In particular, the method works for a single value,
an array of values, and a matrix of values. It is important to note that 
it will only work for float and int values. If one tries to use strings, 
then the method will break down.
"""

def clean_machine_ep(*args):
    machine_epsilon = np.finfo(float).eps
    for i in args:
        if type(i) == float and np.abs(i) < machine_epsilon :
                return 0
        elif type(i) == list:
            for item in range(len(i)):
                if type(i[item]) == list:
                    for j in range(len(i[item])):
                        if np.abs(i[item][j]) < machine_epsilon:
                            i[item][j] = 0
                else:
                    if np.abs(i[item]) < machine_epsilon:
                        i[item] = 0
    return i

"""
This method is used in place of the euclidean distance in the 3D case. Since we are finding the distance 
on a sphere, it is better to be calculating the arc length.
"""    
def arc_length(x1,x2):
    cc = np.dot(x1,x2)
    if cc > 1:
        cc = 1
    if cc < -1:
        cc = -1
    return np.arccos(cc)   
 
    
'''
This method is used to calculate the error. A number of different norms are calculated,
along with giving the condition of the interpolation matrix.
'''    
def error_calculations(array):
    error_linfty_u = []
    error_linfty_v = []
    error_rmse_u = []
    error_rmse_v = []
    error_l2_u = []
    error_l2_v = []
    condition_max = []
    condition_min =[]
    scale_values = []
    for i in array:
        scale_values.append(i[0][2])
        error_l2_u.append(np.sqrt(sum([j[3] for j in i])))
        error_l2_v.append(np.sqrt(sum([j[4] for j in i])))
        error_linfty_u.append(max([j[3] for j in i]))
        error_linfty_v.append(max([j[4] for j in i]))
        error_rmse_u.append(np.sqrt((1/(len(i)))*sum([j[5] for j in i])))
        error_rmse_v.append(np.sqrt((1/(len(i)))*sum([j[6] for j in i])))
        condition_max.append(max([j[7] for j in i]))
        condition_min.append(min([j[7] for j in i]))
    return[scale_values,error_l2_u,error_l2_v,error_rmse_u,error_rmse_v,error_linfty_u,condition_max,condition_min]
    
    
    