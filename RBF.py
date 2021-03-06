# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import numpy.ma as ma

class rbf:
    
    def __init__(self,func = None, eps = None):
        if func == None and eps == None:
            self.function = 'gau'
            self.epsilon = 1
        elif eps == None and func != None:
            self.epsilon = 1
            self.function = func
            
        else:
            self.function = func
            self.epsilon = eps

    '''
    All radial basis functions will be combined into a single method. The desired RBF can be chosen by giving the 
    keyword that is associated with that particular function.
            - 'gau' : gaussian 
            - 'MQ' : multiquadratic
            - 'IQ' : inverse quadratic
            - 'IMQ' : Inverse Multiquadratic
            - 'TPS' : Thin plate spline 
        inputs: 
            - epsilon-float parameter. Called the shape parameter and is used to determine how 'peaked' the RBF
            is. THe closer to 0, the sharper the function is.
            - func- string. Determines the type of RBF to be returned
            - a - float. The independent input variable of the problem
        output: 
            float value
    '''

    def return_func(self):
        print(self.function)
        
    def return_eps(self):
        print(self.epsilon)
        
    def calculate(self,r):
        if self.function == 'gau':
            return np.exp(-((r/self.epsilon)**2))
        elif self.function == 'MQ':
            return np.sqrt(1+(r/self.epsilon)**2)
        elif self.function == 'IQ':
            return 1/(1 +(r/self.epsilon)**2)
        elif self.function == 'IMQ':
            return 1/np.sqrt(1 + (r/self.epsilon)**2)
        elif self.function == 'TPS':
            return ((r/self.epsilon)**2)*np.log(r/self.epsilon) if r > 0 else 0


    """
        This method produces a continuous heat map of the error.
            Input:
                x_values - array of ints
                y_values - array of ints
                error - array of ints
            Output:
                Graph
    
    """
def Produce_error_plot(ftype,evalue,x_values,y_values,error,interpolation_grid):
    x_list = np.array(x_values)
    y_list = np.array(y_values)
    z_list = np.array(error)

    grid_x, grid_y = np.mgrid[x_list.min():x_list.max():1000j, y_list.min():y_list.max():1000j]
    method = 'cubic'
    plt.figure()
    grid_z = scipy.interpolate.griddata(interpolation_grid,z_list,(grid_x, grid_y), method=method)
    # [pcolormesh with missing values?](https://stackoverflow.com/a/31687006/395857)
    plt.pcolormesh(grid_x, grid_y, ma.masked_invalid(grid_z), cmap='RdBu', vmin=np.nanmin(grid_z), vmax=np.nanmax(grid_z))
    plt.title(str('Error'))
    #plt.title('{0} interpolation'.format(method))
    plt.colorbar()
    plt.show()
    #plt.savefig('heatmap_interpolation_{0}.png'.format(method), dpi=300)
    #plt.clf()
    #plt.close()
    