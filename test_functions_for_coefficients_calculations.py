# -*- coding: utf-8 -*-
import unittest
from functions_for_coefficients_calculations import ldl_decomp
from functions_for_coefficients_calculations import clean_machine_ep
import numpy as np
#from np.linalg import cholesky,LinAlgError

class Testcoefficientfunctions(unittest.TestCase):
    def test_hermitian(self):
        A = [[1,3],[2,4]]
        self.assertIsNone(ldl_decomp(A)[0])
    def test_square(self):
        A = [[1,2],[1,1],[1,2]]
        self.assertIsNone(ldl_decomp(A)[0])
    def test_positive_definite(self):
        A = [[-1,1],[2,-3]]
        with self.assertRaises(np.linalg.LinAlgError):
            np.linalg.cholesky(A)
    def test_nonsingular(self):
        A = np.matrix([[1,1],[1,1]])
        self.assertFalse(A.shape[0] == A.shape[1] and np.linalg.matrix_rank(A) == A.shape[0])
    def test_machine_epsilon_scalar(self):
        self.assertEqual(3.2e-16,clean_machine_ep(3.2e-16))
        self.assertNotEqual(1.2e-16,clean_machine_ep(1.2e-16))
    def test_machine_epsilon_list(self):
        self.assertEqual([-2.4, 0, 0],clean_machine_ep([-2.4,-2.4e-17,0]))
        self.assertEqual([[-2.4, 2, 1.53, 3], [0, 2, 0, 0]],clean_machine_ep([[-2.4,2,1.53,3],[-2.4e-17,2,1.53e-16,0]]))
        
       
            
if __name__ == '__main__':
    unittest.main()