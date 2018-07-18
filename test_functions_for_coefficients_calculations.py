# -*- coding: utf-8 -*-
import unittest
from functions_for_coefficients_calculations import ldl_decomp
import numpy as np
#from np.linalg import cholesky,LinAlgError

class Testcoefficientfunctions(unittest.TestCase):
    def test_hermitian(self):
        A = [[1,3],[2,4]]
        self.assertIsNone(ldl_decomp(A)[0])
    def test_square(self):
        A = [[1,2],[1,1],[1,2]]
        self.assertIsNone(ldl_decomp(A)[0])
    def positive_definite(self):
        A = [[-1,1],[2,-3]]
        with self.assertRaises(np.linalg.LinAlgError):
            np.linalg.cholesky(A)
    def nonsingular(self):
        A = [[1,1],[1,1]]
        with self.assertRaises(np.linalg.LinAlgError):
            np.invert(A)
            
if __name__ == '__main__':
    unittest.main()