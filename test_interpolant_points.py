# -*- coding: utf-8 -*-
import unittest
from interpolant_points import random_point_generator

class InterpolantTestCase(unittest.TestCase):
    
    def test_1D(self):
       self.assertIsInstance(random_point_generator(-1,1,1000,1)[0],float)
    def test_2d(self):
        self.assertIsInstance(random_point_generator(-1,1,1000,2)[0],list)
    
    
    
if __name__ == '__main__':
    unittest.main()
    
