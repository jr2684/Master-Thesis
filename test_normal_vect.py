# -*- coding: utf-8 -*-
import unittest
from normal_vect import perpendicular_vector

class NormalTestCase(unittest.TestCase):
    """Tests for normal_vect.py"""
    
    def test_zero_vect(self):
        with self.assertRaises(ValueError):
            perpendicular_vector(Vector(0,0,0))
    def test_z_vect(self):
        self.assertEqual(perpendicular_vector(Vector(2,3,0)).vector,[0,0,1])
    def test_y_vect(self):
        self.assertEqual(perpendicular_vector(Vector(2,0,3)).vector,[0,1,0])
    def test_x_vect(self):
        self.assertEqual(perpendicular_vector(Vector(0,2,3)).vector,[1,0,0])
    def test_arbitrary_vect(self):
        self.assertIsNot(perpendicular_vector(Vector(2,3,4)).vector,[0,0,0])
if __name__ == '__main__':
    unittest.main()