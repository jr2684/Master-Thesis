# -*- coding: utf-8 -*-
import unittest
import grid_creation


class gridtestcase(unittest.TestCase):
    "Tests for grid_creation.py."
    
    def test_kosher_grid_values_0(self):
        self.assertNotEqual(grid_creation.kosher_grid_values(0),0)
    def test_kosher_grid_values_3(self):
        self.assertNotEqual(grid_creation.kosher_grid_values(3),3)
    def test_kosher_grid_values_mult3(self):
        self.assertNotEqual(grid_creation.kosher_grid_values(6),7)
        

if __name__ == '__main__':
    unittest.main()