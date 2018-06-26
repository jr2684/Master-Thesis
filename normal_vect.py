# -*- coding: utf-8 -*-
from vectors import Vector

import numpy as np

def perpendicular_vector(v):
    r""" Finds an arbitrary perpendicular vector to *v*."""
    # for two vectors (x, y, z) and (a, b, c) to be perpendicular,
    # the following equation has to be fulfilled
    #     0 = ax + by + cz

    # x = y = z = 0 is not an acceptable solution
    if v.x == v.y == v.z == 0:
        raise ValueError('zero-vector')

    # If one dimension is zero, this can be solved by setting that to
    # non-zero and the others to zero. Example: (4, 2, 0) lies in the
    # x-y-Plane, so (0, 0, 1) is orthogonal to the plane.
    if v.x == 0:
        return Vector(1, 0, 0)
    if v.y == 0:
        return Vector(0, 1, 0)
    if v.z == 0:
        return Vector(0, 0, 1)

    # arbitrarily set a = b = 1
    # then the equation simplifies to
    #     c = -(x + y)/z
    return Vector(1, 1, -1.0 * (v.x + v.y) / v.z)