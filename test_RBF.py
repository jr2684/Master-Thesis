# -*- coding: utf-8 -*-
from RBF import rbf

a = rbf('gau')
b = rbf()
c = rbf('IQ',.5)

print(a.calculate(1))