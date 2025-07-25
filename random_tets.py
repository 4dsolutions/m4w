#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 12:59:17 2024

@author: kirbyurner

Generate random tetrahedron in the IVM
"""

from qrays import Qvector
from itertools import permutations
from sympy import Integer
from tetravolume import qvolume, Tetrahedron
from random import choice

zero = Integer(0)
one  = Integer(1)
two  = Integer(2)

ivms = [Qvector(v) for v in 
        set(list(permutations((two, one, one, zero))))]

def gen_tet(n):
    "Generate random IVM tetrahedron"
    
    # all turtles start at the same center
    ball0 = Qvector((zero, zero, zero, zero))
    ball1 = Qvector((zero, zero, zero, zero))
    ball2 = Qvector((zero, zero, zero, zero))
    ball3 = Qvector((zero, zero, zero, zero))
    
    for _ in range(n):
        # hop randomly n times in any of 12 ivm directions
        ball0 += choice(ivms)
        ball1 += choice(ivms)
        ball2 += choice(ivms)
        ball3 += choice(ivms)
    
    print(ball0)
    print(ball1)
    print(ball2)
    print(ball3)
    
    a = (ball0 - ball1).length()
    b = (ball0 - ball2).length()
    c = (ball0 - ball3).length()
    d = (ball1 - ball2).length()
    e = (ball2 - ball3).length()
    f = (ball3 - ball1).length()
    
    T = Tetrahedron(a,b,c,d,e,f)
    
    lengths = T.edges(True, 20)
    for edge, length in T.edges().items():
        print("{}: {:14} {}".format(edge, str(length), lengths[edge] ))
        
    print()
    print("Volume:", T.ivm_volume())
    print("Volume:", T.ivm_volume(True, 30))
    print("Volume:", qvolume(ball0, ball1, ball2, ball3))

if __name__ == "__main__":
    gen_tet(1000)
    
