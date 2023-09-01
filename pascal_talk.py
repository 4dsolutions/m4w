#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 09:57:21 2023

@author: kirbyurner
"""

def pascal():
    row = [1]
    while True:
        yield row
        rshift = [0] + row
        lshift = row + [0]
        row = [(a+b) for a,b in zip(rshift, lshift)]

def test():
    yield 1
    yield 2
    yield 3  
    
    