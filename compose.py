#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 15 16:50:47 2023

@author: kirbyurner
"""

def f(x):
    return 2 * x

def g(x):
    return 2 + x

#%%

print("f(g(10))=", f(g(10)))
print("g(f(10))=", g(f(10)))

#%%

print([f(g(k)) for k in range(10)])
print([g(f(k)) for k in range(10)])

#%%

def compose(m, n): 
    """compose callables"""
    def newfunc(x):
        return m(n(x))
    return newfunc

h = compose(f, g)
print("h type:", type(h))
print("f(g(10)): ", h(10))

h = compose(g, f)
print("g(f(10)): ", h(10))

#%%

class Compose:
    
    def __init__(self, func):
        """store callable"""
        self.func = func
    
    def __mul__(self, other):
        """compose callables"""
        return Compose(lambda x: self(other(x)))
        
    def __call__(self, x):
        """call self with argument"""
        return self.func(x)
    
f = Compose(f)
g = Compose(g)

h = f * g
print("h type:", type(h))
print("f(g(10))= ", h(10))

h = g * f
print("g(f(10))= ", h(10))

#%%

@Compose
def f(x):
    return 2 * x

@Compose
def g(x):
    return 2 + x

h = f * g
print("h type:", type(h))
print("f(g(10))= ", h(10))

h = g * f
print("g(f(10))= ", h(10))
