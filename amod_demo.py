#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 04:12:48 2024

@author: kirbyurner
"""

from tetravolume import A, CM, qvolume

amod = A()
# print(amod.edges())
# print(amod.angles())

print("QRAY:", qvolume(amod.amod_E, amod.amod_D, 
              amod.amod_C, amod.amod_F))

print("CM:  ", CM(amod.a, amod.b, amod.c, 
                  amod.d, amod.e, amod.f))