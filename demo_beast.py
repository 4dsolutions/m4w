#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:18:26 2024

@author: Kirby Urner
"""

from sympy import N, Rational, sqrt, Eq, Symbol
from tetravolume import B,E,A,S,T

φ = (1 + sqrt(5))/2
𝛕 = 1/φ 
E_vol = sqrt(2) * Rational(1,8) * 𝜏 ** 3
S_vol = Rational(1,2) * 𝜏 ** 5
T_vol = A_vol = B_vol = Rational(1,24)
sfactor = (S_vol/E_vol).simplify()
    
def demo():
    print("Algebraic using φ...")
    print("E_vol = ", E_vol) # keep algebraic
    print("S_vol = ", S_vol) # keep algebraic
    print()
    print("Numeric evaluation...")
    print("E_vol = ", N(E_vol, 50))  # evaluate
    print("S_vol = ", N(S_vol, 50))  # evaluate
    print()
    print("Algebraic using CM Volume Formula...")
    b,e,a,s,t = B(),E(),A(),S(),T()

    print("B volume =", b.ivm_volume()) # keep algebraic
    print("E volume =", e.ivm_volume()) # keep algebraic
    print("A volume =", a.ivm_volume()) # keep algebraic
    print("S volume =", s.ivm_volume()) # keep algebraic
    print("T volume =", t.ivm_volume()) # keep algebraic
    print()
    
    print("Numeric evaluation...")
    print("E volume =", e.ivm_volume(True, 50))  # evaluate
    print("S volume =", s.ivm_volume(True, 50))  # evaluate
    print()
    print("Testing algebraic equivalence...")
    print("Two E expressions equivalent?: ", Eq(e.ivm_volume(), E_vol))
    print("Two S expressions equivalent?: ", Eq(s.ivm_volume(), S_vol))
    print()
    print("SkewIcosa + 24 S modules...")
    print("sfactor = ", sfactor.simplify(), " = ", N(sfactor.simplify()))

    # IcosaWithin + 24 S mods
    icosa = (Rational(5,2) * sfactor * sfactor).simplify()
    octa = (icosa + 24 * S_vol).simplify()
    
    print("Icosa   = ", icosa, " = ", N(icosa))
    print("Octa 4  = ", octa, " = ", N(octa))
    print("RT5     = ", (120 * T_vol), " = ", N(120 * T_vol))
    print("RT5+    = ", (120 * E_vol), " = ", N(120 * E_vol))
    print("Icosa   = ", 20 * (1/sfactor), " = ", N(20 * (1/sfactor)))
    print("Cubocta = ", 20)
    print("SuperRT = ", (120 * E_vol) * φ**3, " = ", N((120 * E_vol) * φ**3))

if __name__ == "__main__":
    demo()