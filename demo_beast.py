#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 16:18:26 2024

@author: Kirby Urner

The doodle that got the ball rolling:
https://www.flickr.com/photos/kirbyurner/53631820894

"""

from sympy import N, Rational, sqrt, Eq
from tetravolume import B,E,A,S,T

φ = (1 + sqrt(5))/2
𝛕 = 1/φ 
b,e,a,s,t = B(), E(), A(), S(), T()
E_vol = sqrt(2) * Rational(1,8) * 𝜏 ** 3
S_vol = Rational(1,2) * 𝜏 ** 5
T_vol = A_vol = B_vol = Rational(1,24)
sfactor = (S_vol/E_vol).simplify()
Syn3 = sqrt(Rational(9,8))
    
def ivm_vols(n=50):
    print("Algebraic using φ...")
    print("E_vol = ", E_vol) # keep algebraic
    print("S_vol = ", S_vol) # keep algebraic
    print()
    print("Numeric evaluation...")
    print("E_vol = ", N(E_vol, n))  # evaluate
    print("S_vol = ", N(S_vol, n))  # evaluate
    print()
    print("Algebraic using Volume Formula...")
    print("B volume =", b.ivm_volume())
    print("E volume =", e.ivm_volume().simplify())
    print("A volume =", a.ivm_volume())
    print("S volume =", s.ivm_volume().simplify())
    print("T volume =", t.ivm_volume().simplify())
    print()
    
    print("Numeric evaluation...")
    print("B volume =", b.ivm_volume(True, n))
    print("E volume =", e.ivm_volume(True, n))
    print("A volume =", a.ivm_volume(True, n))
    print("S volume =", s.ivm_volume(True, n))
    print("T volume =", t.ivm_volume(True, n)) 
    print()
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

def xyz_vols(n=50):
    print("Algebraic using φ...")
    print("E_vol = ", E_vol * 1/Syn3) # keep algebraic
    print("S_vol = ", S_vol * 1/Syn3) # keep algebraic
    print()
    print("Numeric evaluation...")
    print("E_vol = ", N(E_vol * 1/Syn3, n))  # evaluate
    print("S_vol = ", N(S_vol * 1/Syn3, n))  # evaluate
    print()
    print("Algebraic using Volume Formula...")
    print("B volume =", b.xyz_volume().simplify())
    print("E volume =", e.xyz_volume().simplify())
    print("A volume =", a.xyz_volume().simplify())
    print("S volume =", s.xyz_volume().simplify())
    print("T volume =", t.xyz_volume().simplify())
    print()
    
    print("Numeric evaluation...")
    print("B volume =", b.xyz_volume(True, n))
    print("E volume =", e.xyz_volume(True, n))
    print("A volume =", a.xyz_volume(True, n))
    print("S volume =", s.xyz_volume(True, n))
    print("T volume =", t.xyz_volume(True, n)) 
    print()
    print()
    print("Testing algebraic equivalence...")
    print("Two E expressions equivalent?: ", Eq(e.ivm_volume(), E_vol))
    print("Two S expressions equivalent?: ", Eq(s.ivm_volume(), S_vol))
    
def conversion_constants():
    print("B2E = ", (e.ivm_volume()/b.ivm_volume()).simplify() )
    print("E2B = ", (b.ivm_volume()/e.ivm_volume()).simplify() )
    
    print("E2A = ", (a.ivm_volume()/e.ivm_volume()).simplify() )
    print("A2E = ", (e.ivm_volume()/a.ivm_volume()).simplify() )
    
    print("A2S = ", (s.ivm_volume()/a.ivm_volume()).simplify() )
    print("S2A = ", (a.ivm_volume()/s.ivm_volume()).simplify() )
    
    print("S2T = ", (t.ivm_volume()/s.ivm_volume()).simplify() )
    print("T2S = ", (s.ivm_volume()/t.ivm_volume()).simplify() )
    
    print("T2B = ", (b.ivm_volume()/t.ivm_volume()).simplify() )
    print("B2T = ", (t.ivm_volume()/b.ivm_volume()).simplify() )
    
    print("B2S = ", (s.ivm_volume()/b.ivm_volume()).simplify() )
    print("S2B = ", (b.ivm_volume()/s.ivm_volume()).simplify() )
    
    print("B2A = ", (a.ivm_volume()/b.ivm_volume()).simplify() )
    print("A2B = ", (b.ivm_volume()/a.ivm_volume()).simplify() )
    
    print("E2T = ", (t.ivm_volume()/e.ivm_volume()).simplify() )
    print("T2E = ", (e.ivm_volume()/t.ivm_volume()).simplify() )

    print("E2S = ", (s.ivm_volume()/e.ivm_volume()).simplify() )
    print("S2E = ", (e.ivm_volume()/s.ivm_volume()).simplify() )
    
    print("A2T = ", (t.ivm_volume()/a.ivm_volume()).simplify() )
    print("T2A = ", (a.ivm_volume()/t.ivm_volume()).simplify() ) 
    
def conversion_values(n=50):
    
    print("B2E = ", (e.ivm_volume()/b.ivm_volume()).evalf(n) )
    print("E2B = ", (b.ivm_volume()/e.ivm_volume()).evalf(n))
    
    print("E2A = ", (a.ivm_volume()/e.ivm_volume()).evalf(n) )
    print("A2E = ", (e.ivm_volume()/a.ivm_volume()).evalf(n) )
    
    print("A2S = ", (s.ivm_volume()/a.ivm_volume()).evalf(n) )
    print("S2A = ", (a.ivm_volume()/s.ivm_volume()).evalf(n) )
    
    print("S2T = ", (t.ivm_volume()/s.ivm_volume()).evalf(n) )
    print("T2S = ", (s.ivm_volume()/t.ivm_volume()).evalf(n) )
    
    print("T2B = ", (b.ivm_volume()/t.ivm_volume()).evalf(n) )
    print("B2T = ", (t.ivm_volume()/b.ivm_volume()).evalf(n) )
    
    print("B2S = ", (s.ivm_volume()/b.ivm_volume()).evalf(n) )
    print("S2B = ", (b.ivm_volume()/s.ivm_volume()).evalf(n) )
    
    print("B2A = ", (a.ivm_volume()/b.ivm_volume()).evalf(n) )
    print("A2B = ", (b.ivm_volume()/a.ivm_volume()).evalf(n) )
    
    print("E2T = ", (t.ivm_volume()/e.ivm_volume()).evalf(n) )
    print("T2E = ", (e.ivm_volume()/t.ivm_volume()).evalf(n) )

    print("E2S = ", (s.ivm_volume()/e.ivm_volume()).evalf(n) )
    print("S2E = ", (e.ivm_volume()/s.ivm_volume()).evalf(n) )
    
    print("A2T = ", (t.ivm_volume()/a.ivm_volume()).evalf(n) )
    print("T2A = ", (a.ivm_volume()/t.ivm_volume()).evalf(n) ) 
    
if __name__ == "__main__":
    # ivm_vols(n=20)
    xyz_vols(n=20)
    # conversion_constants()
    # conversion_values()