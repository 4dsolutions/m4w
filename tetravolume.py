"""
tetravolume.py
Kirby Urner (c) MIT License

April  1, 2024: wire up all three volume-from-edges algorithms as options 

March 30, 2024: tighten the unittests to rely more on sympy
March 30, 2024: save edges and angles when initializing the Tetrahedron
March 30, 2024: make_tet now returns Tetrahedon (not volumes tuple)

The tetravolume.py methods make_tet and make_tri 
assume that volume and area use R-edge cubes and 
triangles for XYZ units respectively, and D-edge 
tetrahedrons and triangles for IVM units of volume 
and area (D = 2R).

The tetrahedron of edges D has sqrt(8/9) the 
volume of a cube of edges R, yet each is unit in
its respective matrix.

The triangle of edges D has an XYZ area of sqrt(3)
i.e. an equilateral triangle of edges 2 in R-square
units.  The IVM area of the same triangle is simply 1.

The cube of edges sqrt(2) in R units, has volume 
sqrt(2) to the 3rd power.  One third of that volume
is our unit tetrahedron of edges D (cube face diagonals).

See:
https://coda.io/d/Math4Wisdom_d0SvdI3KSto/ivm-xyz_suqdu#_luR7B
for explanation of quadrays, used for some unit tests

https://flic.kr/p/2pGmvWD (labeling system)

https://github.com/4dsolutions/School_of_Tomorrow/blob/master/VolumeTalk.ipynb
for background on adapting volume formulas and/or using GdJ's

A goal for this version of tetravolume.py + qrays.py
is to keep computations symbolic, with precision 
open-ended. A work in progress.  Mar 5, 2024.

Another goal as of March-April 2024 is to flesh out
the Tetrahedron instance with more specific info as
to which segments are what length.

Apex A goes to base B, C, D, creating edges:
    a: AB
    b: AC
    c: AD
    d: BC
    e: CD
    f: BD
(alphabetized pairs)
    
And faces: 
    ABC ACD ADB BCD (in terms of verts)
or (in terms of lengths): 
    (a, d, b), (b, e, c), (c, f, a), (d, e, f)
    [ tuples any order, elements may cycle ]
    
Therefore we have three angles from each vertex:
A: BAC CAD BAD
B: ABC CBD ABD
C: ACB ACD BCD
D: ADC BDC ADB
(middle letter is vertex angle, left and right letters alphabetized)

When we get the six lengths as inputs, lets 
assign them to edge names and compute the 
twelve angles.

Example resource:
https://www.omnicalculator.com/math/triangle-angle
"""

import sympy as sp
from sympy import Rational, Integer, Matrix, acos, deg, N
from sympy import sqrt as rt2
from qrays import Qvector, Vector
import sys

R = Rational(1,2)
D = Integer(1)

Syn3  = rt2(sp.Rational(9, 8))
root2 = rt2(2)
root3 = rt2(3)
root5 = rt2(5)
root6 = rt2(6)

PHI = (1 + root5)/2

Smod = (PHI **-5)/2  
Emod = (root2/8) * (PHI ** -3)
Amod = Bmod = Tmod = Rational(1,24)

sfactor = Smod/Emod

class Tetrahedron:
    """
    Takes six edges of tetrahedron with faces
    (a,b,d)(b,c,e)(c,a,f)(d,e,f) -- returns volume
    in ivm tvs and xyz cubic units, S3 ratio.
    
        Apex A goes to base B, C, D, creating edges:
            a: AB
            b: AC
            c: AD
            d: BC
            e: CD
            f: DB  
            
    https://flic.kr/p/2pGmvWD (labeling system)
    """

    def __init__(self, a, b, c, d, e, f):
        self.a, self.a2 = a, a**2
        self.b, self.b2 = b, b**2
        self.c, self.c2 = c, c**2
        self.d, self.d2 = d, d**2
        self.e, self.e2 = e, e**2
        self.f, self.f2 = f, f**2
   
        # 2-letter edge labels
        self.AB = a
        self.AC = b
        self.AD = c
        self.BC = d
        self.CD = e
        self.DB = f
        
        # 3-letter face angles
        a2,b2,c2,d2,e2,f2 = self.a2, self.b2, self.c2, self.d2, self.e2, self.f2

        self.BAC = acos( (a2 + b2 - d2)/(2 * a * b) )
        self.CAD = acos( (b2 + c2 - e2)/(2 * b * c) )
        self.BAD = acos( (a2 + c2 - f2)/(2 * a * c) )
        
        self.ABC = acos( (a2 + d2 - b2)/(2 * a * d) )
        self.CBD = acos( (d2 + f2 - e2)/(2 * d * f) )
        self.ABD = acos( (a2 + f2 - c2)/(2 * a * f) )
        
        self.ACB = acos( (b2 + d2 - a2)/(2 * b * d) )
        self.ACD = acos( (b2 + e2 - c2)/(2 * b * e) ) 
        self.BCD = acos( (d2 + e2 - f2)/(2 * d * e) )
        
        self.ADC = acos( (c2 + e2 - b2)/(2 * c * e) ) 
        self.BDC = acos( (e2 + f2 - d2)/(2 * e * f) ) 
        self.ADB = acos( (c2 + f2 - a2)/(2 * c * f) )

    def dump(self):
        return self.a2, self.b2, self.c2
    
    def __mul__(self, other):
        a = self.a * other
        b = self.b * other
        c = self.c * other
        d = self.d * other
        e = self.e * other
        f = self.f * other
        return Tetrahedron(a,b,c,d,e,f)
        
    __rmul__ = __mul__
        
    def edges(self):
        """
            a: AB
            b: AC
            c: AD
            d: BC
            e: CD
            f: BD 
        """
        return {
            "AB": self.a,
            "AC": self.b,
            "AD": self.c,
            "BC": self.d,
            "DC": self.e,
            "BD": self.f,
            }
        
    def angles(self, values=False, prec=15):
        """
        Three angles from each vertex:
        A: BAC CAD BAD
        B: ABC CBD ABD
        C: ACB ACD BCD
        D: ADC BDC ADB
        (middle letter is vertex angle, left and right letters alphabetized)
        """
        if values:
            return {
                "BAC": N(self.BAC, prec),
                "CAD": N(self.CAD, prec),
                "BAD": N(self.BAD, prec),
                "ABC": N(self.ABC, prec),
                "CBD": N(self.CBD, prec),
                "ABD": N(self.ABD, prec),
                "ACB": N(self.ACB, prec),
                "ACD": N(self.ACD, prec),
                "BCD": N(self.BCD, prec),
                "ADC": N(self.ADC, prec),
                "BDC": N(self.BDC, prec),
                "ADB": N(self.ADB, prec)
                }            
        else:
            return {
                "BAC": self.BAC,
                "CAD": self.CAD,
                "BAD": self.BAD,
                "ABC": self.ABC,
                "CBD": self.CBD,
                "ABD": self.ABD,
                "ACB": self.ACB,
                "ACD": self.ACD,
                "BCD": self.BCD,
                "ADC": self.ADC,
                "BDC": self.BDC,
                "ADB": self.ADB
                }

    def degrees(self, values=False, prec=15):
        output = {}
        if values:
            for k,v in self.angles().items():
                output[k] = N(deg(v), prec) 
        else:
            for k,v in self.angles().items():
                output[k] = deg(v)
        return output
            
    def ivm_volume(self, value=False, prec=50):
        """
        Three options to compute volume from edges...
        GdJ: Gerald de Jong, similar to Euler's, lost his notes, works
        PdF: Pierro della Francesca, modified by Syn3 (XYZ->IVM constant)
        CM : Caley-Menger, modified by Syn3 (XYZ->IVM constant) 
        """
        
        #ivmvol = GdJ(self.a, self.b, self.c, self.d, self.e, self.f)
        #ivmvol = PdF(self.a, self.b, self.c, self.d, self.e, self.f)
        ivmvol = CM(self.a, self.b, self.c, self.d, self.e, self.f)
        
        return ivmvol if not value else N(ivmvol, prec)

    def xyz_volume(self, value=False, prec=50):
        xyzvol = (1/Syn3) * self.ivm_volume()
        return xyzvol if not value else N(xyzvol, prec)

def GdJ(a, b, c, d, e, f):
    "Gerald de Jong"
    A,B,C,D,E,F = [x**2 for x in (a, b, c, d, e, f)] # 2nd power us

    _open   = sum((A * B * E, A * B * F, A * C * D,  
                   A * C * E, A * D * E, A * E * F,
                   B * C * D, B * C * F, B * D * F, 
                   B * E * F, C * D * E, C * D * F))
    
    _closed = sum((A * B * D, 
                   A * C * F, 
                   B * C * E, 
                   D * E * F))

    _oppo   = sum((A * E * (A + E),
                   B * F * (B + F),
                   C * D * (C + D)))
    
    return rt2((_open - _closed - _oppo)/2)

def PdF(a,b,c,d,e,f):
    """
    Pierro della Francesca
    """
    
    def adapter(a, e, c, d, f, b):
        "unscramble input's to match GdJ order"
        return a, b, c, d, e, f

    A,B,C,D,E,F = [x**2 for x in 
                   adapter(2*a,2*b,2*c,2*d,2*e,2*f)] 
    
    comp_chunk =  ((A * F) * (-A + B + C + D + E - F)
                 + (B * E) * ( A - B + C + D - E + F)
                 + (C * D) * ( A + B - C - D + E + F)
                 - (A + F) * (B + E) * (C + D)/2
                 - (A - F) * (B - E) * (C - D)/2 )
    
    return rt2(2 * comp_chunk) / 16  # Syn3 is blended in here

def CM(a, b, c, d, e, f):
    """
    Caley-Menger
    """
    A,B,C,D,E,F = [(2*x)**2 for x in (a,b,c,d,e,f)]
    
    # Construct a 5x5 matrix per Caley-Menger
    M = Matrix(
        [[0, 1, 1, 1, 1],
         [1, 0, A, B, C],
         [1, A, 0, D, F],
         [1, B, D, 0, E],
         [1, C, F, E, 0]])
    return rt2(M.det())/16  # Syn3 factored in 

def make_tet(v0,v1,v2):
    """
    three edges from any corner, remaining three edges computed
    """
    return Tetrahedron(v0.length(), v1.length(), v2.length(), 
                      (v0-v1).length(), (v1-v2).length(), (v2-v0).length())

class Triangle:
    
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c  

    def ivm_area(self):
        ivmarea = self.xyz_area() * 1/rt2(3)
        return ivmarea

    def xyz_area(self):
        """
        Heron's Formula, without the 1/4
        """
        a,b,c = self.a, self.b, self.c
        xyzarea = rt2((a+b+c) * (-a+b+c) * (a-b+c) * (a+b-c))
        return xyzarea
    
def make_tri(v0,v1):
    """
    three edges from any corner, remaining three edges computed
    """
    tri = Triangle(v0.length(), v1.length(), (v1-v0).length())
    return tri.ivm_area(), tri.xyz_area()


class E(Tetrahedron):
    
    def __init__(self):
        e0 = D/2
        e1 = root3 * PHI**-1 /2
        e2 = rt2((5 - root5)/2)/2
        e3 = (3 - root5)/2/2
        e4 = rt2(5 - 2*root5)/2
        e5 = 1/PHI/2        
        super().__init__(e0, e1, e2, e3, e4, e5)

class T(Tetrahedron):
    
    def __init__(self):
        E2T = (2**sp.Rational(5,6) * 3**sp.Rational(2,3) * PHI / 6).simplify()
        e0 = D/2
        e1 = root3 * PHI**-1 /2
        e2 = rt2((5 - root5)/2)/2
        e3 = (3 - root5)/2/2
        e4 = rt2(5 - 2*root5)/2
        e5 = 1/PHI/2  
        
        e0 *= E2T
        e1 *= E2T
        e2 *= E2T
        e3 *= E2T
        e4 *= E2T
        e5 *= E2T
 
        super().__init__(e0, e1, e2, e3, e4, e5)
        
class S(Tetrahedron):
    
    def __init__(self):
        e0 = 1/PHI
        e1 = sfactor/2
        e2 = root3 * e1/2
        e3 = (3 - root5)/2
        e4 = e1/2
        e5 = e4
        super().__init__(e0, e1, e2, e3, e4, e5)    
        
import unittest
class Test_Tetrahedron(unittest.TestCase):

    def test_unit_volume(self):
        tet = Tetrahedron(D, D, D, D, D, D)
        self.assertEqual(tet.ivm_volume(), Integer(1), "Volume not 1")

    def test_e_module(self):
        tet = E()
        self.assertTrue(sp.Eq(tet.ivm_volume(), 
                              (root2/8) * (PHI ** -3)))
        
    def test_s_module(self):
        tet = S()
        self.assertTrue(sp.Eq(tet.ivm_volume(), 
                               (PHI ** -5)/2))
        
    def test_unit_volume2(self):
        tet = Tetrahedron(R, R, R, R, R, R)
        self.assertTrue(sp.Eq(tet.xyz_volume(), root2/12))

    def test_unit_volume3(self):
        tet = Tetrahedron(R, R, R, R, R, R)
        self.assertTrue(sp.Eq(tet.ivm_volume(), Rational(1,8)))
        
    def test_phi_edge_tetra(self):
        tet = Tetrahedron(D, D, D, D, D, PHI)
        self.assertTrue(sp.Eq(tet.ivm_volume(), root2/2))

    def test_right_tetra(self):
        e = rt2((root3/2)**2 + (root3/2)**2)  # right tetrahedron
        tet = Tetrahedron(D, D, D, D, D, e)
        self.assertEqual(tet.xyz_volume(), D)

    def test_quadrant(self):
        one = Integer(1)
        qA = Qvector((one,0,0,0))
        qB = Qvector((0,one,0,0))
        qC = Qvector((0,0,one,0))
        tet = make_tet(qA, qB, qC) 
        self.assertEqual(tet.ivm_volume(), Rational(1,4)) 

    def test_octant(self):
        x = Vector((R, 0, 0))
        y = Vector((0, R, 0))
        z = Vector((0, 0, R))
        tet = make_tet(x,y,z)
        self.assertEqual(tet.xyz_volume(), Rational(1,6))

    def test_quarter_octahedron(self):
        a = Vector((D,0,0))
        b = Vector((0,D,0))
        c = Vector((R,R,root2/2))
        tet = make_tet(a, b, c)
        self.assertEqual(tet.ivm_volume(), D) 

    def test_xyz_cube(self):
        a = Vector((R,0,0))
        b = Vector((0,R,0))
        c = Vector((0,0,R))
        R_octa = make_tet(a,b,c) 
        self.assertEqual(6 * R_octa.xyz_volume(), D) # good to 4 places  

    def test_s3(self):
        D_tet = Tetrahedron(D, D, D, D, D, D)
        a = Vector((R,0,0))
        b = Vector((0,R,0))
        c = Vector((0,0,R))
        R_cube = 6 * make_tet(a,b,c).xyz_volume()
        self.assertEqual(D_tet.xyz_volume() * Syn3, R_cube)

    def test_martian(self):
        two = Integer(2)
        one = Integer(1)
        p = Qvector((two,one,0,one))
        q = Qvector((two,one,one,0))
        r = Qvector((two,0,one,one))
        result = make_tet(5*q, 2*p, 2*r)
        self.assertEqual(result.ivm_volume(), Integer(20))
        
    def test_area_martian1(self):
        two = Integer(2)
        one = Integer(1)
        p = Qvector((two,one,0,one))
        q = Qvector((two,one,one,0))
        result = p.area(q)
        self.assertEqual(result, D)        
 
    def test_area_martian2(self):
        two = Integer(2)
        one = Integer(1)
        p = 3 * Qvector((two,one,0,one))
        q = 4 * Qvector((two,one,one,0))
        result = p.area(q)
        self.assertEqual(result, 12)

    def test_area_martian3(self):
        qx = Vector((D,0,0)).quadray()
        qy = Vector((R,rt2(3)/2,0)).quadray()
        result = qx.area(qy)
        self.assertEqual(result, D)
        
    def test_area_earthling1(self):
        vx = Vector((1,0,0))
        vy = Vector((0,1,0))
        result = vx.area(vy)
        self.assertEqual(result, 1)        

    def test_area_earthling2(self):
        vx = Vector((2,0,0))
        vy = Vector((1,rt2(3),0))
        result = vx.area(vy)
        self.assertEqual(result, 2*rt2(3))    
        
    def test_phi_tet(self):
        "edges from common vertex: phi, 1/phi, 1"
        p = Vector((1, 0, 0))
        q = Vector((1, 0, 0)).rotz(60) * PHI # rotation is still a bit fuzzy
        r = Vector((Rational(1,2), root3/6, root6/3)) * 1/PHI
        result = make_tet(p, q, r)
        self.assertAlmostEqual(N(result.ivm_volume()), D)
        
    def test_phi_tet_2(self):
        two = Integer(2)
        one = Integer(1)
        p = Qvector((two,one,0,one))
        q = Qvector((two,one,one,0))
        r = Qvector((two,0,one,one))
        result = make_tet(PHI*q, (1/PHI)*p, r)
        self.assertTrue(sp.Eq(result.ivm_volume(), D))
        
    def test_phi_tet_3(self):
        T = Tetrahedron(PHI, 1/PHI, 1, 
                        root2, root2/PHI, root2)
        result = T.ivm_volume()
        self.assertTrue(sp.Eq(result, D))

    def test_koski(self):
        a = 1 
        b = PHI ** -1
        c = PHI ** -2
        d = (root2) * PHI ** -1 
        e = (root2) * PHI ** -2
        f = (root2) * PHI ** -1       
        T = Tetrahedron(a,b,c,d,e,f)
        result = T.ivm_volume()
        self.assertTrue(sp.Eq(result, PHI ** -3))      


class Test_Koski(unittest.TestCase):
        
    def test_Tetrahedron(self):
        "Tetrahedron =   S6  +    S3 # (volume = 1)"
        S6 = S() * PHI**2
        S3 = S() * PHI
        self.assertTrue(sp.Eq(S6.ivm_volume() + S3.ivm_volume(), D))
    
    def test_Icosahedron(self):
        "Icosahedron =100*E3 + 20*E  # (volume = ~18.51)"
        E3 = E() * PHI
        E0 = E()
        self.assertEqual(N(100*(E3.ivm_volume()) + 20*(E0.ivm_volume())), 
                         N(20 * 1/sfactor) )       

    def test_RT5(self):
        "RhTriac_T   =  5*S6 +  5*S3 # (volume = 5)"
        S6 = S() * PHI**2
        S3 = S() * PHI
        self.assertTrue(sp.Eq(5*(S6.ivm_volume()) + 5*(S3.ivm_volume()), 5))
        
        
class Test_Triangle(unittest.TestCase):
    
    def test_unit_area1(self):
        tri = Triangle(D, D, D)
        self.assertEqual(tri.ivm_area(), 1)
        
    def test_unit_area2(self):
        tri = Triangle(2, 2, 2)
        self.assertEqual(tri.ivm_area(), 4)
        
    def test_xyz_area3(self):
        tri = Triangle(D, D, D)
        self.assertEqual(tri.xyz_area(), rt2(3))
        
    def test_xyz_area4(self):
        v1 = Vector((D, 0, 0))
        v2 = Vector((0, D, 0))
        xyz_area = make_tri(v1, v2)[1]
        self.assertTrue(sp.Eq(xyz_area, 2))

    def test_xyz_area5(self):
        tri = Triangle(R, R, R)
        self.assertEqual(tri.xyz_area(), (root3)/4)
        
def command_line():
    args = sys.argv[1:]
    try:
        args = [float(x) for x in args] # floats
        t = Tetrahedron(*args)
    except TypeError:
        t = Tetrahedron(1,1,1,1,1,1)
        print("defaults used")
    print(t.ivm_volume())
    print(t.xyz_volume())
        
if __name__ == "__main__":
    if len(sys.argv)==7:
        command_line()  
    else:
        unittest.main()
