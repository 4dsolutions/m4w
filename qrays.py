# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 09:07:22 2016

Notes:

Jul 26, 2024: This version was about stripping out
anything that might keep the computation objects 
algebraic (symbolic) versus allowing them to decay
into numerics. Some cool Qvector methods were
sacrificed towards this end, but not in all versions.

In some versions of Qvector, we have some elegant
4x4 rotation matrix coupled with norm0 operations
per the Tom Ace website.

https://minortriad.com/quadray.html

My current workaround is a hack in that I allow 
the Qvector parent class, Vector, to handle rotation 
using conventional XYZ methods.  Not optimum.

For more background on Quadrays in the M4W context,
see: 

https://coda.io/d/Math4Wisdom_d0SvdI3KSto/Synergetics_su5DS#_luecx

Vectors and Qvectors use the same metric i.e. the 
xyz vector and corresponding ivm vector always have 
the same length.

In contrast, the tetravolume.py modules in some cases 
assumes that volume and area use R-edge cubes and triangles
for XYZ units respectively, and D-edge tetrahedrons 
and triangles for IVM units of volume and area.  See
the docstring for more details.

@author:  K. Urner, 4D Solutions, (M) MIT License
 Mar 10, 2024: the m4w.qrays.py fork aims to stay sympy friendly
 Mar  5, 2024: continuing to remove any fixed-precision limitations
 Nov 15, 2023: customized for use with sympy, mpmath in m4w repo (incomplete)
 Oct  8, 2021: remove gmpy2 dependency
 Sep 19, 2021: remove autoconvert to floating point when initializing Vector
 Sep 19, 2021: make xyz Vector a property of Qvector vs. a method
 Sep 06, 2019: have type(self)() instead of Qvector() return outcomes
 May 25, 2019: add area methods based on cross product
 Jun 20, 2018: make Qvectors sortable, hashable
 Jun 11, 2016: refactored to make Qvector and Vector each standalone
 Aug 29, 2000: added extra-class, class dependent methods for
    dot and cross as alternative syntax
 June 6, 2020: spherical coordinates debug, working on blender integration
 July 8,2000: added method for rotation around any axis vector
 May 27,2000: shortend several methods thanks to Peter Schneider-Kamp
 May 24,2000: added unit method, tweaks
 May 8, 2000: slight tweaks re rounding values
 May 7, 2000: enhanced the Qvector subclass with native
    length, dot, cross methods -- keeps coords as a 4-tuple
    -- generalized Vector methods to accommodate 4-tuples
    if Qvector subclass, plus now returns vector of whatever
    type invokes method (i.e. Qvector + Qvector = Qvector)
 Mar 23, 2000:
 added spherical coordinate subclass
 added quadray coordinate subclass
 Mar  5, 2000: added angle function
"""

from sympy import cos, sin, acos, atan, N, Rational, Integer, Sum, S, rad, deg
from sympy import sqrt as rt2
import mpmath
from operator import add, sub, mul, neg
from collections import namedtuple
import sympy as sp

XYZ = namedtuple("xyz_vector", "x y z")
IVM = namedtuple("ivm_vector", "a b c d")

zero       = Integer(0)
one        = Integer(1)
root2      = rt2(2)
root3      = rt2(3)
half       = S.Half
mpmath.dps = 50

class Vector:

    def __init__(self, arg):
        """Initialize a vector at an (x,y,z)"""
        self.xyz = XYZ(*arg)

    def __repr__(self):
        return repr(self.xyz)
    
    @property
    def x(self):
        return self.xyz.x

    @property
    def y(self):
        return self.xyz.y

    @property
    def z(self):
        return self.xyz.z
        
    def __mul__(self, scalar):
        """Return vector (self) * scalar."""
        newcoords = [scalar * dim for dim in self.xyz]
        return type(self)(newcoords)

    __rmul__ = __mul__ # allow scalar * vector

    def __truediv__(self,scalar):
        """Return vector (self) * 1/scalar"""        
        return self.__mul__(1/scalar)
    
    def __add__(self,v1):
        """Add a vector to this vector, return a vector""" 
        newcoords = map(add, v1.xyz, self.xyz)
        return type(self)(newcoords)
        
    def __sub__(self,v1):
        """Subtract vector from this vector, return a vector"""
        return self.__add__(-v1)
    
    def __neg__(self):      
        """Return a vector, the negative of this one."""
        return type(self)(tuple(map(neg, self.xyz)))

    def unit(self):
        return self.__mul__(1/self.length())

    def dot(self, v1):
        """
        Return scalar dot product of this with another vector.
        """
        return (self.x * v1.x) + (self.y * v1.y) + (self.z * v1.z)

    def cross(self,v1):
        """
        Return the vector cross product of this with another vector
        """
        newcoords = (self.y * v1.z - self.z * v1.y, 
                     self.z * v1.x - self.x * v1.z,
                     self.x * v1.y - self.y * v1.x )
        return type(self)(newcoords)
    
    def area(self,v1):
        """
        xyz area of a parallelogram with these edge lengths
        """
        return self.cross(v1).length()
    
    def length(self):
        """Return this vector's length"""
        return rt2(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def angle(self, v1):
        """
        return in degrees
        """
        costheta = self.dot(v1)/(self.length() * v1.length())
        return deg(acos(costheta))

    def rotaxis(self,vAxis,degr):
        """
        Rotate around vAxis by deg
        realign rotation axis with Z-axis, realign self accordingly,
        rotate by deg (counterclockwise) around Z, resume original
        orientation (undo realignment)
        """
        
        r,phi,theta = vAxis.spherical()
        newv  = self.rotz(-theta).roty(phi)
        newv  = newv.rotz(-degr)
        newv  = newv.roty(-phi).rotz(theta)
        return type(self)(newv.xyz)        

    def rotx(self, degr):
        rads   = rad(degr)
        newy   = cos(rads) * self.y - sin(rads) * self.z
        newz   = sin(rads) * self.y + cos(rads) * self.z
        newxyz = (self.x, newy, newz)
        return type(self)(newxyz)
   
    def roty(self, degr):
        rads   = rad(degr)
        newx   = cos(rads) * self.x - sin(rads) * self.z
        newz   = sin(rads) * self.x + cos(rads) * self.z
        newxyz = (newx, self.y, newz)
        return type(self)(newxyz)

    def rotz(self, degr):
        rads   = rad(degr)
        newx   = cos(rads) * self.x - sin(rads) * self.y
        newy   = sin(rads) * self.x + cos(rads) * self.y
        newxyz = (newx, newy, self.z)
        return type(self)(newxyz)
    
    def spherical(self):
        """Return (r,phi,theta) spherical coords based 
        on current (x,y,z)"""
        r = self.length()
        
        if self.x == 0:
            if   self.y ==0: theta =   0.0
            elif self.y < 0: theta = -90.0
            else:            theta =  90.0
            
        else:  
            
            theta = deg(atan(self.y/self.x))
            if   self.x < 0 and self.y == 0:   theta = 180
            # theta is positive so turn more than 180
            elif self.x < 0 and self.y <  0:   theta = 180 + theta
            # theta is negative so turn less than 180
            elif self.x < 0 and self.y >  0:   theta = 180 + theta

        if r == 0: 
            phi=0.0
        else: 
            phi = deg(acos(self.z/r))
        
        return (r, phi, theta)

    def quadray(self):
        """return (a, b, c, d) quadray based on current (x, y, z)"""
        x, y, z = self.xyz
        k = 2/root2
        a = k * (int(bool(x >= 0)) * ( x) + int(bool(y >= 0)) * ( y) + int(bool(z >= 0)) * ( z))
        b = k * (int(bool(x <  0)) * (-x) + int(bool(y <  0)) * (-y) + int(bool(z >= 0)) * ( z))
        c = k * (int(bool(x <  0)) * (-x) + int(bool(y >= 0)) * ( y) + int(bool(z <  0)) * (-z))
        d = k * (int(bool(x >= 0)) * ( x) + int(bool(y <  0)) * (-y) + int(bool(z <  0)) * (-z))
        return Qvector((a, b, c, d))

        
class Qvector(Vector):
    """Quadray vector"""

    def __init__(self, arg):
        """Initialize a vector at an (a, b, c, d)"""
        self.coords = self.norm(arg)

    def __repr__(self):
        return repr(self.coords)

    def norm(self, arg):
        """Normalize such that 4-tuple all non-negative members."""
        minarg = min(arg)
        return IVM(arg[0] - minarg, arg[1] - minarg, arg[2] - minarg, arg[3] - minarg) 
    
    def norm0(self):
        """Normalize such that sum of 4-tuple members = 0"""
        q  = self.coords
        av = (q[0] + q[1] + q[2] + q[3])/4
        return IVM(q[0]-av, q[1]-av, q[2]-av, q[3]-av) 

    @property
    def a(self):
        return self.coords.a

    @property
    def b(self):
        return self.coords.b

    @property
    def c(self):
        return self.coords.c

    @property
    def d(self):
        return self.coords.d
    
    def __eq__(self, other):
        return self.coords == other.coords
        
    def __lt__(self, other):
        return self.coords < other.coords

    def __gt__(self, other):
        return self.coords > other.coords
    
    def __hash__(self):
        return hash(self.coords)
    
    def __mul__(self, scalar):
        """Return vector (self) * scalar."""
        newcoords = [scalar * dim for dim in self.coords]
        return type(self)(newcoords)

    __rmul__ = __mul__ # allow scalar * vector

    def __truediv__(self,scalar):
        """Return vector (self) * 1/scalar"""        
        return self.__mul__(1/scalar)
    
    def __add__(self,v1):
        """Add a vector to this vector, return a vector""" 
        newcoords = tuple(map(add, v1.coords, self.coords))
        return type(self)(newcoords)
        
    def __sub__(self,v1):
        """Subtract vector from this vector, return a vector"""
        return self.__add__(-v1)
    
    def __neg__(self):      
        """Return a vector, the negative of this one."""
        return type(self)(tuple(map(neg, self.coords)))
                  
    def length(self):
        """
        Uses norm0
        """
        t = self.norm0()
        return sp.sqrt(half * (t[0]**2 + t[1]**2 + t[2]**2 + t[3]**2))
        
    def cross(self,v1):
        """Return the cross product of self with another vector.
        return a Qvector"""
        A = type(self)((one,zero,zero,zero))
        B = type(self)((zero,one,zero,zero))
        C = type(self)((zero,zero,one,zero))
        D = type(self)((zero,zero,zero,one))
        a1,b1,c1,d1 = v1.coords
        a2,b2,c2,d2 = self.coords
        k= root2/4
        the_sum =   (A*c1*d2 - A*d1*c2 - A*b1*d2 + A*b1*c2
               + A*b2*d1 - A*b2*c1 - B*c1*d2 + B*d1*c2 
               + b1*C*d2 - b1*D*c2 - b2*C*d1 + b2*D*c1 
               + a1*B*d2 - a1*B*c2 - a1*C*d2 + a1*D*c2
               + a1*b2*C - a1*b2*D - a2*B*d1 + a2*B*c1 
               + a2*C*d1 - a2*D*c1 - a2*b1*C + a2*b1*D)
        return k * the_sum
    
    def area(self, v1):
        """
        area in unit triangles of edges D
        """
        return self.cross(v1).length() * 2/root3

    def angle(self, v1):
        return self.xyz.angle(v1.xyz)
        
    @property
    def xyz(self):
        a,b,c,d     =  self.coords
        k           =  1/(2 * root2)
        xyz         = (k * (a - b - c + d),
                       k * (a - b + c - d),
                       k * (a + b - c - d))
        return Vector(xyz)
        
class Svector(Vector):
    """Subclass of Vector that takes spherical coordinate args."""
    
    def __init__(self,arg):
        # if returning from Vector calc method, spherical is true
        arg = Vector(arg).spherical()
            
        # initialize a vector at an (r,phi,theta) tuple (= arg)
        r     = arg[0]
        phi   = radians(arg[1])
        theta = radians(arg[2])
        self.coords = tuple(map(lambda x:round(x,15),
                      (r * cos(theta) * sin(phi),
                       r * sin(theta) * sin(phi),
                       r * cos(phi))))
        self.xyz = self.coords

    def __repr__(self):
        return "Svector " + str(self.spherical())

def dot(a,b):
    return a.dot(b)

def cross(a,b):
    return a.cross(b)

def angle(a,b):
    return a.angle(b)

def length(a):
    return a.length()

def ivm_basis():
    a = Qvector((one, zero, zero, zero)) 
    b = Qvector((zero, one, zero, zero))
    c = Qvector((zero, zero, one, zero))
    d = Qvector((zero, zero, zero, one))
    return a, b, c, d

def xyz_basis():
    x = Vector((one, zero, zero))
    y = Vector((zero, one, zero))
    z = Vector((zero, zero, one))
    return x, y, z

def test1():
    a, b, c, d = ivm_basis() # unpacking assignment
    print("Qvector length (symbolic ) :", a.length())
    print("Qvector length (evaluated):" , a.length().evalf(50))
    print("Vector length  (symbolic ) :", a.xyz.length())
    print("Vector length  (evaluated) :", a.xyz.length().evalf(50))

def test2():
    b0, b1, b2, b3 = ivm_basis() # unpacking assignment
    print("IVM -> XYZ")
    print("{:32}{:56}{}".format("IVM", "XYZ", "Length"))
    for v in b0, b1, b2, b3:
        xyz = v.xyz
        print("{}\t{}\t{}".format(v, xyz, xyz.length()))

def test3():
    x, y, z = xyz_basis() # unpacking assignment
    print("XYZ -> IVM")
    print("{:32}{:48}{}".format("XYZ", "IVM", "Length"))
    for v in x, y, z:
        q = v.quadray()
        print("{}\t{}\t{}".format(v, q, q.length()))
    
if __name__ == "__main__":
    test1()