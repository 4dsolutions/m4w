"""Here is using Perplexity to research Kirby’s repositories

https://www.perplexity.ai/search/Research-all-of-q7eQENnPT.eaEiaDyYfbZg

These functions purport to provide the following functionality:

1. quadray_to_xyz: Converts a quadray (a, b, c, d) to Cartesian coordinates (x, y, z).
2. xyz_to_quadray: Converts Cartesian coordinates (x, y, z) to a quadray (a, b, c, d).
3. quadray_midpoint: Finds the midpoint between two quadrays.
4. quadray_distance: Calculates the Euclidean distance between two quadrays.
5. generate_fcc_quadrays: Generates the first n layers of the face-centered cubic (FCC) lattice using quadrays.

Analysis:
https://coda.io/d/Math4Wisdom_d0SvdI3KSto/Synergetics-Language-Code-Models_suW8a#_luU8e

Sept 4, 2024:

Looking back, using AI to make quadrays code more perspicacious was 
in lieu of making the original codebase more intuitive.  Rewritten
xyz <-> abcd conversion algorithms clear up a lot of the confusions.
I recommend not doing anything with the code below. The original 
codebase (qrays.py) is what to use from this repo.
"""

from math import sqrt

def quadray_to_xyz(quadray):
    """
    Convert a quadray (a, b, c, d) to Cartesian coordinates (x, y, z).
    """
    a, b, c, d = quadray
    x = (a + b) / 2
    y = (b + c) / 2
    z = (c + d) / 2
    return x, y, z

def xyz_to_quadray(x, y, z):
    """
    Convert Cartesian coordinates (x, y, z) to a quadray (a, b, c, d).
    """
    a = x - y + z
    b = x + y - z
    c = -x + y + z
    d = -x - y - z
    return a, b, c, d

def quadray_midpoint(q1, q2):
    """
    Find the midpoint between two quadrays.
    """
    return tuple((a + b) / 2 for a, b in zip(q1, q2))

def quadray_distance(q1, q2):
    """
    Calculate the Euclidean distance between two quadrays.
    """
    return sqrt(sum((a - b)**2 for a, b in zip(q1, q2)))

def generate_fcc_quadrays(n):
    """
    Generate the first n layers of the face-centered cubic (FCC) lattice using quadrays.
    """
    fcc_quadrays = [(0, 0, 0, 0)]
    for layer in range(1, n + 1):
        for a in range(layer + 1):
            for b in range(layer + 1 - a):
                c = layer - a - b
                d = -layer
                fcc_quadrays.append((a, b, c, d))
                fcc_quadrays.append((a, b, d, c))
                fcc_quadrays.append((a, c, b, d))
                fcc_quadrays.append((a, c, d, b))
                fcc_quadrays.append((a, d, b, c))
                fcc_quadrays.append((a, d, c, b))
    return fcc_quadrays

