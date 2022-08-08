"""Theoretical models for studying polymers in free or confining medium, in crowded or crowder-free media based on the concept of Kuhn segment.
"""
from typing import Callable
import sympy as sp
import numpy as np

def w_ao(
    r: float,
    phi_c: float,
    a_c: float,
    a: float = 1.0, 
    beta: float = 1.0
) -> Callable:
    """Return the symbolic form of the Asakura-Oosawa depletion potential.

    Parameters
    ----------
    r: center to center distance: sympy symbol
    vfrc_c: the volume fraction of crowders: float
    a_c: the crowder diameter: float
    a: the monomer diameter: float
    beta: the inverse of thermal energy - beta= 1/ k_B*T where k_b is the Boltzmann constant 
    and T is absolute temperature: float

    Returns:
    symbolic ao potential

    Dependencies:
    Sympy

    Idea:
    phi_c symbol or float

    """
    ao = 1.0 * phi_c * ((a + a_c)**3 / a_c**3) * (
        1- ((3 * r)/ (2 * (a + a_c))) + (r**3 / (2 * (a + a_c)**3))
        )
    return ao

def v_wca(r

, sigma, rcut, epsilon=1.0):
    """
    Return the symbolic form of the fully-repulsive Weeks-Chandler-Anderson (WCA) potential.


    Parameters:
    r: center to center distance: simpy symbol
    sigma: particle 1 diameter: float
    rcut: the cuttoff distance: float # the WCA cutoff distance
    espsilon: strength of the atteractive well: float

    Return:
    symbolic wca potential

    Dependencies:
    Sympy

    Idea to imporve:
    mechanisms to internally choose sigma and rcut.
    """
    wca = sp.Piecewise((4 * epsilon * ((sigma**12 / r**12) - (sigma**6 / r**6) + 1/4),r<=rcut),(0,r>rcut))
    return wca

def mayer_f(potential, beta=1.0):
    """
    return the sympolic form of the Mayer-f function.

    Parameters:
    potential: the potential of interest: sympy symbolic function 
    beta: the inverse of thermal energy - beta= 1/ k_B*T where k_b is the Boltzmann constant 
    and T is absolute temperature: float

    Return:
    symbolic Mayer-f function

    Dependencies:
    Sympy

    """
    return sp.exp(-1 * beta * potential) - 1.0

def exc_vol(r, mayer_f):

    """
    return the value of the excluded volume of a monomer with Mayer-f function mayer_f

    Parameters:
    mayer_f: the symbolic form of thr Mayer-f function constructed from a pair potential.

    Return:
    a tuple: numerical value of the excluded value (float) and its error (float)

    Dependencies:
    Scipi
    Sympi

    Ideas:
    This is dangerous to converting the symbolic function to numeric function using lambdify.

    """

    integrand = -4.0 * sp.pi * r**2 * mayer_f # the integrand of the excluded volume integral
    integrand = sp.lambdify(r,integrand, "numpy") # converting symbolic form to numeric one
    return quad(integrand,0,np.inf)

# coefficients
def coeffs(a_c, a):
    """
    Return empirical coefficients for computing the excluded volume of a monomer.

    The excluded volume of a monomer is given by a fitting function to a FENE polymer in LJ-crowdered bulk space.

    Reference:
    https://doi.org/10.1039/C6SM01184E

    Parameters: 
    a: float 
    diameter (size) of a monomer
    a_c: float
    diameter (size) of a crowder

    Returns: 
    c: ndarray
    the coefficient of the excluded volume function

    Dependencies:
    numpy package

    """
    q_inverse = a / a_c

    if q_inverse > 1.0: # a_c < a (q<1):
        b1 = -1.60
        b2 = 2.92
        b3 = -2.65
        c = np.array([b1 * q_inverse, b2 * q_inverse**2, b3 * q_inverse**3])
    else : ## a_c > a (q>1):
        c1 = -5.04
        c2 = 11.90
        c3 = -11.40
        c = np.array([c1, c2, c3])
    return c