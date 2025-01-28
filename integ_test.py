from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from math import factorial

hbar = 1
m = 1
C = 1


TRACE_BOUNDARY = 50.  # hartree units


@dataclass
class Vec3:
    x: float
    y: float
    z: float

    def magnitude(self) -> float:
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    def __add__(self, other: "Vec3") -> "Vec3":
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: "Vec3") -> "Vec3":
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float) -> "Vec3":
        return Vec3(self.x * scalar, self.y * scalar, self.z * scalar)

    def __rmul__(self, scalar: float):
        return self.__mul__(scalar)

def norm_term(n: int, l: int) -> float:
    """
    Calculates the normalization term for given n and l.
    """
    return 2.  # todo temp for n=1; QC and put back the rest.

    nf = float(n)

    # Calculate the numerator of the normalization term
    norm_term_num = (2.0 / nf) ** 3 * factorial(n - l - 1)

    # Calculate the denominator of the normalization term
    norm_term_denom = 2 * n * factorial(n + l) ** 3

    # Return the square root of the ratio
    return sqrt(norm_term_num / norm_term_denom)


def sto(posit_sample: Vec3, posit_nuc: Vec3, n: int, l: int) -> float:
    # todo: Radial only; n=1 only
    diff = posit_sample - posit_charge
    r = diff.magnitude()

    nf = float(n)

#     exp_term = exp(-self.xi * r / (nf * A_0))
    exp_term = np.exp(-r / nf)  # todo temp

    # Laguerre polynomial term
    n_L = n - l - 1
    a = 2 * l + 1
    #     L = util.make_laguerre(n_L, a)
    L = lambda x: 1.0

    # Polynomial term
    polynomial_term = 2.0 * r / nf ** l * L(2.0 * r / nf)

    return norm_term(n, l) * polynomial_term * exp_term


# Single-point Coulomb potential, from a single point charge.
def V_coulomb(posit_charge: Vec3, posit_sample: Vec3, charge: float) -> float:
    diff = posit_sample - posit_charge
    r = diff.magnitude()

#     K_C * charge / (r + SOFTENING_FACTOR)
    K_C * charge / r



def local_approx(x: float, x_ref: float, psi: float, psi_pp: float) -> float:
    x_diff = x - x_ref
    lin_term = 1.

#     print(f"x: {x_diff}, psi: {psi}, psi'': {psi_pp}'")

    return 0.5 * psi_pp * x_diff**2 + psi + lin_term * x_diff


# Should be cplx.
# At a given pt.
def lagrangian(psi_pp: float, posit_trace: Vec3, posit_nuc: Vec3) -> float:
    diff = posit_trace - posit_nuc

    V_energy = 1. / diff.magnitude()

    -hbar ** 2 / (2.0*m) * psi_pp - V_energy


# Relativistic variant. See Feynman V2, Chp 19
def lagrangian_rel(psi_pp: float, posit_trace: Vec3, posit_nuc: Vec3) -> float:
    # Phi is the scalar V. A is the vector V(related to magnetism?)
    -m0 * C**2 * np.sqrt(1. - v**2 / c**2) - q * phi(posit_time) - v * A(posit_time)


def action_rel():
    pass
#      -m0 * C**2 * integ(a) - q * integ(b)


def action():
    # todo - important: Is this action approach equivalent to the comparison of psi'' and
    # todo psi you've been using for a while?
    # todo Perhaps even if so, this could help: Use established techniques for variational
    # todo calculus and principle of Least Action that may be novel to your prev approach.

    # todo: You must  figure out what the trace is. A path through, x, y, z, with
    # todo: an associated psi value at each step??
    # Your classical lagrangian is failing here, as is your KE *operator* term.

    #trace = [Vec3(-TRACE_BOUNDARY, -TRACE_BOUNDARY, -TRACE_BOUNDARY)]

    # todo: Try both a 1D trace (A slice of the WF), and a 3D volume covering the whole WF.

    range_min, range_max = -5, 5
    points_per_side = 20

    # Generate linearly spaced points for x, y, and z
    grid = np.linspace(range_min, range_max, points_per_side)

    # Create the 3D array of Vec3 objects
    posits = np.empty((points_per_side, points_per_side, points_per_side), dtype=object)
    psi = np.empty((points_per_side, points_per_side, points_per_side), dtype=np.float64)

    for i, x in enumerate(grid):
        for j, y in enumerate(grid):
            for k, z in enumerate(grid):
                posits[i, j, k] = Vec3(x, y, z)
                psi[i, j, k] = 0.

    posit_nuc = Vec3(0., 0., 0.)

    dx = 1.  # todo

    # S = 0.
    S = Vec3(0., 0., 0.)

    for i in range(posits.size):
        for j in range(posits.size):
            for k in range(posits.size):
                psi_pp = Vec3()
                S += lagrangian(psi_pp, posits[i, j, k], posit_nuc)

    S *= dx

    for posit_trace in trace:
        psi_pp = 0.  # todo
        S += lagrangian(psi_pp, posit_trace, posit_nuc) * dx

    plt.plot(trace)


def integ_test():
        x = np.linspace(0, 100, 1000)
        psi = np.zeros(x.size)

        psi[0] = -2  # initial condition

        for i in range(x.size - 1):
            x_ref = x[i]
            psi_pp = 1 * psi[i]  # Simple Poisson diffeq
            psi[i+1] = local_approx(x[i+1], x_ref, psi[i], psi_pp)

        plt.title("Integration of a second order ODE")

        plt.plot(x, psi)
        plt.show()


def main():
    action()


main()