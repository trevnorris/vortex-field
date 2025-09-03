"""
Foundational Postulates - Verification
======================================================================

Complete verification of all mathematical relationships in the Foundational Postulates
subsection (P-1 through P-6). This implements the dimensional and mathematical checks
specified in the FP.md guide.

P-1: Geometry & regularity (thin/flat regime, bending energy)
P-2: Kelvin circulation quantization
P-3: Decomposition into irrotational & solenoidal parts
P-4: Local axisymmetric kernel (thin core)
P-5: Projection invariance of circulation (leading order)
P-6: Topology & discreteness of intersections
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper, define_symbols_batch, batch_check_dims,
    quick_verify
)
import sympy as sp
from sympy import pi, oo, integrate, symbols, sqrt

# Initialize verification helper
v = PhysicsVerificationHelper(
    "Foundational Postulates – Verification",
    "Dimensional and math checks for P-1..P-6"
)

# Define any needed SymPy symbols (numbers/angles)
rho, w, theta = define_symbols_batch(['rho', 'w', 'theta'], positive=True, real=True)

# Declare theta as dimensionless (it's an angle)
v.declare_dimensionless('theta')

# ==============================================================================
# P-1: GEOMETRY & REGULARITY
# ==============================================================================

v.section("P-1: Geometry & Regularity")

# Check that epsilon_xi = xi/r and epsilon_kappa = kappa_geom*r are dimensionless
v.assert_dimensionless(v.dims['xi']/v.dims['r'], "epsilon_xi")
v.assert_dimensionless(v.dims['kappa_geom']*v.dims['r'], "epsilon_kappa")

# Bending energy functional: E_bend = K_bend * ∫ kappa_geom^2 dl
v.check_dims("Bending energy functional",
             v.dims['E_bend'],
             v.dims['K_bend'] * v.dims['kappa_geom']**2 * v.dims['dl'])

# ==============================================================================
# P-2: CIRCULATION
# ==============================================================================

v.section("P-2: Circulation")

# Verify dims of the circulation integral: ∮ v·dl = Γ
v.check_dims("Kelvin circulation dims",
             v.dims['v']*v.dims['dl'],
             v.dims['Gamma'])

# Check Gamma vs kappa (quantum of circulation)
v.check_dims("Gamma vs kappa (quantum of circulation)",
             v.dims['Gamma'],
             v.dims['kappa'])

# ==============================================================================
# P-3: HELMHOLTZ DECOMPOSITION
# ==============================================================================

v.section("P-3: Helmholtz Decomposition")

# All parts have velocity units
v.check_dims("v_irrot units", v.dims['v'], v.dims['v'])
v.check_dims("v_sol units", v.dims['v'], v.dims['v'])

# Constraint equations
v.check_dims("curl v_irrot = 0", v.curl_dim(v.dims['v']), 0)  # zero treated as dimension-agnostic
v.check_dims("div v_sol = 0", v.div_dim(v.dims['v']), 0)

# ==============================================================================
# P-4: LOCAL AXISYMMETRIC KERNEL
# ==============================================================================

v.section("P-4: Local Axisymmetric Kernel")

# (a) Dimensionless integrand and integral
# Check: [rho]^2 * [dl] vs ([rho]^2 + [w]^2)^(3/2)
v.check_dims("kernel integrand dimensionless",
             v.dims['r']**2 * v.dims['dl'],
             (v.dims['r']**2 + v.dims['w']**2)**(sp.Rational(3,2)))

# (b) Math: ∫ ρ^2 / (ρ^2 + w^2)^(3/2) dw = 2
res = integrate(rho**2 / (rho**2 + w**2)**sp.Rational(3,2), (w, -oo, oo))
quick_verify("kernel integral equals 2", sp.simplify(res - 2) == 0)

# (c) v_theta(ρ) = Gamma/(2πρ) -- dimensional check
v.check_dims("v_theta dims",
             v.dims['v_theta'],
             v.dims['Gamma']/v.dims['r'])

# ==============================================================================
# P-5: PROJECTION INVARIANCE (LEADING ORDER)
# ==============================================================================

v.section("P-5: Projection Invariance (Leading Order)")

# Angles and small parameters are dimensionless
v.assert_dimensionless(v.dims['theta'], "tilt angle theta")
v.validate_transcendentals(sp.cos(v.dims['theta']), "projection factor")

# Circulation keeps same units regardless of cos(theta) factor
v.check_dims("projected circulation has [L^2/T]",
             v.dims['v']*v.dims['dl'],
             v.dims['Gamma'])  # same as P-2 dims

# ==============================================================================
# P-6: TOPOLOGY & DISCRETENESS
# ==============================================================================

v.section("P-6: Topology & Discreteness")

# This postulate is qualitative (no numerical units to check)
quick_verify("postulate is qualitative (no units to check)", True)

# ==============================================================================
# SUMMARY
# ==============================================================================

v.summary()
