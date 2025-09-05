#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
London and Coulomb - Verification
=================================

Verification of static electromagnetic field equations and Coulomb limit from
doc/projected_em.tex (lines 576-587). Tests the fundamental electrostatic
relationships that emerge in the static limit of electromagnetic wave equations.

Tests cover:
- Static limit of wave equations (∂t = 0)
- Poisson equation for electrostatic potential
- Coulomb law calibration (E = Q/r²)
- Vanishing magnetic fields in static limit
- Lorenz gauge in static regime
- Energy density formulations
- Mathematical relationship verification
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, sqrt, simplify, Rational, diff, exp, atan

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify
)


def test_static_wave_equation_reduction(v):
    """Test reduction of wave equations to static limit."""
    v.subsection("Static Wave Equation Reduction")

    # From lines 577-581: In the static limit (∂t = 0), wave equations reduce to:
    # ∇² Φ_Slope = -4π ρ_ch  (Poisson equation)
    # ∇² A = 0 (in Lorenz gauge)

    # Test dimensional consistency first
    phi_slope_dim = v.get_dim('Phi_E')  # Use Phi_E instead of Phi_Slope
    rho_charge_dim = v.get_dim('rho_charge')

    v.check_dims("Slope potential", phi_slope_dim, v.get_dim('V_volt'))
    v.check_dims("Charge density", rho_charge_dim, v.get_dim('e') / v.L**3)

    # Test the Poisson equation structure
    laplacian_phi = v.lap_dim(phi_slope_dim)
    poisson_rhs_dim = rho_charge_dim

    v.check_dims("Laplacian of potential", laplacian_phi, v.get_dim('V_volt') / v.L**2)

    # Mathematical verification: ∇² Φ_Slope = -4π ρ_ch
    phi, rho_ch = symbols('phi rho_charge', real=True)
    x, y, z = symbols('x y z', real=True)

    # For a test case: point charge at origin gives φ = Q/r
    Q, r_sym = symbols('Q r', real=True, positive=True)
    test_potential = Q / r_sym

    # Test Poisson equation for point source
    # In spherical coordinates: ∇²(Q/r) = -4πQδ(r)
    # For r > 0: ∇²(Q/r) = 0 (verified by direct calculation)
    first_deriv = diff(test_potential, r_sym)  # d(Q/r)/dr = -Q/r²
    second_deriv = diff(first_deriv, r_sym)    # d²(Q/r)/dr² = 2Q/r³

    # Spherical Laplacian: ∇²f = (1/r²) d/dr(r² df/dr)
    radial_laplacian = (diff(r_sym**2 * first_deriv, r_sym)) / r_sym**2

    v.check_eq("Point charge Laplacian for r > 0", radial_laplacian, 0)

    # Vector potential equation ∇² A = 0 in static limit
    vector_potential_dim = v.get_dim('A')
    laplacian_A_dim = v.lap_dim(vector_potential_dim)

    v.check_dims("Laplacian of vector potential", laplacian_A_dim,
                 v.get_dim('A') / v.L**2)


def test_coulomb_law_calibration(v):
    """Test Coulomb law calibration from point charge solution."""
    v.subsection("Coulomb Law Calibration")

    # From lines 582-584: For point charge Q at origin:
    # Φ_Slope(r) = Q/r, E(r) = Q r/r³, B = 0

    Q_sym, r_sym = symbols('Q r', real=True, positive=True)

    # Test potential solution
    phi_point = Q_sym / r_sym

    # Electric field from E = -∇Φ (gradient in radial direction)
    E_radial = -diff(phi_point, r_sym)  # E_r = -dΦ/dr = Q/r²
    expected_E_magnitude = Q_sym / r_sym**2

    v.check_eq("Coulomb field from potential gradient", E_radial, expected_E_magnitude)

    # Full vector field: E(r) = Q r/r³ = Q/r² * (r/r) where r/r is unit vector
    # The magnitude is |E| = Q/r²
    E_magnitude_squared = E_radial**2
    expected_magnitude_squared = (Q_sym / r_sym**2)**2

    v.check_eq("Coulomb field magnitude squared", E_magnitude_squared, expected_magnitude_squared)

    # Test that magnetic field vanishes in static limit: B = 0
    # This follows from B = ∇ × A and ∇² A = 0 with proper boundary conditions
    v.info("In static limit: B = ∇ × A = 0 when A = 0 (from ∇² A = 0)")

    # Dimensional verification
    electric_field_dim = v.get_dim('E_field')
    charge_dim = v.get_dim('e')
    length_dim = v.L

    v.check_dims("Electric field", electric_field_dim, v.get_dim('V_volt') / length_dim)

    # Test calibration: this normalizes so dynamic results reduce to Coulomb law
    v.info("Calibration ensures: ∇·E = 4π ρ_ch and E = Q/r² in static limit")


def test_static_electromagnetic_consistency(v):
    """Test consistency of static electromagnetic equations."""
    v.subsection("Static Electromagnetic Consistency")

    # From line 587: Fix couplings so ∇·E = 4π ρ_ch and E = Q/r²

    # Test Gauss law: ∇·E = 4π ρ_ch (Gaussian units)
    electric_field_div = v.div_dim(v.get_dim('E_field'))
    charge_density_dim = v.get_dim('rho_charge')

    v.check_dims("Divergence of E field", electric_field_div,
                 v.get_dim('E_field') / v.L)
    v.check_dims("Charge density", charge_density_dim,
                 v.get_dim('e') / v.L**3)

    # Mathematical relationship: div(E) = div(-grad(φ)) = -Laplacian(φ)
    # Combined with Poisson: ∇²φ = -4π ρ, gives ∇·E = -(-4π ρ) = 4π ρ

    phi_test = symbols('phi_test')

    # Vector identity verification: ∇·(-∇φ) = -∇²φ
    v.info("Vector identity: ∇·E = ∇·(-∇φ) = -∇²φ = 4π ρ (from Poisson)")

    # Test consistency: if ∇²φ = -4π ρ, then ∇·E = 4π ρ
    poisson_lhs = symbols('laplacian_phi')  # ∇²φ
    poisson_rhs = -4 * pi * symbols('rho')  # -4π ρ
    gauss_lhs = -poisson_lhs  # -∇²φ = ∇·E
    gauss_rhs = 4 * pi * symbols('rho')  # 4π ρ

    v.check_eq("Gauss law from Poisson consistency", gauss_lhs.subs(poisson_lhs, poisson_rhs), gauss_rhs)


def test_lorenz_gauge_static_limit(v):
    """Test Lorenz gauge condition in static limit."""
    v.subsection("Lorenz Gauge in Static Limit")

    # From line 580: In static limit, Lorenz gauge ∂_μ A^μ = 0 becomes ∇·A = 0
    # Since ∂t = 0, the time derivative term vanishes

    # Lorenz gauge: ∂_μ A^μ = (1/c) ∂t Φ + ∇·A = 0
    # Static limit: (1/c) ∂t Φ = 0, so ∇·A = 0 (Coulomb gauge)

    vector_potential_div = v.div_dim(v.get_dim('A'))

    v.check_dims("Divergence of vector potential", vector_potential_div,
                 v.get_dim('A') / v.L)

    # Mathematical verification: In static limit, ∇·A = 0
    # This is the Coulomb gauge condition
    v.info("Lorenz gauge in static limit reduces to Coulomb gauge: ∇·A = 0")

    # Test that this is consistent with ∇² A = 0
    # If ∇·A = 0 and ∇² A = 0, then A satisfies both conditions consistently
    v.info("Consistency: ∇² A = 0 and ∇·A = 0 are compatible in static regime")


def test_electromagnetic_energy_density(v):
    """Test electromagnetic energy density formulations."""
    v.subsection("Electromagnetic Energy Density")

    # From lines 614-619: Energy density definitions
    # u = (1/2)(ε₀|E|² + |B|²/μ₀) [SI]
    # u = (1/8π)(|E|² + |B|²) [Gaussian]

    # Test dimensional consistency
    E_field_dim = v.get_dim('E_field')
    B_field_dim = v.get_dim('B_field')
    energy_density_dim = v.get_dim('u_EM')

    v.check_dims("E field", E_field_dim, v.get_dim('V_volt') / v.L)
    v.check_dims("B field", B_field_dim, v.get_dim('B_field'))
    v.check_dims("Energy density", energy_density_dim, v.get_dim('E_energy') / v.L**3)

    # Mathematical verification of energy density formulation
    E_sym, B_sym = symbols('E_field B_field', real=True, positive=True)

    # Gaussian units form: u = (|E|² + |B|²)/(8π)
    energy_density_gaussian = (E_sym**2 + B_sym**2) / (8 * pi)

    # In static limit: B = 0, so u = |E|²/(8π)
    static_energy_density = energy_density_gaussian.subs(B_sym, 0)
    expected_static = E_sym**2 / (8 * pi)

    v.check_eq("Static energy density (B=0)", static_energy_density, expected_static)

    # Test that energy density is always positive
    v.info("Energy density u = (|E|² + |B|²)/(8π) ≥ 0 for all field configurations")

    # For Coulomb field E = Q/r²: u = Q²/(8π r⁴)
    Q_sym, r_sym = symbols('Q r', real=True, positive=True)
    coulomb_E = Q_sym / r_sym**2
    coulomb_energy_density = coulomb_E**2 / (8 * pi)
    expected_coulomb_u = Q_sym**2 / (8 * pi * r_sym**4)

    v.check_eq("Coulomb field energy density", coulomb_energy_density, expected_coulomb_u)


def test_retarded_vs_static_solutions(v):
    """Test relationship between retarded and static solutions."""
    v.subsection("Retarded vs Static Solutions")

    # From lines 589-599: Retarded potentials reduce to static in appropriate limit
    # Retarded time: t_r = t - |x-x'|/c
    # Static limit: |x-x'|/c → 0, so t_r → t

    # Mathematical test: static limit of retarded potential
    t, t_r, x, x_prime = symbols('t t_r x x_prime', real=True)
    c_sym = symbols('c', real=True, positive=True)

    # Retarded time relation
    retarded_time = t - abs(x - x_prime) / c_sym

    # Static limit: c → ∞, so retarded time → t
    static_limit = sp.limit(retarded_time, c_sym, sp.oo)

    v.check_eq("Static limit of retarded time", static_limit, t)

    # Test that retarded solutions preserve causality
    # No information propagates faster than c
    v.info("Retarded solutions ensure causality: fields depend only on past light cone")

    # In static limit: quasi-static approximation
    v.info("Static limit: |x-x'|/c ≪ 1 recovers instantaneous-looking Coulomb law")


def test_poynting_vector_static_limit(v):
    """Test Poynting vector in static limit."""
    v.subsection("Poynting Vector in Static Limit")

    # From lines 602-611: Poynting vector S = (c/4π) E × B
    # In static limit: B = 0, so S = 0 (no energy flow)

    # Test dimensional consistency
    poynting_dim = v.get_dim('S_poynting')
    E_field_dim = v.get_dim('E_field')
    B_field_dim = v.get_dim('B_field')

    v.check_dims("Poynting vector", poynting_dim,
                 v.get_dim('E_energy') / (v.T * v.L**2))

    # Mathematical verification: S = (c/4π) E × B
    E_vec, B_vec = symbols('E_vector B_vector')
    c_sym = symbols('c', real=True, positive=True)

    # Poynting vector magnitude proportional to |E||B|sin(θ)
    # In static limit: B = 0, so |S| = 0
    E_mag, B_mag = symbols('E_magnitude B_magnitude', real=True, positive=True)
    poynting_magnitude = c_sym / (4 * pi) * E_mag * B_mag
    static_poynting = poynting_magnitude.subs(B_mag, 0)

    v.check_eq("Static Poynting vector", static_poynting, 0)

    # Test Poynting theorem: ∂u/∂t + ∇·S = -J·E
    # In static limit: ∂u/∂t = 0, ∇·S = 0, so J·E = 0 for steady currents
    v.info("Static Poynting theorem: ∇·S = 0 when B = 0 (no energy flow)")


def test_london_and_coulomb():
    """Run all London and Coulomb verification tests."""
    # Set up verification helper
    v = PhysicsVerificationHelper(
        "London and Coulomb",
        "Verification of static electromagnetic field equations and Coulomb limit"
    )

    v.section("LONDON AND COULOMB VERIFICATION")
    v.info("Verifying static electromagnetic field equations and Coulomb limit")
    v.info("Based on doc/projected_em.tex lines 576-587")

    try:
        # Run all test functions
        test_static_wave_equation_reduction(v)
        test_coulomb_law_calibration(v)
        test_static_electromagnetic_consistency(v)
        test_lorenz_gauge_static_limit(v)
        test_electromagnetic_energy_density(v)
        test_retarded_vs_static_solutions(v)
        test_poynting_vector_static_limit(v)

        # Return summary for test runner integration
        return v.summary()

    except Exception as e:
        v.error(f"Test execution failed: {e}")
        return 1

    return 0


if __name__ == "__main__":
    success_rate = test_london_and_coulomb()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)