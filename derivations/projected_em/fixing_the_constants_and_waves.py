#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "Fixing the Constants and Waves" subsection.

This module implements dimensional and mathematical verification for the subsection
covering electromagnetic wave equations, speed of light relationships, and the
fixing of electromagnetic constants μ₀ and ε₀.

Based on doc/projected_em.tex, subsection "Fixing the Constants and Waves"
(lines 185-201).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, sqrt

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    verify_wave_equation,
    quick_verify
)


def test_wave_equation_dimensional_consistency(v):
    """
    Test dimensional consistency of electromagnetic wave equations:
    (∇² - (1/c²)∂_tt)E = 0 and (∇² - (1/c²)∂_tt)B = 0
    """
    v.subsection("Wave Equation Dimensional Consistency")

    # Electric field wave equation: (∇² - (1/c²)∂_tt)E = 0
    laplacian_E = v.lap_dim(v.get_dim('E'))                    # ∇²E
    time_term_E = v.get_dim('E') / (v.get_dim('c')**2 * v.T**2)   # (1/c²)∂_tt E

    v.check_dims("E wave equation: ∇²E vs (1/c²)∂_tt E", laplacian_E, time_term_E)

    # Magnetic field wave equation: (∇² - (1/c²)∂_tt)B = 0
    laplacian_B = v.lap_dim(v.get_dim('B'))                    # ∇²B
    time_term_B = v.get_dim('B') / (v.get_dim('c')**2 * v.T**2)   # (1/c²)∂_tt B

    v.check_dims("B wave equation: ∇²B vs (1/c²)∂_tt B", laplacian_B, time_term_B)

    # Both wave equations should have same operator structure
    wave_operator_dim = laplacian_E  # This is the dimensional structure of the wave operator

    v.success("Wave equations are dimensionally consistent")


def test_speed_of_light_relationship(v):
    """
    Test the fundamental relationship c² = 1/(μ₀ε₀) and its dimensional consistency.
    """
    v.subsection("Speed of Light Fundamental Relationship")

    # c² should have dimensions [L² T⁻²]
    c_squared = v.get_dim('c')**2
    expected_c_squared = v.L**2 / v.T**2

    v.check_dims("c² dimensions", c_squared, expected_c_squared)

    # 1/(μ₀ε₀) should have same dimensions as c²
    mu0_eps0_product = v.get_dim('mu_0') * v.get_dim('epsilon_0')
    inverse_mu0_eps0 = 1 / mu0_eps0_product

    v.check_dims("c² = 1/(μ₀ε₀) relationship", c_squared, inverse_mu0_eps0)

    # Individual constant dimensions should be consistent
    # μ₀: [M L Q⁻²], ε₀: [M⁻¹ L⁻³ T² Q²]
    mu0_expected = v.M * v.L / v.Q**2
    eps0_expected = v.Q**2 * v.T**2 / (v.M * v.L**3)

    v.check_dims("μ₀ dimensions", v.get_dim('mu_0'), mu0_expected)
    v.check_dims("ε₀ dimensions", v.get_dim('epsilon_0'), eps0_expected)

    # Product should give [T² L⁻²]
    product_expected = v.T**2 / v.L**2
    v.check_dims("μ₀ε₀ product", mu0_eps0_product, product_expected)

    v.success("Speed of light relationship c² = 1/(μ₀ε₀) verified")


def test_maxwell_to_wave_derivation(v):
    """
    Test the dimensional consistency of the derivation from Maxwell equations to wave equations.
    The derivation involves taking curl of Ampère-Maxwell and using homogeneous equations.
    """
    v.subsection("Maxwell to Wave Equation Derivation")

    # Starting with Ampère-Maxwell: ∇×B - μ₀ε₀∂_t E = μ₀J
    # Taking curl: ∇×(∇×B) - μ₀ε₀∇×(∂_t E) = μ₀∇×J

    # ∇×(∇×B) term - using vector identity ∇×(∇×B) = ∇(∇·B) - ∇²B
    curl_curl_B = v.grad_dim(v.div_dim(v.get_dim('B'))) - v.lap_dim(v.get_dim('B'))

    # Since ∇·B = 0 (homogeneous equation), first term vanishes, leaving -∇²B
    minus_laplacian_B = -v.lap_dim(v.get_dim('B'))

    # ∇×(∂_t E) = ∂_t(∇×E) term
    curl_dt_E = v.curl_dim(v.dt(v.get_dim('E')))
    dt_curl_E = v.dt(v.curl_dim(v.get_dim('E')))

    v.check_dims("∇×(∂_t E) = ∂_t(∇×E) consistency", curl_dt_E, dt_curl_E)

    # Using homogeneous equation ∇×E + ∂_t B = 0, so ∇×E = -∂_t B
    # Therefore: ∂_t(∇×E) = ∂_t(-∂_t B) = -∂_tt B
    dt_curl_E_sub = -v.dtt(v.get_dim('B'))

    v.check_dims("∂_t(∇×E) = -∂_tt B substitution", dt_curl_E, dt_curl_E_sub)

    # Final wave equation term: -∇²B - μ₀ε₀(-∂_tt B) = μ₀∇×J
    # Rearranging: -∇²B + μ₀ε₀∂_tt B = μ₀∇×J
    # In vacuum (J=0): ∇²B - μ₀ε₀∂_tt B = 0
    # Using c² = 1/(μ₀ε₀): ∇²B - (1/c²)∂_tt B = 0

    vacuum_lhs = v.lap_dim(v.get_dim('B')) - v.dtt(v.get_dim('B')) / v.get_dim('c')**2

    # This should be dimensionally consistent (both terms same dimension)
    term1 = v.lap_dim(v.get_dim('B'))
    term2 = v.dtt(v.get_dim('B')) / v.get_dim('c')**2

    v.check_dims("Wave equation terms: ∇²B vs (1/c²)∂_tt B", term1, term2)

    v.success("Maxwell to wave equation derivation verified")


def test_electromagnetic_constant_relationships(v):
    """
    Test relationships between electromagnetic constants and their physical meaning.
    """
    v.subsection("Electromagnetic Constant Relationships")

    # Impedance of free space Z₀ = √(μ₀/ε₀)
    Z0_from_constants = sqrt(v.get_dim('mu_0') / v.get_dim('epsilon_0'))
    Z0_expected = v.get_dim('Z_0')  # Should be defined in helper.py

    v.check_dims("Free space impedance Z₀ = √(μ₀/ε₀)", Z0_from_constants, Z0_expected)

    # Z₀ should also equal √(μ₀/ε₀) = μ₀c = 1/(ε₀c)
    mu0_c = v.get_dim('mu_0') * v.get_dim('c')
    inv_eps0_c = 1 / (v.get_dim('epsilon_0') * v.get_dim('c'))

    v.check_dims("Z₀ = μ₀c relationship", Z0_expected, mu0_c)
    v.check_dims("Z₀ = 1/(ε₀c) relationship", Z0_expected, inv_eps0_c)

    # Energy density relationships
    # Electric energy density: (1/2)ε₀E²
    E_energy_density = v.get_dim('epsilon_0') * v.get_dim('E')**2

    # Magnetic energy density: (1/2)B²/μ₀
    B_energy_density = v.get_dim('B')**2 / v.get_dim('mu_0')

    # Both should have energy density dimensions [M L⁻¹ T⁻²]
    energy_density_expected = v.M / (v.L * v.T**2)

    v.check_dims("Electric energy density (1/2)ε₀E²", E_energy_density, energy_density_expected)
    v.check_dims("Magnetic energy density (1/2)B²/μ₀", B_energy_density, energy_density_expected)

    # In plane waves, electric and magnetic energy densities are equal
    v.check_dims("EM energy density equality", E_energy_density, B_energy_density)

    v.success("Electromagnetic constant relationships verified")


def test_wave_operator_consistency(v):
    """
    Test consistency of the wave operator (∇² - (1/c²)∂_tt) across different fields.
    """
    v.subsection("Wave Operator Consistency")

    # Define the wave operator dimension structure
    # ∇² has dimension [L⁻²], ∂_tt has dimension [T⁻²]
    # So (1/c²)∂_tt has dimension [T⁻²]/[L²T⁻²] = [L⁻²]

    laplacian_dim = v.L**(-2)
    time_factor_dim = 1 / (v.get_dim('c')**2 * v.T**2)

    v.check_dims("Wave operator terms: ∇² vs (1/c²)∂_tt", laplacian_dim, time_factor_dim)

    # The wave operator should work on any field with proper dimensions
    # Test with scalar potential Φ (if it satisfies wave equation)
    phi_wave_term1 = v.lap_dim(v.get_dim('Phi'))
    phi_wave_term2 = v.get_dim('Phi') / (v.get_dim('c')**2 * v.T**2)

    v.check_dims("Scalar potential wave terms", phi_wave_term1, phi_wave_term2)

    # Test with vector potential A
    A_wave_term1 = v.lap_dim(v.get_dim('A'))
    A_wave_term2 = v.get_dim('A') / (v.get_dim('c')**2 * v.T**2)

    v.check_dims("Vector potential wave terms", A_wave_term1, A_wave_term2)

    v.success("Wave operator consistency verified across different fields")


def test_physical_interpretation_and_constraints(v):
    """
    Test the physical interpretation: measuring c and Coulomb force fixes both μ₀ and ε₀.
    """
    v.subsection("Physical Constraints and Interpretation")

    # The speed c is measured experimentally
    # The Coulomb force gives us information about ε₀
    # Coulomb's law: F = (1/4πε₀) * (q₁q₂/r²)

    # Force has dimensions [M L T⁻²]
    # Charges have dimensions [Q]
    # Distance has dimensions [L]
    # So 1/ε₀ should have dimensions to make this work

    coulomb_factor = 1 / (v.get_dim('epsilon_0') * v.L**2)  # (1/4πε₀r²) factor
    force_per_charge_squared = v.M * v.L / (v.T**2 * v.Q**2)  # F/(q₁q₂)

    v.check_dims("Coulomb law factor 1/(4πε₀r²)", coulomb_factor, force_per_charge_squared)

    # Once ε₀ is fixed by Coulomb measurements and c is measured,
    # then μ₀ follows from c² = 1/(μ₀ε₀)
    mu0_from_c_and_eps0 = 1 / (v.get_dim('c')**2 * v.get_dim('epsilon_0'))

    v.check_dims("μ₀ from c² = 1/(μ₀ε₀)", v.get_dim('mu_0'), mu0_from_c_and_eps0)

    # Self-consistency check: the product μ₀ε₀ should give 1/c²
    product_check = v.get_dim('mu_0') * v.get_dim('epsilon_0')
    expected_product = 1 / v.get_dim('c')**2

    v.check_dims("Self-consistency: μ₀ε₀ = 1/c²", product_check, expected_product)

    v.success("Physical interpretation and constraint relationships verified")


def test_fixing_constants_and_waves():
    """
    Main test function for the "Fixing the Constants and Waves" subsection.

    Tests all aspects of electromagnetic wave equations, speed of light relationships,
    and the fixing of electromagnetic constants.
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Fixing the Constants and Waves",
        "Electromagnetic wave equations and constant relationships verification"
    )

    v.section("FIXING THE CONSTANTS AND WAVES VERIFICATION")

    # Test 1: Wave equation dimensional consistency
    v.info("\n--- 1) Wave Equation Dimensional Consistency ---")
    test_wave_equation_dimensional_consistency(v)

    # Test 2: Speed of light fundamental relationship
    v.info("\n--- 2) Speed of Light Relationship c² = 1/(μ₀ε₀) ---")
    test_speed_of_light_relationship(v)

    # Test 3: Maxwell to wave equation derivation
    v.info("\n--- 3) Maxwell Equations to Wave Equations Derivation ---")
    test_maxwell_to_wave_derivation(v)

    # Test 4: Electromagnetic constant relationships
    v.info("\n--- 4) Electromagnetic Constant Relationships ---")
    test_electromagnetic_constant_relationships(v)

    # Test 5: Wave operator consistency
    v.info("\n--- 5) Wave Operator Consistency ---")
    test_wave_operator_consistency(v)

    # Test 6: Physical interpretation and constraints
    v.info("\n--- 6) Physical Interpretation and Constraints ---")
    test_physical_interpretation_and_constraints(v)

    # Final summary
    return v.summary()


if __name__ == "__main__":
    test_fixing_constants_and_waves()
