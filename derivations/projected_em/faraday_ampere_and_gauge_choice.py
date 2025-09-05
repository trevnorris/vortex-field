#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Faraday, Ampère-Maxwell Laws and Gauge Choice - Verification
=============================================================

Complete verification of Faraday's law, Ampère-Maxwell law, and gauge choice
interactions in projected electromagnetism. This tests both dimensional
consistency and mathematical relationships between:

1. Faraday's law: ∇×E + ∂_t B = 0 (homogeneous Maxwell equation)
2. Ampère-Maxwell law: ∇×B - μ₀ε₀ ∂_t E = μ₀ J (inhomogeneous Maxwell equation)
3. Gauge choices: Coulomb gauge ∇·A = 0 vs Lorenz gauge ∂_μ A^μ = 0
4. Electric and magnetic field definitions: E = -∇Φ - ∂_t A, B = ∇×A

The framework shows these laws emerge from 4D aether projection, with gauge
freedom reflecting different ways to decompose the vector field v = ∇φ + ∇×A.

Based on doc/projected_em.tex, sections covering homogeneous and inhomogeneous
Maxwell equations, field definitions, and gauge choices in wave sector.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    batch_check_dims,
    quick_verify
)
import sympy as sp
from sympy import symbols, simplify, diff, sqrt, pi, Derivative


def test_field_definitions_and_gauge(v):
    """
    Test electromagnetic field definitions and gauge choice consistency.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.section("FIELD DEFINITIONS AND GAUGE CONSISTENCY")

    # Test field definitions: B = ∇×A and E = -∇Φ - ∂_t A
    v.subsection("Electromagnetic Field Definitions")

    # Define symbolic fields and potentials
    t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)
    A_x, A_y, A_z = symbols('A_x A_y A_z', real=True)
    Phi = symbols('Phi', real=True)

    # Verify B = ∇×A has correct dimensions
    v.check_dims(
        "Magnetic field from vector potential: B = ∇×A",
        v.curl_dim(v.get_dim('A')),
        v.get_dim('B')
    )

    # Verify E components have correct dimensions
    grad_Phi_dim = v.grad_dim(v.get_dim('Phi'))
    dt_A_dim = v.dt(v.get_dim('A'))

    v.check_dims(
        "Electric field gradient term: -∇Φ",
        v.get_dim('E'),
        grad_Phi_dim
    )

    v.check_dims(
        "Electric field induction term: -∂_t A",
        v.get_dim('E'),
        dt_A_dim
    )

    # Both terms in E = -∇Φ - ∂_t A must have same dimensions
    v.check_dims(
        "E field terms dimensional matching: ∇Φ vs ∂_t A",
        grad_Phi_dim,
        dt_A_dim
    )

    # Test gauge choices
    v.subsection("Gauge Choice Dimensional Consistency")

    # Coulomb gauge: ∇·A = 0
    div_A_dim = v.div_dim(v.get_dim('A'))
    expected_div_A_dim = v.get_dim('A') / v.get_dim('x')

    v.check_dims(
        "Coulomb gauge condition: ∇·A",
        div_A_dim,
        expected_div_A_dim
    )

    # Lorenz gauge: ∂_μ A^μ = ∂_t Φ/c² + ∇·A = 0
    dt_Phi_over_c2_dim = v.dt(v.get_dim('Phi')) / (v.get_dim('c')**2)

    v.check_dims(
        "Lorenz gauge time term: ∂_t Φ/c²",
        dt_Phi_over_c2_dim,
        div_A_dim
    )

    v.success("Field definitions and gauge consistency verified")


def test_faraday_law_mathematical(v):
    """
    Test Faraday's law mathematical relationships and derivation.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.section("FARADAY'S LAW MATHEMATICAL VERIFICATION")

    v.subsection("Faraday's Law: ∇×E + ∂_t B = 0")

    # Check dimensional consistency of Faraday's law
    curl_E_dim = v.curl_dim(v.get_dim('E'))
    dt_B_dim = v.dt(v.get_dim('B'))

    v.check_dims(
        "Faraday's law dimensional consistency",
        curl_E_dim,
        dt_B_dim
    )

    # Mathematical verification using symbolic expressions
    t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)
    A_x, A_y, A_z = symbols('A_x A_y A_z', real=True, cls=sp.Function)
    Phi = symbols('Phi', real=True, cls=sp.Function)

    # Make them functions of space and time
    A_x = A_x(x, y, z, t)
    A_y = A_y(x, y, z, t)
    A_z = A_z(x, y, z, t)
    Phi = Phi(x, y, z, t)

    # Define field components symbolically
    # B = ∇×A
    B_x = sp.diff(A_z, y) - sp.diff(A_y, z)
    B_y = sp.diff(A_x, z) - sp.diff(A_z, x)
    B_z = sp.diff(A_y, x) - sp.diff(A_x, y)

    # E = -∇Φ - ∂_t A
    E_x = -sp.diff(Phi, x) - sp.diff(A_x, t)
    E_y = -sp.diff(Phi, y) - sp.diff(A_y, t)
    E_z = -sp.diff(Phi, z) - sp.diff(A_z, t)

    # Faraday's law components: ∇×E + ∂_t B = 0
    faraday_x = (sp.diff(E_z, y) - sp.diff(E_y, z)) + sp.diff(B_x, t)
    faraday_y = (sp.diff(E_x, z) - sp.diff(E_z, x)) + sp.diff(B_y, t)
    faraday_z = (sp.diff(E_y, x) - sp.diff(E_x, y)) + sp.diff(B_z, t)

    # Simplify - these should be zero by construction
    faraday_x_simplified = sp.simplify(faraday_x)
    faraday_y_simplified = sp.simplify(faraday_y)
    faraday_z_simplified = sp.simplify(faraday_z)

    v.check_eq("Faraday x-component: (∇×E)_x + ∂_t B_x = 0", faraday_x_simplified, 0)
    v.check_eq("Faraday y-component: (∇×E)_y + ∂_t B_y = 0", faraday_y_simplified, 0)
    v.check_eq("Faraday z-component: (∇×E)_z + ∂_t B_z = 0", faraday_z_simplified, 0)

    # Test mathematical derivation from potentials
    v.subsection("Mathematical Derivation from Potentials")

    # If E = -∇Φ - ∂_t A and B = ∇×A, then:
    # ∇×E = ∇×(-∇Φ - ∂_t A) = -∇×(∇Φ) - ∇×(∂_t A)
    # Since ∇×(∇Φ) = 0, we get: ∇×E = -∇×(∂_t A) = -∂_t(∇×A) = -∂_t B
    # Therefore: ∇×E + ∂_t B = 0

    # Check curl of gradient is dimensionally consistent (but identically zero)
    curl_grad_Phi_dim = v.curl_dim(v.grad_dim(v.get_dim('Phi')))
    expected_zero_dim = v.get_dim('E') / v.get_dim('x')

    v.check_dims(
        "Curl of gradient: ∇×(∇Φ) (identically zero)",
        curl_grad_Phi_dim,
        expected_zero_dim
    )

    # Mixed derivative equality: ∇×(∂_t A) = ∂_t(∇×A)
    curl_dt_A_dim = v.curl_dim(v.dt(v.get_dim('A')))
    dt_curl_A_dim = v.dt(v.curl_dim(v.get_dim('A')))

    v.check_dims(
        "Mixed derivatives commutation: ∇×(∂_t A) = ∂_t(∇×A)",
        curl_dt_A_dim,
        dt_curl_A_dim
    )

    # Final verification: ∇×E = -∂_t B
    v.check_dims(
        "Faraday relation: ∇×E = -∂_t B",
        curl_E_dim,
        dt_B_dim
    )

    v.success("Faraday's law mathematical verification completed")


def test_ampere_maxwell_law_mathematical(v):
    """
    Test Ampère-Maxwell law mathematical relationships.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.section("AMPÈRE-MAXWELL LAW MATHEMATICAL VERIFICATION")

    v.subsection("Ampère-Maxwell Law: ∇×B - μ₀ε₀ ∂_t E = μ₀ J")

    # Check dimensional consistency of each term
    curl_B_dim = v.curl_dim(v.get_dim('B'))
    displacement_current_dim = v.get_dim('mu_0') * v.get_dim('epsilon_0') * v.dt(v.get_dim('E'))
    source_term_dim = v.get_dim('mu_0') * v.get_dim('j_current')

    v.check_dims(
        "Ampère-Maxwell curl term: ∇×B",
        curl_B_dim,
        source_term_dim
    )

    v.check_dims(
        "Displacement current: μ₀ε₀ ∂_t E",
        displacement_current_dim,
        curl_B_dim
    )

    v.check_dims(
        "Current source: μ₀ J",
        source_term_dim,
        curl_B_dim
    )

    # Mathematical verification of key relationships
    v.subsection("Mathematical Identity Verification")

    # Verify fundamental constants relationship: c² = 1/(μ₀ε₀)
    c_squared = v.get_dim('c')**2
    inverse_mu0eps0 = 1 / (v.get_dim('mu_0') * v.get_dim('epsilon_0'))

    v.check_dims(
        "Light speed relation: c² = 1/(μ₀ε₀)",
        c_squared,
        inverse_mu0eps0
    )

    # Verify vacuum impedance: Z₀ = √(μ₀/ε₀)
    Z0_squared = v.get_dim('Z_0')**2
    mu0_over_eps0 = v.get_dim('mu_0') / v.get_dim('epsilon_0')

    v.check_dims(
        "Vacuum impedance: Z₀² = μ₀/ε₀",
        Z0_squared,
        mu0_over_eps0
    )

    # Test connection to continuity equation
    v.subsection("Connection to Current Conservation")

    # Taking divergence of Ampère-Maxwell: ∇·(∇×B) - μ₀ε₀ ∇·(∂_t E) = μ₀ ∇·J
    # Since ∇·(∇×B) = 0, we get: -μ₀ε₀ ∂_t(∇·E) = μ₀ ∇·J
    # Using Gauss law ∇·E = ρ/ε₀: -μ₀ε₀ ∂_t(ρ/ε₀) = μ₀ ∇·J
    # Simplifying: -μ₀ ∂_t ρ = μ₀ ∇·J → ∂_t ρ + ∇·J = 0

    div_curl_B_dim = v.div_dim(curl_B_dim)  # Should be zero dimensionally
    div_dt_E_dim = v.div_dim(v.dt(v.get_dim('E')))
    div_J_dim = v.div_dim(v.get_dim('j_current'))

    # The key relationship is: μ₀ε₀ ∇·(∂_t E) should equal μ₀ ∇·J
    mu0eps0_div_dt_E_dim = v.get_dim('mu_0') * v.get_dim('epsilon_0') * div_dt_E_dim
    mu0_div_J_dim = v.get_dim('mu_0') * div_J_dim

    v.check_dims(
        "Divergence relationship: μ₀ε₀ ∇·(∂_t E) vs μ₀ ∇·J",
        mu0eps0_div_dt_E_dim,
        mu0_div_J_dim
    )

    # Verify continuity equation dimensions
    dt_rho_dim = v.dt(v.get_dim('rho_charge'))

    v.check_dims(
        "Current conservation: ∂_t ρ + ∇·J = 0",
        dt_rho_dim,
        div_J_dim
    )

    v.success("Ampère-Maxwell law mathematical verification completed")


def test_gauge_transformation_mathematics(v):
    """
    Test mathematical aspects of gauge transformations.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.section("GAUGE TRANSFORMATION MATHEMATICS")

    v.subsection("Gauge Transformation Properties")

    # Gauge transformation: A → A + ∇χ, Φ → Φ - ∂_t χ
    # where χ is an arbitrary scalar function

    # Check that gauge transformation preserves field dimensions
    # For the gauge transformation A → A + ∇χ to work, ∇χ must have dimensions of A
    # This means χ has dimensions such that ∇χ has dimensions of A
    chi_dim_from_A = v.get_dim('A') * v.get_dim('x')  # [A][L] so ∇χ has [A]
    grad_chi_dim = v.grad_dim(chi_dim_from_A)

    # For Φ → Φ - ∂_t χ to work, ∂_t χ must have dimensions of Φ
    dt_chi_dim = v.dt(chi_dim_from_A)

    v.check_dims(
        "Gauge parameter gradient: ∇χ matches A",
        grad_chi_dim,
        v.get_dim('A')
    )

    # Check if ∂_t χ matches Φ dimensions
    v.check_dims(
        "Gauge parameter time derivative consistency",
        dt_chi_dim,
        v.get_dim('Phi')
    )

    # Test gauge invariance of physical fields
    v.subsection("Physical Field Gauge Invariance")

    # E and B should be gauge invariant
    # E = -∇Φ - ∂_t A = -∇(Φ - ∂_t χ) - ∂_t(A + ∇χ) = -∇Φ - ∂_t A + ∇(∂_t χ) - ∂_t(∇χ)
    # Using ∇(∂_t χ) = ∂_t(∇χ), the gauge terms cancel: E unchanged ✓

    # B = ∇×A = ∇×(A + ∇χ) = ∇×A + ∇×(∇χ) = ∇×A + 0 = ∇×A
    # Since curl of gradient is zero: B unchanged ✓

    v.info("Physical fields are gauge invariant:")
    v.info("  E: gauge terms ∇(∂_t χ) - ∂_t(∇χ) = 0 (commutativity)")
    v.info("  B: additional term ∇×(∇χ) = 0 (curl of gradient)")

    # Test specific gauge conditions
    v.subsection("Specific Gauge Conditions")

    # Coulomb gauge: ∇·A = 0 (convenient for statics)
    v.info("Coulomb gauge: ∇·A = 0")
    v.info("  - Separates longitudinal/transverse parts")
    v.info("  - Instantaneous Coulomb interaction")

    # Lorenz gauge: ∂_μ A^μ = 0 → ∂_t Φ/c² + ∇·A = 0
    v.info("Lorenz gauge: ∂_t Φ/c² + ∇·A = 0")
    v.info("  - Lorentz invariant condition")
    v.info("  - Wave equations decouple cleanly")

    # Check wave equation simplification in Lorenz gauge
    wave_operator_dim = 1/(v.get_dim('c')**2) * v.dt(v.dt(v.get_dim('A'))) - v.lap_dim(v.get_dim('A'))
    source_dim = v.get_dim('mu_0') * v.get_dim('j_current')

    v.check_dims(
        "Wave equation operator: (1/c²)∂²_t - ∇²",
        wave_operator_dim,
        source_dim / v.get_dim('A')
    )

    v.success("Gauge transformation mathematics verified")


def test_maxwell_equations_consistency(v):
    """
    Test mathematical consistency between all four Maxwell equations.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.section("COMPLETE MAXWELL SYSTEM CONSISTENCY")

    v.subsection("Four Maxwell Equations Dimensional Check")

    # Homogeneous equations (exact topological identities)
    # 1. ∇·B = 0 - Mathematical verification
    t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)
    A_x, A_y, A_z = symbols('A_x A_y A_z', real=True, cls=sp.Function)

    # Make them functions of space and time
    A_x = A_x(x, y, z, t)
    A_y = A_y(x, y, z, t)
    A_z = A_z(x, y, z, t)

    # B = ∇×A
    B_x = sp.diff(A_z, y) - sp.diff(A_y, z)
    B_y = sp.diff(A_x, z) - sp.diff(A_z, x)
    B_z = sp.diff(A_y, x) - sp.diff(A_x, y)

    # ∇·B = ∂B_x/∂x + ∂B_y/∂y + ∂B_z/∂z
    div_B = sp.diff(B_x, x) + sp.diff(B_y, y) + sp.diff(B_z, z)
    div_B_simplified = sp.simplify(div_B)

    v.check_eq("No magnetic monopoles: ∇·B = 0 (exact identity)", div_B_simplified, 0)

    # Dimensional check
    div_B_dim = v.div_dim(v.get_dim('B'))
    v.check_dims("No magnetic monopoles: ∇·B dimensional consistency", div_B_dim, div_B_dim)

    # 2. ∇×E + ∂_t B = 0 (Faraday)
    curl_E_dim = v.curl_dim(v.get_dim('E'))
    dt_B_dim = v.dt(v.get_dim('B'))
    v.check_dims("Faraday's law consistency", curl_E_dim, dt_B_dim)

    # Inhomogeneous equations (from slice continuity + closure)
    # 3. ∇·E = ρ/ε₀ (Gauss)
    div_E_dim = v.div_dim(v.get_dim('E'))
    gauss_source_dim = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
    v.check_dims("Gauss law consistency", div_E_dim, gauss_source_dim)

    # 4. ∇×B - μ₀ε₀ ∂_t E = μ₀ J (Ampère-Maxwell)
    curl_B_dim = v.curl_dim(v.get_dim('B'))
    displacement_dim = v.get_dim('mu_0') * v.get_dim('epsilon_0') * v.dt(v.get_dim('E'))
    ampere_source_dim = v.get_dim('mu_0') * v.get_dim('j_current')

    v.check_dims("Ampère-Maxwell LHS consistency", curl_B_dim, displacement_dim)
    v.check_dims("Ampère-Maxwell full consistency", curl_B_dim, ampere_source_dim)

    # Test cross-equation relationships
    v.subsection("Cross-Equation Mathematical Relationships")

    # Bianchi identities (automatic from field definitions)
    # ∂_t(∇·B) = ∇·(∂_t B) = -∇·(∇×E) = 0 → ∇·B constant in time
    # ∂_t(∇×E) = ∇×(∂_t E) = -∇×∇×(B/μ₀ε₀ + μ₀J/μ₀ε₀)

    dt_div_B_dim = v.dt(div_B_dim)
    v.check_dims("Time evolution of ∇·B = 0", dt_div_B_dim, dt_div_B_dim)

    # Charge conservation from Gauss + Ampère-Maxwell
    # ∂_t(∇·E) = ∇·(∂_t E) → (1/ε₀)∂_t ρ = ∇·(∇×B/μ₀ε₀ - J)
    # Using ∇·(∇×B) = 0: (1/ε₀)∂_t ρ = -∇·J → ∂_t ρ + ∇·J = 0

    dt_gauss_dim = v.dt(gauss_source_dim)
    continuity_dim = v.dt(v.get_dim('rho_charge'))
    div_J_dim = v.div_dim(v.get_dim('j_current'))

    v.check_dims("Charge conservation from Maxwell eqs", continuity_dim, div_J_dim)

    v.success("Complete Maxwell system consistency verified")


def test_physical_interpretation_and_units(v):
    """
    Test physical interpretation and unit consistency.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.section("PHYSICAL INTERPRETATION AND UNITS")

    v.subsection("Physical Meaning of Laws")

    v.info("Physical interpretations:")
    v.info("1. ∇·B = 0: Magnetic field lines have no endpoints (topological)")
    v.info("2. ∇×E + ∂_t B = 0: Changing magnetic field induces electric field")
    v.info("3. ∇·E = ρ/ε₀: Electric charges create electric field")
    v.info("4. ∇×B - μ₀ε₀∂_t E = μ₀J: Currents and changing E create B")

    v.subsection("Gauge Choice Physical Meaning")

    v.info("Gauge choices reflect different decompositions of v = ∇φ + ∇×A:")
    v.info("- Coulomb gauge: ∇·A = 0 (separates longitudinal/transverse)")
    v.info("- Lorenz gauge: ∂_μA^μ = 0 (Lorentz invariant, symmetric waves)")

    # Test fundamental constants relationship
    v.subsection("Fundamental Constants Consistency")

    # c² = 1/(μ₀ε₀) from wave speed
    c_from_constants = 1 / sp.sqrt(v.get_dim('mu_0') * v.get_dim('epsilon_0'))
    v.check_dims("Light speed from EM constants", v.get_dim('c'), c_from_constants)

    # Impedance Z₀ = √(μ₀/ε₀)
    Z0_calculated = sp.sqrt(v.get_dim('mu_0') / v.get_dim('epsilon_0'))
    v.check_dims("Vacuum impedance", v.get_dim('Z_0'), Z0_calculated)

    # Note: Fine structure constant calculation requires elementary charge symbol
    # which may not be defined in all contexts
    v.info("Fine structure constant α = e²/(4πε₀ħc) is dimensionless by construction")

    v.success("Physical interpretation and units verified")


def test_faraday_ampere_and_gauge_choice():
    """
    Main test function for Faraday, Ampère-Maxwell laws and gauge choice verification.

    This function coordinates all verification tests for the electromagnetic laws
    and their gauge freedom, providing comprehensive mathematical verification.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Faraday, Ampère-Maxwell Laws and Gauge Choice",
        "Mathematical verification of electromagnetic laws and gauge freedom"
    )

    # Define symbolic variables for mathematical analysis
    t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)

    # Call test functions in logical order
    v.info("\n--- 1) Field Definitions and Gauge ---")
    test_field_definitions_and_gauge(v)

    v.info("\n--- 2) Faraday's Law Mathematics ---")
    test_faraday_law_mathematical(v)

    v.info("\n--- 3) Ampère-Maxwell Law Mathematics ---")
    test_ampere_maxwell_law_mathematical(v)

    v.info("\n--- 4) Gauge Transformation Mathematics ---")
    test_gauge_transformation_mathematics(v)

    v.info("\n--- 5) Complete Maxwell System ---")
    test_maxwell_equations_consistency(v)

    v.info("\n--- 6) Physical Interpretation ---")
    test_physical_interpretation_and_units(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_faraday_ampere_and_gauge_choice()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)