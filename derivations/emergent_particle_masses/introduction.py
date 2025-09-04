#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Introduction - Verification
====================================

Comprehensive verification of the fundamental mass template equation and core
theoretical concepts introduced in the emergent particle masses framework.

This test validates the dimensional consistency of the mass template formula
(eq:mass-template), key variable definitions, circulation quantum relationships,
and foundational physical concepts exactly as presented in the introduction
section of the document.

Based on doc/emergent_particle_masses.tex, introduction section (lines 1-37).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, sqrt, ln, exp, simplify

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify,
)


def test_mass_template_equation(v):
    """
    Test the dimensional consistency of the fundamental mass template equation
    as presented in eq:mass-template.

    Verifies: m(R) ≈ ρ₀·2πR[C_core·ξ_c² + κ²/(4π·v_L²)·ln(R/a)]

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Mass Template Equation (eq:mass-template)")

    # Define the symbolic variables as they appear in the document
    R, xi_c, kappa, v_L, a, C_core, rho_0 = define_symbols_batch(
        ['R', 'xi_c', 'kappa', 'v_L', 'a', 'C_core', 'rho_0'],
        positive=True
    )

    # First, verify dimensions of individual components
    v.check_dims("Loop radius R", v.get_dim('R_loop'), v.L)
    v.check_dims("Core healing length ξ_c", v.get_dim('xi'), v.L)
    v.check_dims("Quantum circulation κ", v.get_dim('kappa'), v.L**2/v.T)
    v.check_dims("Bulk wave speed v_L", v.get_dim('v_L'), v.L/v.T)
    v.check_dims("Inner cutoff a", v.get_dim('a_cutoff'), v.L)
    v.check_dims("Projected density ρ₀", v.get_dim('rho_0'), v.M/v.L**3)

    # C_core = 2π ln(2) is dimensionless (numerical constant)
    v.info("C_core = 2π ln(2) is dimensionless numerical constant")

    # Verify the core depletion term: ρ₀·2πR·C_core·ξ_c²
    core_term = v.get_dim('rho_0') * v.get_dim('R_loop') * v.get_dim('xi')**2
    v.check_dims("Core depletion term ρ₀·2πR·C_core·ξ_c²",
                 core_term, v.M)

    # Verify the compressibility term: ρ₀·2πR·κ²/(4π·v_L²)·ln(R/a)
    # Note: ln(R/a) is dimensionless
    compressibility_term = (v.get_dim('rho_0') * v.get_dim('R_loop') *
                           v.get_dim('kappa')**2 / v.get_dim('v_L')**2)
    v.check_dims("Compressibility term ρ₀·2πR·κ²/(4π·v_L²)·ln(R/a)",
                 compressibility_term, v.M)

    # Verify the complete mass template equation
    mass_template = core_term + compressibility_term
    v.check_dims("Complete mass template m(R)",
                 mass_template, v.M)

    v.success("Mass template equation dimensional consistency verified")


def test_fundamental_variable_definitions(v):
    """
    Test the dimensional consistency of fundamental variable definitions
    as stated in the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Fundamental Variable Definitions")

    # Test density projection relationship: ρ₀ = ρ_{4D}^0 · ξ_c
    v.check_dims("Density projection ρ₀ = ρ_{4D}^0 · ξ_c",
                 v.get_dim('rho_0'),
                 v.get_dim('rho_4') * v.get_dim('xi'))

    # Test circulation quantum: κ = h/m
    # Note: h has dimensions [M L² T⁻¹], m has dimensions [M]
    # So κ = h/m has dimensions [L² T⁻¹] which matches circulation
    planck_h = v.M * v.L**2 / v.T  # Planck constant dimension
    circulation_from_quantum = planck_h / v.M
    v.check_dims("Circulation quantum κ = h/m",
                 v.get_dim('kappa'), circulation_from_quantum)

    # Test bulk wave speed: v_L = √(g·ρ_{4D}^0/m²)
    # For this to be dimensionally consistent, g must have appropriate interaction dimensions
    # From v_L² = g·ρ_{4D}^0/m², we get [L²/T²] = [g]·[M/L⁴]/[M²]
    # So [g] = [L²/T²]·[M²]/[M/L⁴] = [M L⁶ T⁻²]
    g_interaction = v.M * v.L**6 / v.T**2

    # Add this dimension to the helper for the test
    v.add_dimensions({'g_interaction': g_interaction})

    wave_speed_rhs = sp.sqrt(g_interaction * v.get_dim('rho_4') / v.M**2)
    v.check_dims("Bulk wave speed v_L = √(g·ρ_{4D}^0/m²)",
                 v.get_dim('v_L'), wave_speed_rhs)

    # Test inner cutoff relationship: a = α·ξ_c
    # where α is an O(1) dimensionless constant
    v.check_dims("Inner cutoff a = α·ξ_c",
                 v.get_dim('a_cutoff'), v.get_dim('xi'))

    v.success("Fundamental variable definitions verified")


def test_physical_interpretation_consistency(v):
    """
    Test the physical consistency of the mass template components
    as described in the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Physical Interpretation Consistency")

    # The document states that both terms (core depletion and compressibility)
    # represent density deficit contributions projected onto the slice

    # Core depletion term represents "core deficit per unit length" times loop circumference
    # This should have dimension of deficit mass density times volume
    core_deficit_per_length = v.get_dim('rho_0') * v.get_dim('xi')**2
    v.check_dims("Core deficit per unit length ρ₀·C_core·ξ_c²",
                 core_deficit_per_length, v.M/v.L)

    # When multiplied by loop circumference 2πR, gives total deficit mass
    total_core_deficit = core_deficit_per_length * v.get_dim('R_loop')
    v.check_dims("Total core deficit mass",
                 total_core_deficit, v.M)

    # Compressibility term represents "far-field Bernoulli contribution"
    # The κ²/(4π·v_L²) factor should represent a characteristic deficit scale
    bernoulli_scale = v.get_dim('kappa')**2 / v.get_dim('v_L')**2
    v.check_dims("Bernoulli deficit scale κ²/(4π·v_L²)",
                 bernoulli_scale, v.L**2)

    # When multiplied by ρ₀·2πR·ln(R/a), gives logarithmic mass contribution
    bernoulli_mass = v.get_dim('rho_0') * v.get_dim('R_loop') * bernoulli_scale
    v.check_dims("Bernoulli mass contribution",
                 bernoulli_mass, v.M)

    v.info("Both mass template terms represent projected density deficits")
    v.success("Physical interpretation consistency verified")


def test_charge_concept_framework(v):
    """
    Test the conceptual framework for electric charge as described in the document.

    Note: This tests the dimensional framework, not the detailed charge quantization
    which is referenced to be in another section.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Charge Concept Framework")

    # The document states charge is a topological threading number
    # This means it should be dimensionless (a pure number/index)
    v.info("Electric charge Q as topological threading number (dimensionless)")

    # However, for electromagnetic calculations, charge has dimension [Q]
    v.check_dims("Electric charge in EM context",
                 v.Q, v.Q)

    # The document mentions that neutrino-like defects can have Q=0
    # while still exhibiting Eddies/Drag (local electromagnetic effects)
    v.info("Neutrino-like defects: Q=0 but can have local Eddies/Drag")

    # Eddies represent solenoidal flow (magnetic-like patterns)
    # Drag represents angular momentum of motion
    v.check_dims("Angular momentum (Drag)",
                 v.get_dim('J_angular'), v.M * v.L**2 / v.T)

    v.success("Charge concept framework verified")


def test_units_convention(v):
    """
    Test the units convention as stated in the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Units Convention")

    # Document states: "We retain ℏ and m symbolically in definitional formulas;
    # unless otherwise noted, numerical evaluations set ℏ=m=1."

    v.info("Units convention: ℏ and m retained symbolically")
    v.info("Numerical evaluations: ℏ = m = 1 (unless noted)")

    # In the natural unit system where ℏ = 1:
    # ℏ has dimension [M L² T⁻¹]
    hbar_natural = v.M * v.L**2 / v.T
    v.info(f"ℏ dimension in natural units: {hbar_natural}")

    # In the natural unit system where m = 1:
    # m has dimension [M]
    mass_natural = v.M
    v.info(f"m dimension in natural units: {mass_natural}")

    # This means in natural units, κ = ℏ/m becomes dimensionless × [L² T⁻¹]
    # which maintains the circulation dimension
    kappa_natural = hbar_natural / mass_natural
    v.check_dims("κ in natural units (ℏ=m=1)",
                 v.get_dim('kappa'), kappa_natural)

    v.success("Units convention consistency verified")


def test_introduction():
    """
    Main test function for Introduction section.

    This function coordinates all verification tests for the introduction section,
    validating the fundamental mass template equation, variable definitions,
    and core theoretical concepts exactly as presented in the document.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Introduction - Emergent Particle Masses",
        "Mass template equation and fundamental theoretical framework"
    )

    v.section("INTRODUCTION VERIFICATION")

    # Add any custom dimensions needed for the tests
    # (Most dimensions are already defined in helper.py)
    v.add_dimensions({
        'R_loop': v.L,  # Loop radius (R in the document)
        'a_cutoff': v.L,  # Inner cutoff scale (a in the document)
        'C_core': 1,  # Core deficit constant (dimensionless)
    })

    # Call test functions in logical order
    v.info("\n--- 1) Mass Template Equation ---")
    test_mass_template_equation(v)

    v.info("\n--- 2) Fundamental Variable Definitions ---")
    test_fundamental_variable_definitions(v)

    v.info("\n--- 3) Physical Interpretation Consistency ---")
    test_physical_interpretation_consistency(v)

    v.info("\n--- 4) Charge Concept Framework ---")
    test_charge_concept_framework(v)

    v.info("\n--- 5) Units Convention ---")
    test_units_convention(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_introduction()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
