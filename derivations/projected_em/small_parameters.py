#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "Small Parameters and Regime of Validity" subsection.

This module implements dimensional and mathematical verification for the small
parameters and regime of validity as outlined in the projected EM section.
Tests the four dimensionless expansion parameters and their physical consistency.

Based on doc/projected_em.tex, subsection "Small Parameters and Regime of Validity".
"""

import os
import sys
import sympy as sp
from sympy import symbols, sqrt, simplify, pi

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    quick_verify
)


def test_custom_dimensions_setup(v):
    """Add custom dimensions needed for small parameter tests."""
    # Add dimensions not already in helper.py
    custom_dims = {
        'ell': v.L,                    # On-slice length scale of interest
        'l_planck': v.L,               # Planck length scale
        'l_TP': v.L,                   # Alternative notation for Planck length
        'kappa': v.L**(-1),            # Local curvature of slice (L^-1)
    }

    v.add_dimensions(custom_dims, allow_overwrite=True)
    v.success("Custom dimensions added for small parameters")


def test_fundamental_dimensions(v):
    """Verify that fundamental quantities have correct dimensions."""
    v.subsection("Fundamental Dimension Verification")

    # Length scales should have length dimension
    v.check_dims("xi (healing length)", v.get_dim('xi'), v.L)
    v.check_dims("ell (on-slice scale)", v.get_dim('ell'), v.L)
    v.check_dims("l_TP (Planck length)", v.get_dim('l_TP'), v.L)

    # Velocities should have velocity dimension
    expected_velocity = v.L / v.T
    v.check_dims("v (characteristic velocity)", v.get_dim('v'), expected_velocity)
    v.check_dims("c (speed of light)", v.get_dim('c'), expected_velocity)

    # Curvature should have inverse length dimension
    v.check_dims("kappa (slice curvature)", v.get_dim('kappa'), v.L**(-1))

    v.success("Fundamental dimensions verified")


def test_small_parameter_epsilon_rho(v):
    """Test thinness parameter ε_ρ = ξ/ℓ is dimensionless."""
    v.subsection("Thinness Parameter ε_ρ")

    # ε_ρ := ξ/ℓ should be dimensionless
    epsilon_rho = v.get_dim('xi') / v.get_dim('ell')
    v.assert_dimensionless(epsilon_rho, "epsilon_rho = xi/ell")

    v.success("Thinness parameter ε_ρ verified as dimensionless")


def test_small_parameter_epsilon_v(v):
    """Test velocity parameter ε_v = v/c is dimensionless."""
    v.subsection("Velocity Parameter ε_v")

    # ε_v := v/c should be dimensionless
    epsilon_v = v.get_dim('v') / v.get_dim('c')
    v.assert_dimensionless(epsilon_v, "epsilon_v = v/c")

    v.success("Velocity parameter ε_v verified as dimensionless")


def test_small_parameter_epsilon_xi(v):
    """Test Planck scale parameter ε_ξ = ℓ_TP/ℓ is dimensionless."""
    v.subsection("Planck Scale Parameter ε_ξ")

    # ε_ξ := ℓ_TP/ℓ should be dimensionless
    epsilon_xi = v.get_dim('l_TP') / v.get_dim('ell')
    v.assert_dimensionless(epsilon_xi, "epsilon_xi = l_TP/ell")

    v.success("Planck scale parameter ε_ξ verified as dimensionless")


def test_small_parameter_epsilon_kappa(v):
    """Test curvature parameter ε_κ = κℓ is dimensionless."""
    v.subsection("Curvature Parameter ε_κ")

    # ε_κ := κℓ should be dimensionless
    epsilon_kappa = v.get_dim('kappa') * v.get_dim('ell')
    v.assert_dimensionless(epsilon_kappa, "epsilon_kappa = kappa*ell")

    v.success("Curvature parameter ε_κ verified as dimensionless")


def test_all_small_parameters_together(v):
    """Test all four small parameters together."""
    v.subsection("Combined Small Parameter Verification")

    # Define all four small parameters
    epsilon_rho = v.get_dim('xi') / v.get_dim('ell')
    epsilon_v = v.get_dim('v') / v.get_dim('c')
    epsilon_xi = v.get_dim('l_TP') / v.get_dim('ell')
    epsilon_kappa = v.get_dim('kappa') * v.get_dim('ell')

    # Test each individually (redundant but clear)
    v.assert_dimensionless(epsilon_rho, "ε_ρ in combined test")
    v.assert_dimensionless(epsilon_v, "ε_v in combined test")
    v.assert_dimensionless(epsilon_xi, "ε_ξ in combined test")
    v.assert_dimensionless(epsilon_kappa, "ε_κ in combined test")

    # Test combined remainder term O(ε_ρ² + ε_v² + ε_ξ² + ε_κ²)
    remainder_term = epsilon_rho**2 + epsilon_v**2 + epsilon_xi**2 + epsilon_kappa**2
    v.assert_dimensionless(remainder_term, "O(ε_ρ² + ε_v² + ε_ξ² + ε_κ²) remainder")

    v.success("All four small parameters verified together")


def test_notation_distinction(v):
    """Test distinction between ε_ρ (geometric) and δρ/ρ_0 (density contrast)."""
    v.subsection("Notation Distinction Verification")

    # ε_ρ is geometric thinness ξ/ℓ (already tested)
    epsilon_rho_geometric = v.get_dim('xi') / v.get_dim('ell')
    v.assert_dimensionless(epsilon_rho_geometric, "ε_ρ geometric thinness")

    # δρ/ρ_0 is density contrast amplitude (when needed in wave-sector)
    density_contrast = v.get_dim('rho') / v.get_dim('rho_0')  # Both should be M/L^3
    v.assert_dimensionless(density_contrast, "δρ/ρ_0 density contrast")

    v.success("Notation distinction between ε_ρ and δρ/ρ_0 verified")


def test_regime_validity_conditions(v):
    """Test that the thin-slow-flat slice limit makes physical sense."""
    v.subsection("Regime Validity Conditions")

    # In the valid regime, all small parameters should be << 1
    # This is more of a physical constraint than mathematical, but we can test structure

    # Test that the regime assumption is dimensionally consistent:
    # We assume ξ << ℓ, v << c, ℓ_TP << ℓ, κℓ << 1

    # These are ordering relations, so we test their dimensional consistency
    length_scale_ordering = v.get_dim('xi') - v.get_dim('ell')  # Should have length dimension
    v.check_dims("Length scale difference ξ - ℓ", length_scale_ordering, v.L)

    velocity_scale_ordering = v.get_dim('v') - v.get_dim('c')  # Should have velocity dimension
    v.check_dims("Velocity scale difference v - c", velocity_scale_ordering, v.L/v.T)

    planck_scale_ordering = v.get_dim('l_TP') - v.get_dim('ell')  # Should have length dimension
    v.check_dims("Planck scale difference ℓ_TP - ℓ", planck_scale_ordering, v.L)

    v.success("Regime validity dimensional structure verified")


def test_leading_order_approximation(v):
    """Test that field equations hold to leading order with controlled remainders."""
    v.subsection("Leading Order Approximation Structure")

    # The text states: "All field equations in this section hold to leading order
    # in these small numbers; we indicate remainders schematically as
    # O(ε_ρ² + ε_v² + ε_ξ² + ε_κ²)"

    # Test that the remainder structure is dimensionally consistent
    epsilon_rho = v.get_dim('xi') / v.get_dim('ell')
    epsilon_v = v.get_dim('v') / v.get_dim('c')
    epsilon_xi = v.get_dim('l_TP') / v.get_dim('ell')
    epsilon_kappa = v.get_dim('kappa') * v.get_dim('ell')

    # Individual squared terms
    v.assert_dimensionless(epsilon_rho**2, "ε_ρ² term")
    v.assert_dimensionless(epsilon_v**2, "ε_v² term")
    v.assert_dimensionless(epsilon_xi**2, "ε_ξ² term")
    v.assert_dimensionless(epsilon_kappa**2, "ε_κ² term")

    # Full remainder expression
    full_remainder = epsilon_rho**2 + epsilon_v**2 + epsilon_xi**2 + epsilon_kappa**2
    v.assert_dimensionless(full_remainder, "Full O(ε²) remainder")

    v.success("Leading order approximation structure verified")


def test_small_parameters_validity():
    """
    Main test function for Small Parameters and Regime of Validity subsection.

    Comprehensive verification of:
    A) Custom dimension setup
    B) Fundamental dimension verification
    C) Individual small parameter tests (4 parameters)
    D) Combined parameter verification
    E) Notation distinction tests
    F) Regime validity conditions
    G) Leading order approximation structure
    """
    v = PhysicsVerificationHelper(
        "Small Parameters and Regime of Validity",
        "Verification of dimensionless expansion parameters in thin-slow-flat limit"
    )

    v.section("SMALL PARAMETERS AND REGIME OF VALIDITY VERIFICATION")

    # A) Setup custom dimensions
    v.section_header("A) Custom Dimension Setup")
    test_custom_dimensions_setup(v)

    # B) Verify fundamental dimensions
    v.section_header("B) Fundamental Dimension Verification")
    test_fundamental_dimensions(v)

    # C) Individual small parameter tests
    v.section_header("C) Individual Small Parameter Tests")
    test_small_parameter_epsilon_rho(v)
    test_small_parameter_epsilon_v(v)
    test_small_parameter_epsilon_xi(v)
    test_small_parameter_epsilon_kappa(v)

    # D) Combined parameter verification
    v.section_header("D) Combined Parameter Verification")
    test_all_small_parameters_together(v)

    # E) Notation distinction
    v.section_header("E) Notation Distinction")
    test_notation_distinction(v)

    # F) Regime validity conditions
    v.section_header("F) Regime Validity Conditions")
    test_regime_validity_conditions(v)

    # G) Leading order approximation
    v.section_header("G) Leading Order Approximation")
    test_leading_order_approximation(v)

    # Final summary
    return v.summary()


if __name__ == "__main__":
    test_small_parameters_validity()
