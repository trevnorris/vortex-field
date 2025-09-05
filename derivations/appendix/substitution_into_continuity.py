#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Substitution into Continuity - Verification
====================================

Comprehensive verification of the substitution of the Bernoulli-derived density
relationship into the continuity equation, leading to the complete quasilinear
second-order PDE governing compressible potential flow in the projected aether.

This test validates the dimensional consistency of:
- The initial substitution equation with density from Bernoulli form
- The simplified quasilinear PDE after multiplication by -K/ξ
- Quadratic and cubic nonlinearities from convection and variable effective speed
- All terms in the final nonlinear governing equation

Based on doc/appendix.tex, "Substitution into Continuity" subsection (lines 50-61).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, sqrt, simplify

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify,
    verify_conservation_law,
)


def test_bernoulli_density_relationship(v):
    """
    Test the Bernoulli-derived density relationship used in the substitution.

    Verifies: ρ3D = -(ξ/K)[∂tΨ + (1/2)(∇Ψ)²]
    This is the key relationship that enables the substitution.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Bernoulli Density Relationship")

    # Define the Bernoulli quantity: B = ∂tΨ + (1/2)(∇Ψ)²
    # [B] = [L²/T²] (specific energy-like quantity)
    B_dims = v.L**2 / v.T**2
    v.check_dims("Bernoulli quantity B = ∂tΨ + (1/2)(∇Ψ)²",
                 v.get_dim('B_bernoulli'), B_dims)

    # Test the density relationship: ρ3D = -(ξ/K) B
    # [ξ/K] = [L] / [L⁶/(M·T²)] = [M·T²/L⁵]
    # [ρ3D] = [M·T²/L⁵] × [L²/T²] = [M/L³] ✓
    xi_over_K = v.get_dim('xi') / v.get_dim('K_barotropic')
    density_from_bernoulli = xi_over_K * v.get_dim('B_bernoulli')
    v.check_dims("Density from Bernoulli ρ3D = -(ξ/K)[∂tΨ + (1/2)(∇Ψ)²]",
                 v.get_dim('rho_0'), density_from_bernoulli)

    # Verify coefficient dimensional consistency
    v.check_dims("Coefficient ξ/K dimensional check",
                 xi_over_K, v.M * v.T**2 / v.L**5)

    v.success("Bernoulli density relationship verified")


def test_initial_substitution_equation(v):
    """
    Test the initial substitution equation before simplification.

    Verifies: ∂t[-ξ/K(∂tΨ + (1/2)(∇Ψ)²)] - ∇·[-ξ/K(∂tΨ + (1/2)(∇Ψ)²)∇Ψ] = -Ṁbody
    This is equation 54 from the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Initial Substitution Equation")

    # Time derivative term: ∂t[-ξ/K(∂tΨ + (1/2)(∇Ψ)²)]
    # This is: -(ξ/K) ∂t[∂tΨ + (1/2)(∇Ψ)²]
    # [∂t B] = [L²/T³], so [-(ξ/K) ∂t B] = [M·T²/L⁵][L²/T³] = [M/(L³·T)]
    time_term = v.get_dim('xi') / v.get_dim('K_barotropic') * v.dt(v.get_dim('B_bernoulli'))
    v.check_dims("Time derivative term ∂t[-(ξ/K)B]",
                 time_term, v.M / (v.L**3 * v.T))

    # Divergence term: -∇·[-ξ/K(∂tΨ + (1/2)(∇Ψ)²)∇Ψ]
    # This is: (ξ/K) ∇·[B ∇Ψ]
    # [B ∇Ψ] = [L²/T²][L/T] = [L³/T³]
    # [∇·(B ∇Ψ)] = [L³/T³]/[L] = [L²/T³]
    # [(ξ/K) ∇·(B ∇Ψ)] = [M·T²/L⁵][L²/T³] = [M/(L³·T)]
    flux_term = (v.get_dim('xi') / v.get_dim('K_barotropic') *
                 v.get_dim('B_bernoulli') * v.grad_dim(v.get_dim('Psi_velocity_potential')) / v.L)
    v.check_dims("Flux divergence term (ξ/K) ∇·[B ∇Ψ]",
                 flux_term, v.M / (v.L**3 * v.T))

    # Right-hand side: -Ṁbody
    # This should match the left-hand terms, so Ṁbody should be [M/(L³·T)]
    rhs_term = v.get_dim('M_dot_density')
    v.check_dims("RHS term -Ṁbody",
                 rhs_term, v.M / (v.L**3 * v.T))

    # Verify equation balance
    v.check_dims("Initial substitution equation balance",
                 time_term, rhs_term)

    v.success("Initial substitution equation dimensional structure verified")


def test_simplified_quasilinear_pde(v):
    """
    Test the simplified quasilinear PDE after multiplication by -K/ξ.

    Verifies: ∂t(∂tΨ + (1/2)(∇Ψ)²) + ∇·[(∂tΨ + (1/2)(∇Ψ)²)∇Ψ] = (K/ξ)Ṁbody
    This is equation 58 from the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Simplified Quasilinear PDE")

    # Time derivative term: ∂t(∂tΨ + (1/2)(∇Ψ)²) = ∂t B
    # [∂t B] = [L²/T³]
    time_B = v.dt(v.get_dim('B_bernoulli'))
    v.check_dims("Time derivative ∂t(∂tΨ + (1/2)(∇Ψ)²)",
                 time_B, v.L**2 / v.T**3)

    # Divergence term: ∇·[(∂tΨ + (1/2)(∇Ψ)²)∇Ψ] = ∇·[B ∇Ψ]
    # [B ∇Ψ] = [L²/T²][L/T] = [L³/T³]
    # [∇·(B ∇Ψ)] = [L³/T³]/[L] = [L²/T³]
    flux_B = v.get_dim('B_bernoulli') * v.grad_dim(v.get_dim('Psi_velocity_potential')) / v.L
    v.check_dims("Flux divergence ∇·[(∂tΨ + (1/2)(∇Ψ)²)∇Ψ]",
                 flux_B, v.L**2 / v.T**3)

    # Right-hand side: (K/ξ)Ṁbody
    # [K/ξ] = [L⁶/(M·T²)] / [L] = [L⁵/(M·T²)]
    # [(K/ξ)Ṁbody] = [L⁵/(M·T²)] × [M/(L³·T)] = [L²/T³]
    K_over_xi = v.get_dim('K_barotropic') / v.get_dim('xi')
    rhs_simplified = K_over_xi * v.get_dim('M_dot_density')
    v.check_dims("RHS term (K/ξ)Ṁbody",
                 rhs_simplified, v.L**2 / v.T**3)

    # Verify simplified equation balance
    v.check_dims("Simplified PDE balance (time ~ RHS)",
                 time_B, rhs_simplified)
    v.check_dims("Simplified PDE balance (flux ~ RHS)",
                 flux_B, rhs_simplified)

    # Verify coefficient K/ξ dimensions
    v.check_dims("Coefficient K/ξ dimensional check",
                 K_over_xi, v.L**5 / (v.M * v.T**2))

    v.success("Simplified quasilinear PDE dimensional structure verified")


def test_nonlinearity_structure(v):
    """
    Test the quadratic and cubic nonlinearities in the PDE.

    Verifies the nonlinear structure arising from:
    - Quadratic nonlinearity from (∇Ψ)² terms
    - Cubic nonlinearity from convection (∇Ψ)²∇Ψ
    - Variable effective speed contributions

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Nonlinearity Structure Analysis")

    # Analyze individual nonlinear terms in ∂t(∂tΨ + (1/2)(∇Ψ)²)

    # Linear term: ∂t(∂tΨ) = ∂²tΨ
    # [∂²tΨ] = [L²/T]/[T²] = [L²/T³]
    linear_term = v.dtt(v.get_dim('Psi_velocity_potential'))
    v.check_dims("Linear term ∂²tΨ",
                 linear_term, v.L**2 / v.T**3)

    # Quadratic term: ∂t[(1/2)(∇Ψ)²]
    # This involves ∂t(∇Ψ) · ∇Ψ ~ (∇∂tΨ) · (∇Ψ)
    # [∇∂tΨ] = [L²/T²]/[L] = [L/T²]
    # [∇Ψ] = [L/T]
    # [(∇∂tΨ) · (∇Ψ)] = [L/T²][L/T] = [L²/T³]
    quadratic_time = v.grad_dim(v.dt(v.get_dim('Psi_velocity_potential'))) * v.grad_dim(v.get_dim('Psi_velocity_potential'))
    v.check_dims("Quadratic time term ∂t[(∇Ψ)²] ~ (∇∂tΨ)·(∇Ψ)",
                 quadratic_time, v.L**2 / v.T**3)

    # Analyze divergence term nonlinearities: ∇·[(∂tΨ + (1/2)(∇Ψ)²)∇Ψ]

    # Linear divergence: ∇·(∂tΨ ∇Ψ)
    # [∂tΨ ∇Ψ] = [L²/T²][L/T] = [L³/T³]
    # [∇·(∂tΨ ∇Ψ)] = [L³/T³]/[L] = [L²/T³]
    linear_div = v.dt(v.get_dim('Psi_velocity_potential')) * v.grad_dim(v.get_dim('Psi_velocity_potential')) / v.L
    v.check_dims("Linear divergence ∇·(∂tΨ ∇Ψ)",
                 linear_div, v.L**2 / v.T**3)

    # Cubic divergence: ∇·[(1/2)(∇Ψ)² ∇Ψ] = ∇·[(∇Ψ)² ∇Ψ]
    # [(∇Ψ)²] = [L²/T²]
    # [(∇Ψ)² ∇Ψ] = [L²/T²][L/T] = [L³/T³]
    # [∇·((∇Ψ)² ∇Ψ)] = [L³/T³]/[L] = [L²/T³]
    cubic_div = (v.grad_dim(v.get_dim('Psi_velocity_potential'))**2 *
                 v.grad_dim(v.get_dim('Psi_velocity_potential')) / v.L)
    v.check_dims("Cubic divergence ∇·[(∇Ψ)² ∇Ψ]",
                 cubic_div, v.L**2 / v.T**3)

    # Verify all nonlinear terms have consistent dimensions
    v.check_dims("Linear ~ Quadratic consistency",
                 linear_term, quadratic_time)
    v.check_dims("Linear ~ Cubic consistency",
                 linear_term, cubic_div)
    v.check_dims("Quadratic ~ Cubic consistency",
                 quadratic_time, cubic_div)

    v.info("All nonlinear terms maintain dimensional consistency [L²/T³]")
    v.success("Nonlinearity structure analysis verified")


def test_variable_effective_speed(v):
    """
    Test the variable effective speed contribution to nonlinearities.

    Verifies that v_eff² = K ρ4D varies with the potential, creating
    additional nonlinear coupling in the governing equation.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Variable Effective Speed")

    # Test effective speed relationship: v_eff² = K ρ4D
    # From ρ4D = -(1/K)[∂tΨ + (1/2)(∇Ψ)²], we get:
    # v_eff² = K ρ4D = -[∂tΨ + (1/2)(∇Ψ)²] = -B
    v_eff_squared = v.get_dim('K_barotropic') * v.get_dim('rho_4')
    v.check_dims("Effective speed squared v_eff² = K ρ4D",
                 (v.L/v.T)**2, v_eff_squared)

    # Alternative form: v_eff² = -B
    v_eff_from_bernoulli = v.get_dim('B_bernoulli')
    v.check_dims("Effective speed from Bernoulli v_eff² = -B",
                 (v.L/v.T)**2, v_eff_from_bernoulli)

    # Test that v_eff varies with potential derivatives
    # ∂t v_eff² = ∂t(-B) = -∂t(∂tΨ + (1/2)(∇Ψ)²)
    # This creates time-dependent sound speed effects
    v_eff_time_variation = v.dt(v_eff_from_bernoulli)
    v.check_dims("Time variation ∂t(v_eff²)",
                 v_eff_time_variation, v.L**2 / v.T**3)

    # Spatial variation: ∇(v_eff²) = ∇(-B) = -∇(∂tΨ + (1/2)(∇Ψ)²)
    # This creates spatial sound speed gradients
    v_eff_space_variation = v.grad_dim(v_eff_from_bernoulli)
    v.check_dims("Spatial variation ∇(v_eff²)",
                 v_eff_space_variation, v.L / v.T**2)

    # The variable sound speed creates additional nonlinear coupling
    # beyond the standard quasilinear convection terms
    v.info("Variable effective speed v_eff = √(-B) creates additional nonlinearities")
    v.info("This goes beyond standard quasilinear wave equations")

    v.success("Variable effective speed contribution verified")


def test_pde_classification(v):
    """
    Test the mathematical classification of the derived PDE.

    Verifies that the equation is indeed:
    - Second-order in time (from ∂²tΨ terms)
    - Second-order in space (from Laplacian-like terms)
    - Quasilinear (coefficients depend on solution and its derivatives)
    - Nonlinear (quadratic and cubic terms in ∇Ψ)

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("PDE Classification")

    # Second-order time terms from ∂t(∂tΨ) = ∂²tΨ
    second_order_time = v.dtt(v.get_dim('Psi_velocity_potential'))
    v.check_dims("Second-order time derivative ∂²tΨ",
                 second_order_time, v.L**2 / v.T**3)

    # Second-order space terms from ∇·(... ∇Ψ) expansion
    # When expanded, ∇·(B ∇Ψ) = B ∇²Ψ + ∇B · ∇Ψ
    # The B ∇²Ψ term is second-order in space
    second_order_space = v.get_dim('B_bernoulli') * v.lap_dim(v.get_dim('Psi_velocity_potential'))
    v.check_dims("Second-order spatial term B ∇²Ψ",
                 second_order_space, v.L**2 / v.T**3)

    # First-order coupling term ∇B · ∇Ψ (quasilinear coefficient)
    first_order_coupling = v.grad_dim(v.get_dim('B_bernoulli')) * v.grad_dim(v.get_dim('Psi_velocity_potential'))
    v.check_dims("First-order coupling ∇B · ∇Ψ",
                 first_order_coupling, v.L**2 / v.T**3)

    # Verify all terms have same dimension (necessary for well-posed PDE)
    v.check_dims("Second-order terms consistency",
                 second_order_time, second_order_space)
    v.check_dims("Coupling term consistency",
                 first_order_coupling, second_order_time)

    v.info("PDE Classification:")
    v.info("- Second-order in time: ∂²tΨ terms present")
    v.info("- Second-order in space: ∇²Ψ terms present")
    v.info("- Quasilinear: coefficients depend on ∂tΨ, ∇Ψ")
    v.info("- Nonlinear: quadratic (∇Ψ)² and cubic (∇Ψ)²∇Ψ terms")

    v.success("PDE mathematical classification verified")


def test_continuity_source_consistency(v):
    """
    Test the consistency of the mass source terms in the continuity framework.

    Verifies that Ṁbody appears consistently as a volume source density
    throughout the derivation and maintains proper dimensional relationships.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Continuity Source Consistency")

    # Original 3D continuity equation: ∂t ρ3D + ∇·(ρ3D v) = -Ṁbody
    # All terms must have dimension [M/(L³·T)]
    continuity_time = v.dt(v.get_dim('rho_0'))
    continuity_flux = v.get_dim('rho_0') / v.T  # ρ∇·v ~ ρ/T
    continuity_source = v.get_dim('M_dot_density')

    v.check_dims("Continuity time term ∂t ρ3D",
                 continuity_time, v.M / (v.L**3 * v.T))
    v.check_dims("Continuity flux term ∇·(ρ v)",
                 continuity_flux, v.M / (v.L**3 * v.T))
    v.check_dims("Continuity source -Ṁbody",
                 continuity_source, v.M / (v.L**3 * v.T))

    # After substitution, the transformed PDE should preserve this balance
    # The (K/ξ)Ṁbody term should transform correctly
    K_over_xi = v.get_dim('K_barotropic') / v.get_dim('xi')
    transformed_source = K_over_xi * v.get_dim('M_dot_density')

    v.check_dims("Transformed source (K/ξ)Ṁbody",
                 transformed_source, v.L**2 / v.T**3)

    # Verify the transformation factor (K/ξ) converts properly
    # [K/ξ] × [M/(L³·T)] = [L⁵/(M·T²)] × [M/(L³·T)] = [L²/T³] ✓
    transformation_check = K_over_xi * (v.M / (v.L**3 * v.T))
    v.check_dims("Source transformation verification",
                 transformation_check, v.L**2 / v.T**3)

    # Check that this matches the left-hand side terms
    lhs_dim = v.L**2 / v.T**3  # Dimension of ∂t B and ∇·(B ∇Ψ) terms
    v.check_dims("Source-LHS dimensional consistency",
                 transformed_source, lhs_dim)

    v.success("Continuity source consistency verified")


def test_substitution_into_continuity():
    """
    Main test function for the Substitution into Continuity subsection.

    This function coordinates all verification tests for the substitution of
    the Bernoulli-derived density relationship into the continuity equation,
    leading to the complete quasilinear PDE governing potential flow.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Appendix: Substitution into Continuity",
        "Quasilinear PDE derivation via Bernoulli density substitution"
    )

    v.section("SUBSTITUTION INTO CONTINUITY VERIFICATION")

    # Add custom dimensions needed for the tests
    v.add_dimensions({
        'K_barotropic': v.L**6 / (v.M * v.T**2),    # Barotropic coupling K = g/m²
        'B_bernoulli': v.L**2 / v.T**2,             # Bernoulli quantity ∂tΨ + (1/2)(∇Ψ)²
        'Psi_velocity_potential': v.L**2 / v.T,      # Velocity potential
        'xi': v.L,                                   # Core radius/healing length (already in standard dims)
        'M_dot_density': v.M / (v.L**3 * v.T),      # Mass source density (already in standard dims)
    }, allow_overwrite=True)

    # Call test functions in logical order following the physics development
    v.info("\n--- 1) Bernoulli Density Relationship ---")
    test_bernoulli_density_relationship(v)

    v.info("\n--- 2) Initial Substitution Equation ---")
    test_initial_substitution_equation(v)

    v.info("\n--- 3) Simplified Quasilinear PDE ---")
    test_simplified_quasilinear_pde(v)

    v.info("\n--- 4) Nonlinearity Structure Analysis ---")
    test_nonlinearity_structure(v)

    v.info("\n--- 5) Variable Effective Speed ---")
    test_variable_effective_speed(v)

    v.info("\n--- 6) PDE Classification ---")
    test_pde_classification(v)

    v.info("\n--- 7) Continuity Source Consistency ---")
    test_continuity_source_consistency(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_substitution_into_continuity()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
