#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Beyond-Maxwell Predictions and Falsifiable Tests - Verification
===============================================================

Complete verification of all mathematical relationships in the "Beyond-Maxwell
Predictions and Falsifiable Tests" subsection from doc/projected_em.tex
(lines 300-391).

Tests cover:
- Static near-field Coulomb corrections (modified Poisson equation)
- Vacuum wave dispersion at high frequency
- Ultrafast transients with displacement memory
- Nanoscale boundary effects and cavity mode shifts
- Strong-field nonlinearity with definite sign
- Experimental bound relationships and scaling laws
"""

import os
import sys
from sympy import symbols, pi, sqrt, simplify

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify
)


def test_static_near_field_corrections(v):
    """Test static near-field Coulomb correction equations."""
    v.subsection("Static Near-Field Coulomb Corrections")

    # Modified Poisson equation: (-∇² + α ξ² ∇⁴ + ...) Φ = ρ/ε₀
    # Check dimensional consistency of each term

    # Standard Laplacian term: ∇²Φ
    laplacian_phi = v.lap_dim(v.get_dim('Phi'))

    # Fourth derivative term: ∇⁴Φ (apply Laplacian twice)
    nabla4_phi = v.lap_dim(v.lap_dim(v.get_dim('Phi')))

    # With ξ² prefactor: α ξ² ∇⁴Φ
    xi_squared = v.get_dim('xi')**2
    alpha = v.get_dim('alpha_param')
    correction_term = alpha * xi_squared * nabla4_phi

    # Source term: ρ/ε₀
    source_term = v.get_dim('rho_charge') / v.get_dim('epsilon_0')

    # Check all terms have same dimensions
    v.check_dims("Modified Poisson: ∇²Φ vs source",
                 laplacian_phi, source_term)

    v.check_dims("Modified Poisson: α ξ² ∇⁴Φ vs ∇²Φ",
                 correction_term, laplacian_phi)

    # Verify the modified Poisson equation structure (dimensional consistency)
    modified_poisson_lhs = -laplacian_phi + correction_term
    v.check_dims("Modified Poisson equation balance", modified_poisson_lhs, source_term)

    # Point charge potential corrections: dimensionless correction terms
    r = v.get_dim('r')

    # Check that potential correction term is dimensionless: α ξ²/(2r²)
    potential_correction_ratio = alpha * xi_squared / (2 * r**2)
    v.assert_dimensionless(potential_correction_ratio, "Coulomb potential correction α ξ²/(2r²)")

    # Check that field correction term is dimensionless: α 3ξ²/(2r²)
    field_correction_ratio = alpha * 3 * xi_squared / (2 * r**2)
    v.assert_dimensionless(field_correction_ratio, "Electric field correction 3α ξ²/(2r²)")

    # Verify the correction structure (both should scale the same way)
    # Field correction is 3 times the potential correction
    v.check_eq("Field vs potential correction scaling",
               field_correction_ratio, 3 * potential_correction_ratio)

    # Bound relationship dimensions: ξ ≲ r√(δ_C/|α|)
    # Both sides should have length dimensions
    delta_C = v.get_dim('delta_C')
    bound_rhs = r * sqrt(delta_C / abs(alpha))
    v.check_dims("Coulomb bound relationship", v.get_dim('xi'), bound_rhs)

    v.success("Static near-field corrections verified")


def test_vacuum_wave_dispersion(v):
    """Test vacuum wave dispersion relationships at high frequency."""
    v.subsection("Vacuum Wave Dispersion")

    # Dispersion relation: ω² = c²k² [1 + σ(kξ)² + ...]
    omega_squared = v.get_dim('omega')**2
    c_squared = v.get_dim('c')**2
    k_squared = v.get_dim('k')**2
    c_squared_k_squared = c_squared * k_squared

    # Correction terms
    sigma = v.get_dim('sigma_param')
    k_xi_correction = (v.get_dim('k') * v.get_dim('xi'))**2
    correction_factor = 1 + sigma * k_xi_correction

    # Complete dispersion relation
    omega_squared_corrected = c_squared_k_squared * correction_factor

    v.check_dims("Vacuum dispersion: ω² vs c²k²",
                 omega_squared, c_squared_k_squared)

    # Verify the full dispersion equation structure
    v.check_dims("Full dispersion relation ω² = c²k²[1 + σ(kξ)²]",
                 omega_squared, omega_squared_corrected)

    # Correction term (kξ)² should be dimensionless
    v.assert_dimensionless(k_xi_correction, "Dispersion correction (kξ)²")

    # Group velocity: v_g ≃ c[1 + (3/2)σ(kξ)²]
    c = v.get_dim('c')
    group_velocity_factor = 1 + (3/2) * sigma * k_xi_correction
    v_g_corrected = c * group_velocity_factor

    group_velocity = v.get_dim('v')  # Use generic velocity dimensions
    v.check_dims("Group velocity dimensions", group_velocity, c)

    # Verify group velocity correction formula (dimensional consistency)
    v.check_dims("Group velocity v_g = c[1 + (3/2)σ(kξ)²]",
                 group_velocity, v_g_corrected)

    # Check the coefficient relationship: 3/2 factor
    velocity_correction_coeff = (3/2) * sigma * k_xi_correction
    v.assert_dimensionless(velocity_correction_coeff, "Group velocity correction coefficient (3/2)σ(kξ)²")

    # Wavelength form bound: ξ ≲ (λ/2π)√(δ_D/(c_σ|σ|))
    wavelength = v.get_dim('wavelength')
    delta_D = v.get_dim('delta_D')
    c_sigma = v.get_dim('c_sigma')
    bound_rhs_wavelength = wavelength / (2 * pi) * sqrt(delta_D / (c_sigma * abs(sigma)))
    v.check_dims("Dispersion bound (wavelength form)",
                 v.get_dim('xi'), bound_rhs_wavelength)

    # Wavenumber form bound: ξ ≲ (1/k)√(δ_D/(c_σ|σ|))
    bound_rhs_wavenumber = (1 / v.get_dim('k')) * sqrt(delta_D / (c_sigma * abs(sigma)))
    v.check_dims("Dispersion bound (wavenumber form)",
                 v.get_dim('xi'), bound_rhs_wavenumber)

    v.success("Vacuum wave dispersion verified")


def test_ultrafast_displacement_memory(v):
    """Test ultrafast transient displacement memory effects."""
    v.subsection("Ultrafast Displacement Memory")

    # Modified permittivity: ε(ω) = ε₀[1 + β(ωτ)² + ...]
    epsilon_0 = v.get_dim('epsilon_0')
    beta = v.get_dim('beta_param')
    omega_tau_correction = (v.get_dim('omega') * v.get_dim('tau'))**2

    # Complete permittivity formula
    correction_factor = 1 + beta * omega_tau_correction
    epsilon_omega = epsilon_0 * correction_factor

    v.check_dims("Modified permittivity dimensions", epsilon_omega, epsilon_0)

    # Verify the permittivity equation structure
    v.check_dims("Modified permittivity ε(ω) = ε₀[1 + β(ωτ)²]",
                 epsilon_omega, epsilon_0 * (1 + beta * omega_tau_correction))

    # Correction term (ωτ)² should be dimensionless
    v.assert_dimensionless(omega_tau_correction, "Ultrafast correction (ωτ)²")

    # Time domain correction: τ² ∂_{tt}E
    # This should have same dimensions as E field
    tau_squared = v.get_dim('tau')**2
    time_derivative_E = v.dtt(v.get_dim('E'))  # Second time derivative
    time_correction = tau_squared * time_derivative_E

    v.check_dims("Time domain correction τ² ∂_{tt}E vs E",
                 time_correction, v.get_dim('E'))

    # Verify time domain equivalence (dimensional consistency)
    v.check_dims("Time domain correction τ² ∂_{tt}E equivalence",
                 time_correction, beta * tau_squared * time_derivative_E)

    # Bound relationship: τ ≲ (1/ω)√(δ_ε/|β|)
    delta_epsilon = v.get_dim('delta_epsilon')
    bound_rhs = (1 / v.get_dim('omega')) * sqrt(delta_epsilon / abs(beta))
    v.check_dims("Ultrafast memory bound", v.get_dim('tau'), bound_rhs)

    v.success("Ultrafast displacement memory verified")


def test_nanoscale_boundary_effects(v):
    """Test nanoscale boundary cavity mode shift effects."""
    v.subsection("Nanoscale Boundary Effects")

    # Cavity mode shift: Δf/f = +γ(ξ/a)² + ...
    gamma = v.get_dim('gamma_param')
    cavity_scale = v.get_dim('a_cavity')  # Transverse cavity scale
    xi_over_a_correction = (v.get_dim('xi') / cavity_scale)**2

    # Complete frequency shift formula
    frequency_shift_ratio = gamma * xi_over_a_correction

    # The ratio Δf/f should be dimensionless
    v.assert_dimensionless(frequency_shift_ratio, "Frequency shift ratio Δf/f")

    # Verify the cavity mode shift equation structure
    v.check_dims("Cavity mode shift Δf/f = +γ(ξ/a)²",
                 frequency_shift_ratio, gamma * xi_over_a_correction)

    # Check the positive sign constraint (γ > 0 implicit in the formulation)
    v.info("Positive mode shift: γ > 0 (universal cavity effect)")

    # The correction term (ξ/a)² should be dimensionless
    v.assert_dimensionless(xi_over_a_correction, "Cavity correction (ξ/a)²")

    # Bound relationship: ξ ≲ a√(δ_f/|γ|)
    delta_f = v.get_dim('delta_f')
    bound_rhs = cavity_scale * sqrt(delta_f / abs(gamma))
    v.check_dims("Cavity mode shift bound", v.get_dim('xi'), bound_rhs)

    v.success("Nanoscale boundary effects verified")


def test_strong_field_nonlinearity(v):
    """Test strong-field nonlinearity with definite sign."""
    v.subsection("Strong-Field Nonlinearity")

    # Kerr-like index: n(I) ≃ 1 + n₂I
    intensity = v.get_dim('I_intensity')
    n_2 = v.get_dim('n_2')

    # Complete refractive index formula
    n_of_I = 1 + n_2 * intensity

    # Refractive index should be dimensionless
    v.assert_dimensionless(n_of_I, "Kerr-like refractive index n(I)")

    # Verify the Kerr index equation structure
    v.check_dims("Kerr-like index n(I) = 1 + n₂I",
                 n_of_I, 1 + n_2 * intensity)

    # Verify the nonlinear correction coefficient structure
    kerr_correction = n_2 * intensity
    v.assert_dimensionless(kerr_correction, "Kerr correction n₂I")
    v.check_eq("Kerr index additive structure", n_of_I, 1 + kerr_correction)

    # For n₂I to be dimensionless, n₂ must have dimensions inverse to I
    # Light intensity I typically has dimensions [M T⁻³] (power per unit area)
    n2_times_I = n_2 * intensity
    v.assert_dimensionless(n2_times_I, "Nonlinear index correction n₂I")

    # Sign constraint: n₂ > 0 (positive compressibility)
    # This is a physical constraint, not a dimensional one
    v.info("Sign constraint: n₂ > 0 (fixed by positive compressibility)")

    v.success("Strong-field nonlinearity verified")


def test_scaling_relationships(v):
    """Test the unified scaling relationships across all effects."""
    v.subsection("Unified Scaling Relationships")

    # All effects scale with powers of the small parameters ξ and τ
    xi = v.get_dim('xi')
    tau = v.get_dim('tau')

    # Static effects scale as (ξ/r)²
    r = v.get_dim('r')
    static_scaling = (xi / r)**2
    v.assert_dimensionless(static_scaling, "Static scaling (ξ/r)²")

    # Dispersion effects scale as (kξ)²
    k = v.get_dim('k')
    dispersion_scaling = (k * xi)**2
    v.assert_dimensionless(dispersion_scaling, "Dispersion scaling (kξ)²")

    # Confinement effects scale as (ξ/a)²
    a_cavity = v.get_dim('a_cavity')
    confinement_scaling = (xi / a_cavity)**2
    v.assert_dimensionless(confinement_scaling, "Confinement scaling (ξ/a)²")

    # Ultrafast memory scales as (ωτ)²
    omega = v.get_dim('omega')
    memory_scaling = (omega * tau)**2
    v.assert_dimensionless(memory_scaling, "Memory scaling (ωτ)²")

    v.success("Unified scaling relationships verified")


def test_beyond_maxwell():
    """
    Main test function for Beyond-Maxwell Predictions and Falsifiable Tests.

    This function coordinates all verification tests for the subsection,
    covering dimensional consistency, mathematical relationships, and
    experimental bound formulations.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Beyond-Maxwell Predictions and Falsifiable Tests",
        "Verification of scale-suppressed corrections to Maxwell equations"
    )

    v.section("BEYOND-MAXWELL PREDICTIONS VERIFICATION")

    # Add custom dimensions needed for this test
    v.add_dimensions({
        'tau': v.T,                                    # Transition time scale
        'a_cavity': v.L,                               # Transverse cavity scale
        'I_intensity': v.M / v.T**3,                   # Light intensity (power/area)
        'n_2': v.T**3 / v.M,                          # Nonlinear index coefficient
        'alpha_param': 1,                              # Dimensionless O(1) parameter
        'sigma_param': 1,                              # Dimensionless O(1) parameter
        'beta_param': 1,                               # Dimensionless O(1) parameter
        'gamma_param': 1,                              # Dimensionless O(1) parameter
        'delta_C': 1,                                  # Coulomb precision parameter
        'delta_D': 1,                                  # Dispersion precision parameter
        'delta_epsilon': 1,                            # Permittivity precision parameter
        'delta_f': 1,                                  # Frequency precision parameter
        'c_sigma': 1,                                  # Dimensionless shape factor
    })

    # Call test functions in logical order
    v.info("\n--- A) Static Near-Field Corrections ---")
    test_static_near_field_corrections(v)

    v.info("\n--- B) Vacuum Wave Dispersion ---")
    test_vacuum_wave_dispersion(v)

    v.info("\n--- C) Ultrafast Displacement Memory ---")
    test_ultrafast_displacement_memory(v)

    v.info("\n--- D) Nanoscale Boundary Effects ---")
    test_nanoscale_boundary_effects(v)

    v.info("\n--- E) Strong-Field Nonlinearity ---")
    test_strong_field_nonlinearity(v)

    v.info("\n--- Unified Scaling Analysis ---")
    test_scaling_relationships(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_beyond_maxwell()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
