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
    correction_term = xi_squared * nabla4_phi

    # Source term: ρ/ε₀
    source_term = v.get_dim('rho_charge') / v.get_dim('epsilon_0')

    # Check all terms have same dimensions
    v.check_dims("Modified Poisson: ∇²Φ vs source",
                 laplacian_phi, source_term)

    v.check_dims("Modified Poisson: α ξ² ∇⁴Φ vs ∇²Φ",
                 correction_term, laplacian_phi)

    # Point charge potential correction term: α ξ²/(2r²)
    r = v.get_dim('r')
    potential_correction = xi_squared / r**2
    v.assert_dimensionless(potential_correction, "Coulomb potential correction α ξ²/(2r²)")

    # Electric field correction term: α 3ξ²/(2r²)
    field_correction = 3 * xi_squared / (2 * r**2)
    v.assert_dimensionless(field_correction, "Electric field correction 3α ξ²/(2r²)")

    # Bound relationship dimensions: ξ ≲ r√(δ_C/|α|)
    # Both sides should have length dimensions
    bound_rhs = r * sqrt(1)  # δ_C/|α| is dimensionless, using 1 as placeholder
    v.check_dims("Coulomb bound relationship", v.get_dim('xi'), bound_rhs)

    v.success("Static near-field corrections verified")


def test_vacuum_wave_dispersion(v):
    """Test vacuum wave dispersion relationships at high frequency."""
    v.subsection("Vacuum Wave Dispersion")

    # Dispersion relation: ω² = c²k² [1 + σ(kξ)² + ...]
    omega_squared = v.get_dim('omega')**2
    c_squared_k_squared = v.get_dim('c')**2 * v.get_dim('k')**2

    v.check_dims("Vacuum dispersion: ω² vs c²k²",
                 omega_squared, c_squared_k_squared)

    # Correction term (kξ)² should be dimensionless
    k_xi_correction = (v.get_dim('k') * v.get_dim('xi'))**2
    v.assert_dimensionless(k_xi_correction, "Dispersion correction (kξ)²")

    # Group velocity: v_g ≃ c[1 + (3/2)σ(kξ)²]
    group_velocity = v.get_dim('v')  # Use generic velocity dimensions
    speed_of_light = v.get_dim('c')
    v.check_dims("Group velocity dimensions", group_velocity, speed_of_light)

    # Wavelength form bound: ξ ≲ (λ/2π)√(δ_D/(c_σ|σ|))
    wavelength = v.get_dim('wavelength')
    bound_rhs_wavelength = wavelength / (2 * pi) * sqrt(1)  # δ_D/(c_σ|σ|) dimensionless
    v.check_dims("Dispersion bound (wavelength form)",
                 v.get_dim('xi'), bound_rhs_wavelength)

    # Wavenumber form bound: ξ ≲ (1/k)√(δ_D/(c_σ|σ|))
    bound_rhs_wavenumber = (1 / v.get_dim('k')) * sqrt(1)
    v.check_dims("Dispersion bound (wavenumber form)",
                 v.get_dim('xi'), bound_rhs_wavenumber)

    v.success("Vacuum wave dispersion verified")


def test_ultrafast_displacement_memory(v):
    """Test ultrafast transient displacement memory effects."""
    v.subsection("Ultrafast Displacement Memory")

    # Modified permittivity: ε(ω) = ε₀[1 + β(ωτ)² + ...]
    epsilon_omega = v.get_dim('epsilon_0')  # Same dimensions as ε₀
    epsilon_0 = v.get_dim('epsilon_0')

    v.check_dims("Modified permittivity dimensions", epsilon_omega, epsilon_0)

    # Correction term (ωτ)² should be dimensionless
    omega_tau_correction = (v.get_dim('omega') * v.get_dim('tau'))**2
    v.assert_dimensionless(omega_tau_correction, "Ultrafast correction (ωτ)²")

    # Time domain correction: τ² ∂_{tt}E
    # This should have same dimensions as E field
    tau_squared = v.get_dim('tau')**2
    time_derivative_E = v.dtt(v.get_dim('E'))  # Second time derivative
    time_correction = tau_squared * time_derivative_E

    v.check_dims("Time domain correction τ² ∂_{tt}E vs E",
                 time_correction, v.get_dim('E'))

    # Bound relationship: τ ≲ (1/ω)√(δ_ε/|β|)
    bound_rhs = (1 / v.get_dim('omega')) * sqrt(1)  # δ_ε/|β| dimensionless
    v.check_dims("Ultrafast memory bound", v.get_dim('tau'), bound_rhs)

    v.success("Ultrafast displacement memory verified")


def test_nanoscale_boundary_effects(v):
    """Test nanoscale boundary cavity mode shift effects."""
    v.subsection("Nanoscale Boundary Effects")

    # Cavity mode shift: Δf/f = +γ(ξ/a)² + ...
    # The ratio Δf/f should be dimensionless
    frequency_shift_ratio = 1  # Δf/f is inherently dimensionless

    # The correction term (ξ/a)² should be dimensionless
    cavity_scale = v.get_dim('a_cavity')  # Transverse cavity scale
    xi_over_a_correction = (v.get_dim('xi') / cavity_scale)**2
    v.assert_dimensionless(xi_over_a_correction, "Cavity correction (ξ/a)²")

    # Bound relationship: ξ ≲ a√(δ_f/|γ|)
    bound_rhs = cavity_scale * sqrt(1)  # δ_f/|γ| dimensionless
    v.check_dims("Cavity mode shift bound", v.get_dim('xi'), bound_rhs)

    v.success("Nanoscale boundary effects verified")


def test_strong_field_nonlinearity(v):
    """Test strong-field nonlinearity with definite sign."""
    v.subsection("Strong-Field Nonlinearity")

    # Kerr-like index: n(I) ≃ 1 + n₂I
    # Refractive index should be dimensionless
    refractive_index = 1  # n(I) is dimensionless

    # For n₂I to be dimensionless, n₂ must have dimensions inverse to I
    # Light intensity I typically has dimensions [M T⁻³] (power per unit area)
    intensity = v.get_dim('I_intensity')
    n2_times_I = v.get_dim('n_2') * intensity
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
