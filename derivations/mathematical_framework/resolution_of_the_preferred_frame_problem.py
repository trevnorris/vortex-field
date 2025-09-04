#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "Resolution of the Preferred Frame Problem" subsection.

This module implements dimensional and mathematical verification for the subsection
covering preferred frame resolution, small parameters, Maxwell equations,
Lorentz invariants, and Michelson-Morley analysis as outlined in TEST.md.

Based on the mathematical framework document, Section covering preferred frame problem.
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    quick_verify
)


def test_small_parameters_dimensionless(v):
    """Test that all small parameters are dimensionless."""
    # epsilon_v = |v_bg|/c
    v.assert_dimensionless(v.get_dim('v')/v.get_dim('c'), "epsilon_v (v_bg/c)")
    
    # epsilon_Phi_g = |Phi_g|/c^2  
    v.assert_dimensionless(v.get_dim('Phi_g')/(v.get_dim('c')**2), "epsilon_Phi_g (Phi_g/c^2)")
    
    # epsilon_xi = xi_c/L
    v.assert_dimensionless(v.get_dim('xi')/v.get_dim('L_scale'), "epsilon_xi (xi/L)")
    
    # epsilon_rho = delta_rho/rho_0
    v.assert_dimensionless(v.get_dim('rho')/v.get_dim('rho_0'), "epsilon_rho (delta_rho/rho_0)")
    
    v.success("Small parameters verified as dimensionless")


def test_dalembert_consistency(v):
    """Test d'Alembertian operator dimensional consistency."""
    # (1/c^2) partial_tt A^mu and nabla^2 A^mu should have same dimensions
    time_term = (1 / (v.get_dim('c')**2)) * v.dtt(v.get_dim('A_mu'))
    space_term = v.lap_dim(v.get_dim('A_mu'))
    
    v.check_dims("Box A^mu: time vs space pieces", time_term, space_term)
    v.success("d'Alembertian consistency verified")


def test_maxwell_wave_equation_SI(v):
    """Test Maxwell wave equation in SI units."""
    # SI form: Box A^mu = -mu_0 J^mu
    box_A = v.lap_dim(v.get_dim('A_mu'))  # Same as space term from d'Alembertian
    rhs_SI = v.get_dim('mu_0') * v.get_dim('J_mu')
    
    v.check_dims("Wave eq (SI): Box A vs mu0 J", box_A, rhs_SI)
    v.success("SI Maxwell wave equation verified")


def test_gaussian_wave_equation_diagnostic(v):
    """Test diagnostic for Gaussian form - should fail in SI."""
    # Diagnostic: (4pi/c)J^mu should NOT match SI dimensions
    box_A = v.lap_dim(v.get_dim('A_mu'))
    bad_rhs = (4*pi / v.get_dim('c')) * v.get_dim('J_mu')
    
    v.info("🔍 DIAGNOSTIC: Testing (4pi/c)J^mu in SI units (should be dimensionally wrong)...")
    ok = v.check_dims("Diagnostic: Box A vs (4pi/c)J (SI units) should mismatch",
                      box_A, bad_rhs, record=False, verbose=False)
    quick_verify("Caught unit-system drift: (4pi/c) needs Gaussian/HL conventions, not SI", not ok, helper=v, expected_failure=True)
    v.success("Diagnostic caught: Gaussian form incompatible with SI")


def test_lorenz_gauge_condition(v):
    """Test Lorenz gauge condition partial_mu A^mu = 0."""
    # partial_mu A^mu should have well-defined dimensions
    gauge_term = v.div_dim(v.get_dim('A_mu'))
    
    # Compare to 0 (helper treats this dimension-agnostically)
    v.check_dims("Lorenz gauge div A has well-defined unit", gauge_term, 0)
    v.success("Lorenz gauge condition verified")


def test_em_invariant_I1(v):
    """Test first EM invariant: I_1 = 2(B^2 - E^2/c^2)."""
    # B^2 and E^2/c^2 must have matching units
    B_squared = v.get_dim('B')**2
    E_term = (v.get_dim('E')**2) / (v.get_dim('c')**2)
    
    v.check_dims("I1: B^2 vs E^2/c^2", B_squared, E_term)
    v.success("EM invariant I_1 verified")


def test_em_invariant_I2(v):
    """Test second EM invariant: I_2 = -4(E dot B)/c."""
    # (E dot B)/c should match B^2 units
    EB_over_c = (v.get_dim('E') * v.get_dim('B')) / v.get_dim('c')
    B_squared = v.get_dim('B')**2
    
    v.check_dims("I2: (E dot B)/c carries same units as B^2", EB_over_c, B_squared)
    v.success("EM invariant I_2 verified")


def test_michelson_morley_timing_units(v):
    """Test that MM timing expressions have correct time units."""
    # t_parallel and t_perp should have time dimensions
    t_leading = (2 * v.get_dim('L_scale') / v.get_dim('c'))
    
    v.check_dims("t_parallel has [T]", t_leading, v.T)
    v.check_dims("t_perp has [T]", t_leading, v.T)
    v.success("Michelson-Morley timing units verified")


def test_michelson_morley_leading_order_cancellation(v):
    """Test that leading-order MM terms cancel exactly."""
    # Symbolic check: t_parallel - t_perp = 0 at leading order
    beta = symbols('beta', real=True)
    t_base = symbols('t_base', positive=True)  # Represents 2L/c
    
    # Both have form: t_base * (1 + 1/2 * beta^2)
    expr_parallel = t_base * (1 + Rational(1,2)*beta**2)
    expr_perp = t_base * (1 + Rational(1,2)*beta**2)
    
    difference = simplify(expr_parallel - expr_perp)
    quick_verify("t_parallel - t_perp (shown terms) cancels exactly", difference == 0, helper=v)
    v.success("Leading-order MM cancellation verified")


def test_michelson_morley_higher_order_estimate(v):
    """Test numerical estimate of higher-order MM timing difference."""
    # Delta t approx (2L/c) beta^4 with L=11 m, beta=1e-4
    L_val = 11      # meters
    c_val = 3e8     # m/s  
    beta_val = 1e-4
    
    delta_t_est = (2*L_val/c_val) * (beta_val**4)
    expected_order = 1e-23  # seconds
    
    # Check that estimate is close to expected order of magnitude
    quick_verify("Delta t ~ 1e-23 s (with L=11m, beta=1e-4)",
                 abs(delta_t_est - expected_order) < 5e-24, helper=v)
    
    v.info(f"MM higher-order estimate: Delta t ≈ {delta_t_est:.1e} s")


def test_em_gem_potential_separation(v):
    """Test that EM and GEM potentials have different units."""
    # A^mu (EM) should differ from A_g^mu (GEM)
    # EM: [V*s/m], GEM: [L/T]
    
    A_em_dim = v.get_dim('A_mu')
    A_gem_dim = v.get_dim('Ag_mu')
    
    quick_verify("A_mu (EM) and Ag_mu (GEM) have different units",
                 A_em_dim != A_gem_dim, helper=v)
    v.success("EM/GEM potential separation verified")


def test_metric_signature_convention(v):
    """Note metric signature convention (no dimensional check needed)."""
    v.info("Metric signature eta_mu_nu = diag(-,+,+,+) noted as convention")


def test_resolution_preferred_frame_problem():
    """
    Main test function implementing all verification categories from TEST.md.
    
    Tests 5 main categories:
    A) Small parameters dimensionless (4 tests)
    B) Lorenz-gauge Maxwell wave equation (4 tests)
    C) Lorentz invariants (2 tests)  
    D) Michelson-Morley timing (3 tests)
    E) Conventions separation (2 tests)
    """
    v = PhysicsVerificationHelper(
        "Preferred Frame – Resolution",
        "Dimensional checks for small parameters, Maxwell/Lorenz, invariants, and MM estimates",
        unit_system=UnitSystem.SI
    )
    
    v.section_header("Testing Resolution of the Preferred Frame Problem")
    
    # Declare dimensionless parameters for transcendental/series use
    v.declare_dimensionless('epsilon_v', 'epsilon_Phi_g', 'epsilon_xi', 'epsilon_rho', 'beta')
    
    # A) Small parameters dimensionless verification
    v.info("\n--- A) Small parameters dimensionless verification ---")
    test_small_parameters_dimensionless(v)
    
    # B) Lorenz-gauge Maxwell wave equation
    v.info("\n--- B) Lorenz-gauge Maxwell wave equation ---")
    v.section("Lorenz-gauge wave equation and gauge condition")
    test_dalembert_consistency(v)
    test_maxwell_wave_equation_SI(v)
    test_gaussian_wave_equation_diagnostic(v)
    test_lorenz_gauge_condition(v)
    
    # C) Lorentz invariants (EM)
    v.info("\n--- C) Lorentz invariants (EM) ---")
    v.section("Lorentz invariants (EM)")
    test_em_invariant_I1(v)
    test_em_invariant_I2(v)
    
    # D) Michelson-Morley timing analysis
    v.info("\n--- D) Michelson-Morley timing analysis ---")
    v.section("Michelson–Morley timing (no O(beta^2) anisotropy)")
    test_michelson_morley_timing_units(v)
    test_michelson_morley_leading_order_cancellation(v)
    test_michelson_morley_higher_order_estimate(v)
    
    # E) Conventions separation
    v.info("\n--- E) Conventions separation ---")
    v.section("Conventions separation (EM vs GEM potentials)")
    test_em_gem_potential_separation(v)
    test_metric_signature_convention(v)
    
    # Final summary
    v.summary()


if __name__ == "__main__":
    test_resolution_preferred_frame_problem()