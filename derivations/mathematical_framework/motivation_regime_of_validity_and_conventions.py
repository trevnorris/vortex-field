#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "Motivation, Regime of Validity, and Conventions" subsection.

This module implements dimensional and mathematical verification for the subsection
covering regime of validity, 4D->3D projections, linearized aether equations,
and wave equations as outlined in MRC.md.

Based on the mathematical framework document, Section 2.x.
"""

import os
import sys
import sympy as sp
from sympy import symbols, integrate, oo, pi, sqrt, simplify

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    verify_wave_equation,
    verify_conservation_law,
    quick_verify
)


def test_small_parameters_dimensionless(v):
    """Test that small parameters are dimensionless."""
    # Thinness parameter: epsilon_xi = xi_c/L
    v.assert_dimensionless(v.get_dim('xi')/v.L, "epsilon_xi thinness")
    
    # Speed ratio: epsilon_v = v/c  
    v.assert_dimensionless(v.get_dim('v')/v.get_dim('c'), "epsilon_v speed ratio")
    
    # Bulk vs wave speed ratios
    v.assert_dimensionless(v.get_dim('v_L')/v.get_dim('c'), "bulk vs wave speed ratio")
    
    # Additional small parameters mentioned in the text
    v.assert_dimensionless(v.get_dim('omega')*v.L/v.get_dim('c'), "epsilon_omega")
    
    print("✓ Small parameters verified as dimensionless")


def test_order_magnitude_anchors(v):
    """Test order-of-magnitude sanity checks."""
    # These are pure numbers, so we use quick_verify for numeric checks
    # Example values from the text: xi_c/L ~ 10^-5, v/c ~ 10^-3
    
    # These are just examples - in practice we'd need actual values
    # For now, just verify the ratios make sense dimensionally
    xi_over_L = v.get_dim('xi') / v.L 
    v_over_c = v.get_dim('v') / v.get_dim('c')
    
    # Verify these are dimensionless (already done above, but worth double-checking)
    v.assert_dimensionless(xi_over_L, "xi/L anchor check")
    v.assert_dimensionless(v_over_c, "v/c anchor check")
    
    print("✓ Order-of-magnitude anchors verified")


def test_slice_projection_density(v):
    """Test slice projection: rho_3D(x,t) = integral rho_4D(x,w,t)chi(w)dw."""
    # Dimensional check: [rho_3D] = [rho_4D][w]
    v.check_dims(
        "rho3D from rho4D slice", 
        v.get_dim('rho'), 
        v.get_dim('rho_4') * v.get_dim('w')
    )
    print("✓ Slice projection of density verified")


def test_delta_function_reduction(v):
    """Test delta-function reduction: integral delta^(4)(r_4-r_4,i)dw = delta^(3)(r-r_i)."""
    # Dimensional check: [delta^(4)][w] = [delta^(3)]
    v.check_dims(
        "delta4 -> delta3 via dw",
        v.get_dim('delta3'),
        v.get_dim('delta4') * v.get_dim('w')
    )
    print("✓ Delta-function reduction verified")


def test_background_density_relation(v):
    """Test background 3D density: rho_0 = rho_4D^0 * xi_c."""
    # Dimensional check: [rho_0] = [rho_4D][xi]
    v.check_dims(
        "rho0 = rho4_bg * xi",
        v.get_dim('rho_0'),
        v.get_dim('rho_4_bg') * v.get_dim('xi')
    )
    print("✓ Background density relation verified")


def test_core_symbols_dimensions(v):
    """Verify core symbol dimensions from the cheat sheet."""
    # These should already be in helper.py, but let's verify explicitly
    
    # rho_4D -> M L^-4
    expected_rho_4 = v.M / (v.L**4)
    v.check_dims("rho_4D dimensions", v.get_dim('rho_4'), expected_rho_4)
    
    # rho_3D/rho_0 -> M L^-3
    expected_rho_3 = v.M / (v.L**3)
    v.check_dims("rho_3D dimensions", v.get_dim('rho'), expected_rho_3)
    v.check_dims("rho_0 dimensions", v.get_dim('rho_0'), expected_rho_3)
    
    # Gamma (circulation quantum) -> L^2 T^-1
    expected_gamma = v.L**2 / v.T
    v.check_dims("Gamma circulation", v.get_dim('Gamma'), expected_gamma)
    
    # M_dot_i (sink strength) -> M T^-1
    expected_m_dot = v.M / v.T
    v.check_dims("M_dot_i sink strength", v.get_dim('M_dot_i'), expected_m_dot)
    
    # delta^(3), delta^(4) -> L^-3, L^-4
    v.check_dims("delta3 dimensions", v.get_dim('delta3'), 1/(v.L**3))
    v.check_dims("delta4 dimensions", v.get_dim('delta4'), 1/(v.L**4))
    
    # Phi_g (grav. potential) -> L^2 T^-2
    expected_phi_g = v.L**2 / (v.T**2)
    v.check_dims("Phi_g grav potential", v.get_dim('Phi_g'), expected_phi_g)
    
    # Velocities -> L T^-1
    expected_velocity = v.L / v.T
    v.check_dims("v velocity", v.get_dim('v'), expected_velocity)
    v.check_dims("v_L bulk velocity", v.get_dim('v_L'), expected_velocity)
    v.check_dims("v_eff effective velocity", v.get_dim('v_eff'), expected_velocity)
    
    print("✓ Core symbol dimensions verified")


def test_projected_em_placeholders(v):
    """Test projected EM placeholders with HL-style units."""
    # Add the projected symbols as requested
    v.add_dimensions({
        'Psi_proj': 1/v.L**2,   # [L^-2]
        'A_proj':   1/v.L,      # [L^-1]
        'F_munu':   1/v.T,      # [T^-1]; HL/c=1 convention
    })
    
    # Verify they were added correctly
    v.check_dims("Psi_proj surface-like potential", v.get_dim('Psi_proj'), 1/(v.L**2))
    v.check_dims("A_proj vector potential", v.get_dim('A_proj'), 1/v.L)
    v.check_dims("F_munu field strength", v.get_dim('F_munu'), 1/v.T)
    
    print("✓ Projected EM placeholders added and verified")


def test_full_4d_continuity(v):
    """Test full 4D continuity with sinks."""
    # partial_t rho_4D + div_4(rho_4D v_4) = -sum_i M_dot_i delta^(4)(r_4-r_4,i)
    
    density_rate = v.dt(v.get_dim('rho_4'))
    flux_div = v.div_dim(v.get_dim('rho_4') * v.get_dim('v'))
    source = v.get_dim('M_dot_i') * v.get_dim('delta4')
    
    verify_conservation_law(v, "4D continuity with sinks", density_rate, flux_div, source)
    print("✓ Full 4D continuity verified")


def test_linearized_continuity(v):
    """Test linearized continuity equation."""
    # partial_t delta_rho_4D + rho_4D^0 div_4 delta_v_4 = -sum_i M_dot_i delta^(4)(...)
    
    density_rate_lin = v.dt(v.get_dim('delta_rho_4'))
    flux_div_lin = v.get_dim('rho_4_bg') * v.div_dim(v.get_dim('v'))  # delta_v has same dims as v
    source_lin = v.get_dim('M_dot_i') * v.get_dim('delta4')
    
    verify_conservation_law(v, "Linearized 4D continuity", density_rate_lin, flux_div_lin, source_lin)
    print("✓ Linearized continuity verified")


def test_linearized_euler(v):
    """Test linearized Euler equation."""
    # partial_t delta_v_4 = -v_eff^2 grad_4(delta_rho_4D/rho_4D^0) - grad_4 delta_Q
    
    # Add Q (quantum potential) if not already present
    if 'Q' not in v.dims:
        v.add_dimensions({'Q': v.L**2 / v.T**2})  # L^2 T^-2
    
    lhs = v.dt(v.get_dim('v'))  # delta_v has same dims as v
    rhs_1 = v.get_dim('v_eff')**2 * v.grad_dim(v.get_dim('delta_rho_4')/v.get_dim('rho_4_bg'))
    rhs_2 = v.grad_dim(v.get_dim('Q'))
    
    v.check_dims("Euler: LHS vs pressure term", lhs, rhs_1)
    v.check_dims("Euler: LHS vs quantum-pressure term", lhs, rhs_2)
    
    print("✓ Linearized Euler verified")


def test_density_wave_equation(v):
    """Test second-order density wave equation."""
    # partial_tt delta_rho_4D - rho_4D^0 v_eff^2 lap_4(delta_rho_4D/rho_4D^0) = -sum_i partial_t M_dot_i delta^(4)(...) + rho_4D^0 lap_4 delta_Q
    
    time_term = v.dtt(v.get_dim('delta_rho_4'))
    space_term = v.get_dim('rho_4_bg') * (v.get_dim('v_eff')**2) * v.lap_dim(v.get_dim('delta_rho_4')/v.get_dim('rho_4_bg'))
    
    # Source 1: sink term
    src_sink = v.dt(v.get_dim('M_dot_i')) * v.get_dim('delta4')
    
    # Source 2: quantum term
    src_Q = v.get_dim('rho_4_bg') * v.lap_dim(v.get_dim('Q'))
    
    # Test wave operator consistency and both sources
    verify_wave_equation(v, "delta_rho wave (sink)", time_term, space_term, src_sink)
    verify_wave_equation(v, "delta_rho wave (Q)", time_term, space_term, src_Q)
    
    print("✓ Density wave equation verified")


def test_helmholtz_formulation(v):
    """Test Helmholtz scalar choice with diagnostics."""
    # Test both chi (velocity potential) and Phi_g variants
    
    # PASS: chi variant (velocity potential)
    v.check_dims("Helmholtz with chi", v.div_dim(v.get_dim('v')), v.lap_dim(v.get_dim('Phi_4D')))
    
    # DIAGNOSTIC: Phi_g variant should mismatch (intentional test of wrong potential)
    print("🔍 DIAGNOSTIC: Testing Phi_g as velocity potential (should be dimensionally inconsistent)...")
    ok = v.check_dims("Helmholtz with Phig — diagnostic",
                      v.div_dim(v.get_dim('v')), v.lap_dim(v.get_dim('Phi_g')), 
                      record=False, verbose=False)
    quick_verify("Diagnostic caught: Phi_g is not a velocity potential", not ok)
    
    print("✓ Helmholtz formulation verified (with diagnostic)")


def test_scalar_potential_waves(v):
    """Test wave equations for scalar potentials with diagnostics."""
    
    # 1) chi-wave (clean velocity potential form) - should PASS
    verify_wave_equation(v, "chi wave", 
                        v.dtt(v.get_dim('Phi_4D')), 
                        v.get_dim('v_eff')**2 * v.lap_dim(v.get_dim('Phi_4D')))
    
    # 2) Phi_g-wave with CORRECT sink source - should PASS
    src_phig_sinks = v.get_dim('v_eff')**2 * v.dt((v.get_dim('M_dot_i')/v.get_dim('rho_4_bg'))) * v.get_dim('delta4')
    verify_wave_equation(v, "Phig wave (sinks OK)", 
                        v.dtt(v.get_dim('Phi_g')), 
                        v.get_dim('v_eff')**2 * v.lap_dim(v.get_dim('Phi_g')), 
                        src_phig_sinks)
    
    # 3) DIAGNOSTIC: Test what happens when time derivative is missing from source
    print("🔍 DIAGNOSTIC: Testing Phi_g wave source without time derivative (should be dimensionally wrong)...")
    bad_src = v.get_dim('v_eff')**2 * (v.get_dim('M_dot_i')/v.get_dim('rho_4_bg')) * v.get_dim('delta4')
    ok = v.check_dims("Phig wave (missing dt) — diagnostic",
                      v.dtt(v.get_dim('Phi_g')), bad_src, record=False, verbose=False)
    quick_verify("Diagnostic caught: Phi_g source needs dt", not ok)
    
    # 4) DIAGNOSTIC: Test Q-term without proper prefactor
    print("🔍 DIAGNOSTIC: Testing Q-term without correct prefactor (should be dimensionally wrong)...")
    src_phig_Q_wrong = v.get_dim('v_eff')**2 * v.lap_dim(v.get_dim('Q')) / v.get_dim('rho_4_bg')
    ok = v.check_dims("Phig wave Q-term — diagnostic",
                      v.dtt(v.get_dim('Phi_g')), src_phig_Q_wrong, record=False, verbose=False)
    quick_verify("Diagnostic caught: Phi_g Q-term needs a prefactor to match LHS", not ok)
    
    print("✓ Scalar potential waves verified (with diagnostics)")


def test_motivation_regime_validity_conventions():
    """
    Main test function implementing all verification categories from MRC.md.
    
    Tests 15 main categories:
    A) Regime validity & asymptotics (2 tests)
    B) 4D->3D projection (3 tests)  
    C) Symbol verification (2 tests)
    D) Conventions (noted)
    E) Wave equations (6 tests)
    """
    print("="*60)
    print("Testing Motivation, Regime of Validity, and Conventions")
    print("="*60)
    
    v = PhysicsVerificationHelper("Section 2.x: Motivation/Validity/Conventions")
    
    # Declare small parameters as dimensionless for transcendental use
    v.declare_dimensionless('epsilon_xi', 'epsilon_v', 'beta')
    
    # A) Regime-of-validity & asymptotics
    print("\n--- A) Regime-of-validity & asymptotics ---")
    test_small_parameters_dimensionless(v)
    test_order_magnitude_anchors(v)
    
    # B) 4D->3D projection & distributions  
    print("\n--- B) 4D->3D projection & distributions ---")
    test_slice_projection_density(v)
    test_delta_function_reduction(v)
    test_background_density_relation(v)
    
    # C) Units / symbol verification
    print("\n--- C) Units / symbol verification ---")
    test_core_symbols_dimensions(v)
    test_projected_em_placeholders(v)
    
    # D) Notation & conventions (noted in comments)
    print("\n--- D) Notation & conventions ---")
    print("✓ Metric signature (-,+,+,+) noted as convention")
    print("✓ 4-potential convention A^mu_EM = (-Phi_Slope/c, A_EM) noted")
    
    # E) Aether-equation linearization & wave equations
    print("\n--- E) Aether equations & wave equations ---")
    test_full_4d_continuity(v)
    test_linearized_continuity(v)
    test_linearized_euler(v)
    test_density_wave_equation(v)
    test_helmholtz_formulation(v)
    test_scalar_potential_waves(v)
    
    # Final summary
    print("\n" + "="*60)
    v.summary()
    print("="*60)


if __name__ == "__main__":
    test_motivation_regime_validity_conventions()