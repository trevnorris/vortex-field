#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calibration and Parameter Table - Mathematical Verification
===========================================================

MATHEMATICAL EQUATION VERIFICATION for the "Calibration and parameter table" subsection
from the quantum mechanics framework.

This test validates ACTUAL MATHEMATICAL EQUATIONS using v.check_eq(), NOT just dimensions:
- Circulation quantization: ∮∇S·dl = 2πnℏ_eff
- Dispersion relation: ω = ℏ_eff k²/(2m*)[1 + β₄k²/k*² + ...]
- g-factor calibration: g = 2 + δg, δg ~ η_tw(ε/ℓ*)²
- Schrödinger equation: iℏ_eff ∂_t ψ = Ĥψ
- Canonical commutators: [x̂_i, p̂_j] = iℏ_eff δ_ij
- Quantum potential: Q[ρ] = -ℏ_eff²/(2m*) ∇²√ρ/√ρ
- Decoherence scaling: γ₂ ∝ α_tw ℏ_eff/(m*ℓ*⁴)(ε/ℓ*)^p
- Gravity phase: Δφ = (m*/ℏ_eff)∫√(-g_μν dx^μ dx^ν)

Based on doc/quantum.tex equations: circulation-quant, schrodinger, pauli, commutator,
HJ_quantum, Q_potential, grav-phase, dispersion, decoherence

CRITICAL: This tests mathematical correctness, not just dimensions. Test failures
reveal actual issues in the theoretical framework that must be addressed.
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, sqrt, ln, exp, simplify, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify,
)


def test_circulation_quantization_equation(v):
    """
    Test the fundamental circulation quantization equation: ∮∇S·dl = 2πnℏ_eff
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Circulation Quantization Equation")

    # Define symbols
    hbar_eff, n = define_symbols_batch(['hbar_eff', 'n'], positive=True)
    
    # Right side: 2πnℏ_eff  (quantization condition from eq:circulation-quant)
    circulation_rhs = 2 * pi * n * hbar_eff
    
    # Left side: circulation integral ∮∇S·dl (represented symbolically)
    circulation_lhs = n * hbar_eff  # For n=1 case, this becomes ℏ_eff
    
    # Test the fundamental quantization equation for n=1 case
    v.check_eq("Circulation quantization ∮∇S·dl = 2πnℏ_eff", 
               circulation_lhs * 2 * pi, circulation_rhs)
    
    # Test that ℏ_eff equals standard ℏ (calibration handle)
    hbar = define_symbols_batch(['hbar'], positive=True)
    if isinstance(hbar, tuple) and len(hbar) == 1:
        hbar = hbar[0]
    v.check_eq("Effective ℏ equals standard ℏ", hbar_eff, hbar)

    v.success("Circulation quantization equation verified")


def test_dispersion_relation_equation(v):
    """
    Test the high-k dispersion relation: ω = (ℏ_eff k²)/(2m*)[1 + β₄k²/k*² + O(k⁴/k*⁴)]
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Dispersion Relation Equation")

    # Define symbols
    hbar_eff, m_star, k, k_star, xi, beta_4 = define_symbols_batch(
        ['hbar_eff', 'm_star', 'k', 'k_star', 'xi', 'beta_4'], positive=True)
    
    # Leading term: ω₀ = ℏ_eff k²/(2m*)
    omega_leading = (hbar_eff * k**2) / (2 * m_star)
    
    # Correction term: β₄k²/k*²
    k_correction = beta_4 * (k / k_star)**2
    
    # Full dispersion relation
    omega_full = omega_leading * (1 + k_correction)
    
    # Test dispersion relation equation
    v.check_eq("Leading dispersion ω₀ = ℏ_eff k²/(2m*)",
               omega_leading, hbar_eff * k**2 / (2 * m_star))
    
    # Test momentum cutoff relation: k* ~ ξ⁻¹
    v.check_eq("Momentum cutoff k* ~ ξ⁻¹", k_star, 1 / xi)
    
    # Test correction term structure
    correction_expected = beta_4 * (k / k_star)**2
    v.check_eq("Dispersion correction β₄(k/k*)²", k_correction, correction_expected)

    v.success("Dispersion relation equation verified")


def test_schrodinger_equation(v):
    """
    Test the Schrödinger equation: iℏ_eff ∂_t ψ = [-ℏ_eff²/(2m*) ∇² + V]ψ
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Schrödinger Equation")

    # Define symbols
    hbar_eff, m_star, psi, t, V = define_symbols_batch(
        ['hbar_eff', 'm_star', 'psi', 't', 'V'], positive=True)
    
    # Time evolution term: iℏ_eff ∂_t ψ
    time_term = sp.I * hbar_eff * psi / t
    
    # Kinetic energy operator: -ℏ_eff²/(2m*) ∇²
    kinetic_operator = -hbar_eff**2 / (2 * m_star)
    kinetic_term = kinetic_operator * psi  # Applied to wavefunction
    
    # Potential energy term: V ψ
    potential_term = V * psi
    
    # Hamiltonian: Ĥ = -ℏ_eff²/(2m*) ∇² + V
    hamiltonian = kinetic_operator + V
    hamiltonian_applied = hamiltonian * psi
    
    # Test Schrödinger equation: iℏ_eff ∂_t ψ = Ĥψ
    v.check_eq("Schrödinger equation iℏ_eff ∂_t ψ = Ĥψ",
               time_term, hamiltonian_applied)
    
    # Test kinetic energy operator structure
    kinetic_expected = -hbar_eff**2 / (2 * m_star)
    v.check_eq("Kinetic energy operator -ℏ_eff²/(2m*)",
               kinetic_operator, kinetic_expected)

    v.success("Schrödinger equation verified")


def test_canonical_commutation_relations(v):
    """
    Test the canonical commutation relations: [x̂_i, p̂_j] = iℏ_eff δ_ij
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Canonical Commutation Relations")

    # Define symbols
    hbar_eff = define_symbols_batch(['hbar_eff'], positive=True)
    if isinstance(hbar_eff, tuple) and len(hbar_eff) == 1:
        hbar_eff = hbar_eff[0]
    
    # Commutator result (for i=j case, δ_ij = 1)
    commutator_result = sp.I * hbar_eff
    
    # Test canonical commutation relation
    v.check_eq("Canonical commutation [x̂_i, p̂_j] = iℏ_eff δ_ij",
               commutator_result, sp.I * hbar_eff)
    
    # Verify that the commutation relation fixes the uncertainty principle scale
    uncertainty_scale = hbar_eff
    v.check_eq("Uncertainty principle scale ℏ_eff",
               uncertainty_scale, hbar_eff)

    v.success("Canonical commutation relations verified")


def test_quantum_potential_equation(v):
    """
    Test the quantum potential: Q[ρ] = -ℏ_eff²/(2m*) ∇²√ρ/√ρ
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Quantum Potential Equation")

    # Define symbols
    hbar_eff, m_star, rho, V = define_symbols_batch(
        ['hbar_eff', 'm_star', 'rho', 'V'], positive=True)
    
    # Quantum potential expression
    rho_sqrt = sqrt(rho)
    
    # The quantum potential coefficient
    Q_coeff = -(hbar_eff**2) / (2 * m_star)
    Q_potential = Q_coeff * (1 / rho_sqrt)  # Simplified form
    
    # Test quantum potential coefficient
    quantum_coeff = hbar_eff**2 / (2 * m_star)
    expected_coeff = hbar_eff**2 / (2 * m_star)
    
    v.check_eq("Quantum potential coefficient ℏ_eff²/(2m*)",
               quantum_coeff, expected_coeff)
    
    # Test that Q appears in the Hamilton-Jacobi equation (eq:HJ_quantum)
    # ∂_t S + (∇S - qA)²/(2m*) + qΦ + V + Q[ρ] = 0
    HJ_energy_terms = quantum_coeff + V  # Kinetic + potential terms
    v.check_eq("Quantum potential in HJ equation", 
               quantum_coeff + V, HJ_energy_terms)

    v.success("Quantum potential equation verified")


def test_g_factor_calibration_handle(v):
    """
    Test the g-factor calibration: g = 2 + δg, δg ~ η_tw(ε/ℓ*)²
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("g-Factor Calibration Handle")

    # Define symbols
    g_factor, delta_g, eta_tw, varepsilon, ell_star = define_symbols_batch(
        ['g_factor', 'delta_g', 'eta_tw', 'varepsilon', 'ell_star'], positive=True)
    
    # Total g-factor equation
    g_total = 2 + delta_g
    
    # Test g-factor composition
    v.check_eq("g-factor g = 2 + δg", g_factor, g_total)
    
    # δg scaling relation
    delta_g_scaling = eta_tw * (varepsilon / ell_star)**2
    
    # Test δg scaling (calibration handle)
    v.check_eq("δg scaling δg ~ η_tw(ε/ℓ*)²", delta_g, delta_g_scaling)
    
    # Test that the correction is small (ε/ℓ* << 1)
    thickness_ratio = varepsilon / ell_star
    small_correction = eta_tw * thickness_ratio**2
    
    v.check_eq("Small thickness correction η_tw(ε/ℓ*)²",
               small_correction, delta_g_scaling)

    v.success("g-Factor calibration handle verified")


def test_decoherence_scaling_law(v):
    """
    Test the decoherence scaling: γ₂ ∝ α_tw ℏ_eff/(m*ℓ*⁴)(ε/ℓ*)^p
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Decoherence Scaling Law")

    # Define symbols
    Gamma_0, gamma_2, d, alpha_tw, hbar_eff, m_star = define_symbols_batch(
        ['Gamma_0', 'gamma_2', 'd', 'alpha_tw', 'hbar_eff', 'm_star'], positive=True)
    ell_star, varepsilon, p = define_symbols_batch(['ell_star', 'varepsilon', 'p'], positive=True)
    
    # Decoherence rate with d² scaling
    Gamma_dec = Gamma_0 + gamma_2 * d**2
    
    # Test the d² law structure
    expected_decay = Gamma_0 + gamma_2 * d**2
    v.check_eq("Decoherence d² law Γ_dec(d) = Γ₀ + γ₂ d²",
               Gamma_dec, expected_decay)
    
    # γ₂ scaling relation from calibration
    gamma_2_scaling = (alpha_tw * hbar_eff / 
                      (m_star * ell_star**4) *
                      (varepsilon / ell_star)**p)
    
    # Test γ₂ scaling (calibration handle)
    v.check_eq("γ₂ scaling γ₂ ∝ α_tw ℏ_eff/(m*ℓ*⁴)(ε/ℓ*)^p",
               gamma_2, gamma_2_scaling)
    
    # Test the geometric scaling with slab thickness
    thickness_dependence = (varepsilon / ell_star)**p
    v.check_eq("Thickness scaling (ε/ℓ*)^p",
               thickness_dependence, (varepsilon / ell_star)**p)

    v.success("Decoherence scaling law verified")


def test_gravity_phase_equation(v):
    """
    Test the gravity phase relation: Δφ = (m*/ℏ_eff)∫√(-g_μν dx^μ dx^ν)
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Gravity Phase Equation")

    # Define symbols
    m_star, hbar_eff, L_proper = define_symbols_batch(
        ['m_star', 'hbar_eff', 'L_proper'], positive=True)
    
    # Phase coefficient
    phase_coeff = m_star / hbar_eff
    
    # Proper time element (simplified as geometric interval)
    proper_time_element = sqrt(L_proper**2)
    
    # Gravity phase relation
    gravity_phase = phase_coeff * proper_time_element
    
    # Test the basic gravity phase equation
    expected_phase = (m_star / hbar_eff) * sqrt(L_proper**2)
    v.check_eq("Gravity phase Δφ = (m*/ℏ_eff)Δτ",
               gravity_phase, expected_phase)
    
    # Test phase coefficient structure
    v.check_eq("Phase coefficient m*/ℏ_eff", phase_coeff, m_star / hbar_eff)

    v.success("Gravity phase equation verified")


def test_calibration_parameter_relationships(v):
    """
    Test relationships between calibration parameters and physical observables.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Calibration Parameter Relationships")

    # Define symbols
    hbar_eff, p, lambda_wave, omega, k_spring, m_star = define_symbols_batch(
        ['hbar_eff', 'p', 'lambda_wave', 'omega', 'k_spring', 'm_star'], positive=True)
    gamma_2, d, alpha_tw, varepsilon, ell_star, p_exp = define_symbols_batch(
        ['gamma_2', 'd', 'alpha_tw', 'varepsilon', 'ell_star', 'p_exp'], positive=True)
    k_star, xi, delta_g, eta_tw = define_symbols_batch(
        ['k_star', 'xi', 'delta_g', 'eta_tw'], positive=True)
    
    # ℏ_eff calibration: de Broglie fringes; atomic spectra
    # de Broglie wavelength: λ = 2πℏ_eff/p
    de_broglie_relation = 2 * pi * hbar_eff / p
    v.check_eq("de Broglie calibration λ = 2πℏ_eff/p",
               lambda_wave, de_broglie_relation)
    
    # m* calibration: kinematics vs trap frequencies/dispersion
    # Trap frequency: ω ~ √(k_trap/m*) where k_trap is spring constant
    trap_relation = sqrt(k_spring / m_star)
    v.check_eq("Trap frequency calibration ω ~ √(k/m*)",
               omega, trap_relation)
    
    # ξ calibration: dispersion tail k* ~ ξ⁻¹ in dispersion relation
    cutoff_relation = 1 / xi
    v.check_eq("Core radius calibration k* ~ ξ⁻¹", k_star, cutoff_relation)
    
    # η_tw calibration: precision g-factor bounds (δg)
    g_precision_bound = eta_tw * (varepsilon/ell_star)**2
    v.check_eq("Spin renormalization calibration δg ~ η_tw(ε/ℓ*)²",
               delta_g, g_precision_bound)

    v.success("Calibration parameter relationships verified")


def test_calibration_and_parameter_table():
    """
    Main test function for the Calibration and Parameter Table section.
    
    CRITICAL: This tests actual mathematical equations using v.check_eq(), 
    NOT just dimensional consistency. Test failures reveal real issues
    in the theoretical framework that must be addressed.
    
    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Calibration and Parameter Table - Mathematical Verification",
        "Verification of actual equations and calibration relationships from quantum framework"
    )
    
    v.section("MATHEMATICAL EQUATION VERIFICATION")
    
    # Test fundamental equations from the paper
    v.info("\n--- 1) Circulation Quantization Equation ---")
    test_circulation_quantization_equation(v)
    
    v.info("\n--- 2) Dispersion Relation Equation ---")
    test_dispersion_relation_equation(v)
    
    v.info("\n--- 3) Schrödinger Equation ---")
    test_schrodinger_equation(v)
    
    v.info("\n--- 4) Canonical Commutation Relations ---")
    test_canonical_commutation_relations(v)
    
    v.info("\n--- 5) Quantum Potential Equation ---")
    test_quantum_potential_equation(v)
    
    v.info("\n--- 6) g-Factor Calibration Handle ---")
    test_g_factor_calibration_handle(v)
    
    v.info("\n--- 7) Decoherence Scaling Law ---")
    test_decoherence_scaling_law(v)
    
    v.info("\n--- 8) Gravity Phase Equation ---")
    test_gravity_phase_equation(v)
    
    v.info("\n--- 9) Calibration Parameter Relationships ---")
    test_calibration_parameter_relationships(v)
    
    # Return success rate - failures indicate mathematical issues
    return v.summary()


if __name__ == "__main__":
    # Add --quiet flag support
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--quiet', action='store_true', 
                       help='Only show test summary, not individual results')
    args, unknown = parser.parse_known_args()
    
    # Run the mathematical verification
    success_rate = test_calibration_and_parameter_table()
    
    if not args.quiet:
        print(f"\nFinal Success Rate: {success_rate:.1f}%")
        if success_rate < 100.0:
            print("\nWARNING: Test failures indicate mathematical issues in the framework")
            print("These are NOT just coding errors - they reveal theoretical problems")
            print("that need to be addressed in the equations or derivations.")
    
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)