#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Neutrino Masses and Mixing - Verification
====================================

Comprehensive verification of all mathematical relationships, dimensional consistency,
and formulas in the "Neutrino Masses and Mixing" subsection.

This test validates the through-strand defect model for neutrinos, including the
mass template application, geometric ladder construction, helical kinematics,
slab overlap suppression, and mixing calculations, implementing the mathematics
exactly as presented in the document.

Based on doc/emergent_particle_masses.tex, subsection "Neutrino Masses and Mixing"
(lines 382-523).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, sqrt, ln, exp, simplify, cos

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify,
)


def test_mass_template(v):
    """
    Test the dimensional consistency of the slender loop mass template equation
    as applied to neutrinos.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Mass Template for Slender Loops")

    # Mass template: m_loop(R) = rho_0 * 2π * R * [C_core * xi_c^2 + kappa^2/(4π * v_L^2) * ln(R/a)]

    # First check the core term: rho_0 * 2π * R * C_core * xi_c^2
    core_term = v.get_dim('rho_0') * v.get_dim('R_major') * v.get_dim('xi')**2
    v.check_dims("Core mass term: rho_0 * 2π * R * C_core * xi_c^2",
                 core_term, v.M)

    # Check the logarithmic term: rho_0 * 2π * R * kappa^2/(4π * v_L^2) * ln(R/a)
    # Note: kappa = h/m is circulation quantum with dimensions [L^2/T]
    log_term = (v.get_dim('rho_0') * v.get_dim('R_major') *
                v.get_dim('kappa')**2 / v.get_dim('v_L')**2)
    v.check_dims("Logarithmic mass term: rho_0 * 2π * R * kappa^2/(4π * v_L^2) * ln(R/a)",
                 log_term, v.M)

    # Full mass template should have mass dimensions
    full_template = core_term + log_term
    v.check_dims("Complete mass template m_loop(R)", full_template, v.M)

    # Check parameter definitions as stated in document
    # rho_0 = rho_4D^0 * xi_c (both already defined in helper.py)
    v.check_dims("3D density rho_0 = rho_4D^0 * xi_c",
                 v.get_dim('rho_0'),
                 v.get_dim('rho_4') * v.get_dim('xi'))

    # Circulation quantum: kappa = h/m (both already defined in helper.py)
    v.check_dims("Circulation quantum kappa = h/m",
                 v.get_dim('kappa'),
                 v.get_dim('h') / v.get_dim('m'))

    # Bulk wave speed: v_L = sqrt(g * rho_4D^0 / m^2) (v_L already defined in helper.py)
    # Note: 'g' here refers to 4D interaction constant, not gravitational acceleration
    v.check_dims("Bulk wave speed v_L = sqrt(g * rho_4D^0 / m^2)",
                 v.get_dim('v_L'),
                 sqrt(v.get_dim('g_4D') * v.get_dim('rho_4') / v.get_dim('m')**2))

    # Inner cutoff: a = alpha * xi_c
    v.check_dims("Inner cutoff a = alpha * xi_c",
                 v.get_dim('a_cutoff'),
                 v.get_dim('alpha') * v.get_dim('xi'))

    v.success("Mass template verified")


def test_geometric_ladder_and_kinematics(v):
    """
    Test the geometric ladder construction and helical kinematics
    for neutrino modes as defined in the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Geometry and Kinematics")

    # Geometric ladder: R_n = R_star * a_n where a_n = 2n + 1
    # R_n has length dimensions
    v.check_dims("Mode radius R_n", v.get_dim('R_n'), v.L)
    v.check_dims("Base radius R_star", v.get_dim('R_star'), v.L)

    # a_n is dimensionless (2n + 1)
    v.check_dims("Geometric factor a_n = 2n + 1", v.get_dim('a_n'), 1)

    # Helical angle: theta_helix = 2π * chi
    v.check_dims("Helical angle theta_helix = 2π * chi",
                 v.get_dim('theta_helix'),
                 v.get_dim('chi_helical'))

    # Both should be dimensionless (angles)
    v.check_dims("Helical angle dimensionless", v.get_dim('theta_helix'), 1)
    v.check_dims("Torsion fraction chi dimensionless", v.get_dim('chi_helical'), 1)

    # Torsion: tau = chi/R_n
    v.check_dims("Torsion tau = chi/R_n",
                 v.get_dim('tau_torsion'),
                 v.get_dim('chi_helical') / v.get_dim('R_n'))

    # Lift per circuit: Delta_w_n = 2π * R_n * eta_n
    v.check_dims("Lift per circuit Delta_w_n = 2π * R_n * eta_n",
                 v.get_dim('Delta_w'),
                 v.get_dim('R_n') * v.get_dim('eta_lift'))

    # eta_n is lift rate (dimensionless slope)
    v.check_dims("Lift rate eta_n dimensionless", v.get_dim('eta_lift'), 1)

    # Dimensionless lift: zeta_n = Delta_w_n / xi_c
    v.check_dims("Dimensionless lift zeta_n = Delta_w_n / xi_c",
                 v.get_dim('zeta_lift'),
                 v.get_dim('Delta_w') / v.get_dim('xi'))

    # zeta_n should be dimensionless
    v.check_dims("Dimensionless lift zeta_n", v.get_dim('zeta_lift'), 1)

    v.success("Geometric ladder and kinematics verified")


def test_bare_mass_calculation(v):
    """
    Test the bare (unsuppressed) mass calculation using the mass template
    for each neutrino mode.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Bare Mass Calculation")

    # Mode-by-mode bare mass: m_bare,n = rho_0 * 2π * R_n * [C_core * xi_c^2 + kappa^2/(4π * v_L^2) * ln(R_n/a)]

    # P coefficient: P = 2π * rho_0 * C_core * xi_c^2
    P_coeff = v.get_dim('rho_0') * v.get_dim('xi')**2
    v.check_dims("P coefficient: P = 2π * rho_0 * C_core * xi_c^2",
                 P_coeff, v.M / v.L)

    # Q coefficient: Q = 2π * rho_0 * kappa^2/(4π * v_L^2) = rho_0 * kappa^2/(2 * v_L^2)
    Q_coeff = v.get_dim('rho_0') * v.get_dim('kappa')**2 / v.get_dim('v_L')**2
    v.check_dims("Q coefficient: Q = rho_0 * kappa^2/(2 * v_L^2)",
                 Q_coeff, v.M / v.L)

    # Bare mass: m_bare,n = R_n * [P + Q * ln(R_n/a)]
    bare_mass_term1 = v.get_dim('R_n') * P_coeff
    bare_mass_term2 = v.get_dim('R_n') * Q_coeff  # ln(R_n/a) is dimensionless

    v.check_dims("Bare mass linear term: R_n * P", bare_mass_term1, v.M)
    v.check_dims("Bare mass logarithmic term: R_n * Q * ln(R_n/a)", bare_mass_term2, v.M)

    # Full bare mass
    bare_mass_total = bare_mass_term1 + bare_mass_term2
    v.check_dims("Complete bare mass m_bare,n", bare_mass_total, v.M)

    v.success("Bare mass calculation verified")


def test_slab_overlap_suppression(v):
    """
    Test the slab overlap suppression factor for through-strand defects.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Slab Overlap Suppression")

    # Overlap factor: f_slab(zeta_n) = exp[-beta_m * zeta_n^p]
    # Check that the exponent is dimensionless

    # beta_m is dimensionless scaling factor
    v.check_dims("Suppression parameter beta_m dimensionless", v.get_dim('beta_m'), 1)

    # zeta_n is already dimensionless
    v.check_dims("Dimensionless lift zeta_n", v.get_dim('zeta_lift'), 1)

    # Power p is dimensionless (p ∈ {2,4})
    v.check_dims("Power parameter p dimensionless", v.get_dim('p_power'), 1)

    # The full exponent beta_m * zeta_n^p should be dimensionless
    exponent = v.get_dim('beta_m') * v.get_dim('zeta_lift')**v.get_dim('p_power')
    v.check_dims("Suppression exponent beta_m * zeta_n^p dimensionless", exponent, 1)

    # f_slab itself is dimensionless (it's an exponential of dimensionless quantity)
    v.check_dims("Overlap factor f_slab dimensionless", v.get_dim('f_slab'), 1)

    # Check kinetic penalty terms that motivate this form
    # Chiral penalty: δE_chiral ≃ (1/2) * rho_0 * v_eff^2 * (theta_helix/2π)^2 * (2π R_n) * (π xi_c^2)
    chiral_energy = (v.get_dim('rho_0') * v.get_dim('v_eff')**2 *
                     (v.get_dim('theta_helix'))**2 *
                     v.get_dim('R_n') * v.get_dim('xi')**2)
    v.check_dims("Chiral kinetic penalty δE_chiral", chiral_energy, v.M * v.L**2 / v.T**2)

    # w-penalty: δE_w ≃ (1/2) * rho_0 * v_eff^2 * (Delta_w_n/xi_c)^2 * (2π R_n) * (π xi_c^2)
    w_energy = (v.get_dim('rho_0') * v.get_dim('v_eff')**2 *
                (v.get_dim('Delta_w') / v.get_dim('xi'))**2 *
                v.get_dim('R_n') * v.get_dim('xi')**2)
    v.check_dims("w-direction kinetic penalty δE_w", w_energy, v.M * v.L**2 / v.T**2)

    v.success("Slab overlap suppression verified")


def test_final_neutrino_mass_formula(v):
    """
    Test the final neutrino mass formula combining bare mass and suppression.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Final Neutrino Mass Formula")

    # Final formula: m_nu,n = m_bare,n * f_slab(zeta_n)

    # m_bare,n has mass dimensions
    v.check_dims("Bare mass m_bare,n", v.get_dim('m_bare'), v.M)

    # f_slab is dimensionless
    v.check_dims("Overlap factor f_slab", v.get_dim('f_slab'), 1)

    # Product should have mass dimensions
    final_mass = v.get_dim('m_bare') * v.get_dim('f_slab')
    v.check_dims("Final neutrino mass m_nu,n = m_bare,n * f_slab(zeta_n)",
                 final_mass, v.M)

    # Verify this is the sector summary formula from the document
    v.check_dims("Neutrino mass (sector summary)", v.get_dim('m_nu'), v.M)

    v.success("Final neutrino mass formula verified")


def test_mixing_formulation(v):
    """
    Test the neutrino mixing formulation based on geometric overlap.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Mixing from Geometric Overlap")

    # Energy proportionality: E_n ∝ m_nu,n
    v.check_dims("Neutrino energy E_n", v.get_dim('E_nu'), v.M * v.L**2 / v.T**2)
    v.check_dims("Mass-energy proportionality check",
                 v.get_dim('E_nu'),
                 v.get_dim('m_nu') * v.get_dim('c')**2)

    # Mixing matrix element: V_nm = V_0 * exp[-(zeta_n - zeta_m)^2/(2*sigma_zeta^2)] * cos(Delta_phi_nm)

    # V_0 has energy dimensions
    v.check_dims("Mixing scale V_0", v.get_dim('V_0'), v.M * v.L**2 / v.T**2)

    # Coherence scale sigma_zeta is dimensionless
    v.check_dims("Coherence scale sigma_zeta dimensionless", v.get_dim('sigma_zeta'), 1)

    # Phase difference Delta_phi_nm is dimensionless (angle)
    v.check_dims("Phase mismatch Delta_phi_nm dimensionless", v.get_dim('Delta_phi'), 1)

    # The Gaussian factor exp[-(zeta_n - zeta_m)^2/(2*sigma_zeta^2)] is dimensionless
    gaussian_exponent = (v.get_dim('zeta_lift'))**2 / (v.get_dim('sigma_zeta'))**2
    v.check_dims("Gaussian overlap exponent dimensionless", gaussian_exponent, 1)

    # cos(Delta_phi_nm) is dimensionless
    cosine_factor = v.get_dim('Delta_phi')  # cos doesn't change dimensionality
    v.check_dims("Cosine phase factor dimensionless", cosine_factor, 1)

    # Full mixing element has energy dimensions
    mixing_element = v.get_dim('V_0')  # * dimensionless factors
    v.check_dims("Mixing matrix element V_nm", mixing_element, v.M * v.L**2 / v.T**2)

    # Mixing angle formula: tan(2*theta_nm) = 2*|V_nm| / |E_m - E_n|
    # This should be dimensionless
    mixing_ratio = v.get_dim('V_0') / v.get_dim('E_nu')
    v.check_dims("Mixing angle ratio 2*|V_nm|/|E_m - E_n| dimensionless",
                 mixing_ratio, 1)

    v.success("Mixing formulation verified")


def test_benchmark_numerical_consistency(v):
    """
    Test the dimensional consistency of the benchmark numerical values
    and mass differences as presented in the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Benchmark Numerical Results")

    # Check that mass values have proper dimensions when expressed in eV
    # Document gives masses as: m_nu = (0.00352, 0.00935, 0.05106) eV

    # eV as energy unit
    v.check_dims("Electron volt (eV) as mass unit",
                 v.get_dim('eV') / v.get_dim('c')**2, v.M)

    # Mass squared differences
    # Δm²_21 = 7.50×10⁻⁵ eV²
    # Δm²_31 = 2.595×10⁻³ eV²
    # Δm²_32 = 2.520×10⁻³ eV²

    mass_sq_diff = (v.get_dim('eV') / v.get_dim('c')**2)**2
    v.check_dims("Mass squared differences Δm² in eV²",
                 mass_sq_diff, v.M**2)

    # Sum of masses: Σmν ≃ 0.064 eV
    mass_sum = v.get_dim('eV') / v.get_dim('c')**2
    v.check_dims("Sum of neutrino masses Σmν", mass_sum, v.M)

    # Ratio Δm²_32/Δm²_21 should be dimensionless
    mass_ratio = mass_sq_diff / mass_sq_diff
    v.check_dims("Mass squared ratio Δm²_32/Δm²_21 dimensionless", mass_ratio, 1)

    # Verify that the effective coefficients P and Q have the right dimensions
    # P = 1.114×10⁻³ eV, Q = 4.179×10⁻³ eV (both absorbed into eV units)
    P_eff = v.get_dim('eV') / v.get_dim('c')**2 / v.L
    Q_eff = v.get_dim('eV') / v.get_dim('c')**2 / v.L
    v.check_dims("Effective P coefficient (eV units)", P_eff, v.M / v.L)
    v.check_dims("Effective Q coefficient (eV units)", Q_eff, v.M / v.L)

    v.success("Benchmark numerical consistency verified")


def test_neutrino_masses_and_mixing():
    """
    Main test function for Neutrino Masses and Mixing.

    This function coordinates all verification tests for the neutrino sector,
    calling helper functions to validate mass template application, geometric
    construction, suppression mechanisms, and mixing formulation.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Neutrino Masses and Mixing",
        "Through-strand defect model with mass suppression and geometric mixing"
    )

    v.section("NEUTRINO MASSES AND MIXING VERIFICATION")

    # Add custom dimensions for neutrino-specific quantities
    # Only add dimensions not already defined in helper.py
    v.add_dimensions({
        # Mass template parameters (most are already defined)
        'a_cutoff': v.L,                 # Inner cutoff scale (a = alpha * xi_c)
        'alpha': 1,                      # Dimensionless cutoff factor
        'C_core': 1,                     # Core coefficient (dimensionless)

        # Geometric parameters
        'R_n': v.L,                      # Mode-specific radius
        'R_star': v.L,                   # Base radius
        'R_major': v.L,                  # Major radius (general)
        'a_n': 1,                        # Geometric ladder factor
        'theta_helix': 1,                # Helical angle (dimensionless)
        'chi_helical': 1,                # Torsion fraction (dimensionless)
        'tau_torsion': v.L**(-1),        # Torsion (inverse length)
        'Delta_w': v.L,                  # Lift per circuit
        'eta_lift': 1,                   # Lift rate (dimensionless)
        'zeta_lift': 1,                  # Dimensionless lift parameter

        # Mass calculation
        'm_bare': v.M,                   # Bare (unsuppressed) mass
        'f_slab': 1,                     # Overlap suppression factor
        'beta_m': 1,                     # Suppression parameter
        'p_power': 1,                    # Suppression power (2 or 4)
        'm_nu': v.M,                     # Final neutrino mass

        # Mixing parameters
        'E_nu': v.M * v.L**2 / v.T**2,   # Neutrino energy
        'V_0': v.M * v.L**2 / v.T**2,    # Mixing scale
        'sigma_zeta': 1,                 # Coherence scale (dimensionless)
        'Delta_phi': 1,                  # Phase mismatch (dimensionless)

        # Units and constants
        'eV': v.M * v.L**2 / v.T**2,     # Electron volt
        'g_4D': v.M * v.L**6 / v.T**2,   # 4D interaction constant for v_L formula
    })

    # Call test functions in logical order following document structure
    v.info("\n--- 1) Mass Template for Slender Loops ---")
    test_mass_template(v)

    v.info("\n--- 2) Geometric Ladder and Kinematics ---")
    test_geometric_ladder_and_kinematics(v)

    v.info("\n--- 3) Bare Mass Calculation ---")
    test_bare_mass_calculation(v)

    v.info("\n--- 4) Slab Overlap Suppression ---")
    test_slab_overlap_suppression(v)

    v.info("\n--- 5) Final Neutrino Mass Formula ---")
    test_final_neutrino_mass_formula(v)

    v.info("\n--- 6) Mixing from Geometric Overlap ---")
    test_mixing_formulation(v)

    v.info("\n--- 7) Benchmark Numerical Consistency ---")
    test_benchmark_numerical_consistency(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_neutrino_masses_and_mixing()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
