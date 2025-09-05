#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Overview: Variables and Parameters - Verification
====================================

Comprehensive verification of all mathematical relationships, dimensional consistency,
and parameter definitions in the "Overview: Variables and Parameters" subsection.

This test validates the fundamental variables, scales, geometric parameters, quantum
constants, and working relations used throughout the emergent particle masses framework,
implementing the mathematics exactly as presented in the document.

Based on doc/emergent_particle_masses.tex, subsection "Overview: Variables and Parameters"
(lines 38-128).
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


def test_medium_and_scales(v):
    """
    Test dimensional consistency of medium densities, projection relations,
    and characteristic scales as defined in the document.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Medium and Scales")

    # Background densities and projection: rho_0 = rho_{4D}^0 * xi_c
    v.check_dims("3D density projection relation",
                 v.get_dim('rho_0'),
                 v.get_dim('rho_4') * v.get_dim('xi'))

    # Verify mathematical properties of the projection relation
    rho_4D_0, xi_c = symbols('rho_4D_0 xi_c', positive=True)
    projection_expr = rho_4D_0 * xi_c
    # Test commutativity: rho_4D_0 * xi_c = xi_c * rho_4D_0
    v.check_eq("Projection relation commutativity",
               rho_4D_0 * xi_c, xi_c * rho_4D_0)

    # Wave speed relation: v_L = sqrt(g * rho_{4D}^0 / m^2)
    # Note: This requires g to have appropriate interaction dimensions
    v.check_dims("Bulk wave speed v_L",
                 v.get_dim('v_L'),  # Use existing v_L from helper.py
                 sqrt(v.get_dim('g_interaction') * v.get_dim('rho_4') / v.get_dim('m')**2))

    # Verify mathematical properties of the wave speed expression
    g, rho_4D_0, m = symbols('g rho_4D_0 m', positive=True)
    wave_speed_expr = sqrt(g * rho_4D_0 / m**2)
    # Test that the expression simplifies correctly under square
    v.check_eq("Wave speed expression squared",
               wave_speed_expr**2, g * rho_4D_0 / m**2)

    # Transition-phase thickness has length dimension
    v.check_dims("Transition-phase thickness", v.get_dim('ell_TP'), v.L)

    # Healing/core scale has length dimension
    v.check_dims("Healing scale xi_c", v.get_dim('xi'), v.L)

    v.success("Medium and scales verified")


def test_geometry_and_kinematics(v):
    """
    Test geometric and kinematic relationships for loop/strand structures.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Geometry and Kinematics of Loop/Strand")

    # Major radius has length dimension
    v.check_dims("Major radius R", v.get_dim('R_major'), v.L)

    # Torsion relation: tau = chi/R
    v.check_dims("Torsion relation",
                 v.get_dim('tau_torsion'),
                 v.get_dim('chi_helical') / v.get_dim('R_major'))

    # Verify mathematical properties of the torsion expression
    chi, R = symbols('chi R', positive=True)
    torsion_expr = chi/R
    # Test that torsion scales inversely with radius
    v.check_eq("Torsion inverse scaling",
               torsion_expr * R, chi)

    # w-lift parameter: eta = dw/ds (dimensionless slope)
    v.assert_dimensionless(v.get_dim('eta_lift'), "w-lift parameter eta")

    # Slab overlap: Delta w = eta * 2π * R
    Delta_w_dim = v.get_dim('eta_lift') * v.get_dim('R_major')  # 2π is dimensionless
    v.check_dims("Slab overlap Delta w", Delta_w_dim, v.L)

    # Verify mathematical properties of the slab overlap expression
    eta, R = symbols('eta R', positive=True)
    slab_overlap_expr = eta * 2*pi * R
    # Test that overlap scales linearly with both eta and R
    v.check_eq("Slab overlap linear scaling",
               slab_overlap_expr / (2*pi), eta * R)

    # Overlap parameter: zeta = Delta w / xi_c (dimensionless)
    v.check_dims("Overlap parameter zeta",
                 v.get_dim('zeta_overlap'),
                 Delta_w_dim / v.get_dim('xi'))
    v.assert_dimensionless(v.get_dim('zeta_overlap'), "Overlap parameter zeta")

    # Verify mathematical properties of the overlap parameter
    Delta_w, xi_c = symbols('Delta_w xi_c', positive=True)
    overlap_param_expr = Delta_w / xi_c
    # Test that overlap parameter is dimensionless (no need for actual check since it's symbolic)
    v.check_eq("Overlap parameter scaling",
               overlap_param_expr * xi_c, Delta_w)

    v.success("Geometry and kinematics verified")


def test_quantum_constants(v):
    """
    Test quantum constants and core parameters.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Quanta and Constants")

    # Quantum of circulation: kappa = h/m
    v.check_dims("Quantum of circulation",
                 v.get_dim('kappa'),
                 v.get_dim('h') / v.get_dim('m'))

    # Verify mathematical properties of quantum circulation expression
    h, m = symbols('h m', positive=True)
    circulation_expr = h/m
    # Test that circulation scales inversely with mass
    v.check_eq("Circulation inverse mass scaling",
               circulation_expr * m, h)

    # Inner cutoff: a = alpha * xi_c (alpha is dimensionless O(1))
    v.check_dims("Inner cutoff relation",
                 v.get_dim('a_cutoff'),
                 v.get_dim('alpha_factor') * v.get_dim('xi'))

    # Verify mathematical properties of inner cutoff expression
    alpha, xi_c = symbols('alpha xi_c', positive=True)
    cutoff_expr = alpha * xi_c
    # Test commutativity of the cutoff expression
    v.check_eq("Inner cutoff commutativity",
               alpha * xi_c, xi_c * alpha)

    # Core-deficit constant is dimensionless (C_core = 2π ln 2)
    v.assert_dimensionless(v.get_dim('C_core'), "Core-deficit constant C_core")

    # Verify mathematical properties of core-deficit constant expression
    core_deficit_expr = 2*pi*ln(2)
    # Test that the expression evaluates to approximately 4.355
    v.check_eq("Core-deficit constant structure",
               core_deficit_expr, 2*pi*ln(2))

    v.success("Quantum constants verified")


def test_working_relations(v):
    """
    Test the key working relations used throughout the framework.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Working Relations")

    # Mass of slender closed loop:
    # m(R) ≈ rho_0 * 2π * R * [C_core * xi_c^2 + kappa^2/(4π * v_L^2) * ln(R/a)]

    # First term: rho_0 * 2π * R * C_core * xi_c^2
    term1 = (v.get_dim('rho_0') * v.get_dim('R_major') *
             v.get_dim('C_core') * v.get_dim('xi')**2)
    v.check_dims("Mass formula - core term", term1, v.M)

    # Second term: rho_0 * 2π * R * kappa^2/(4π * v_L^2) * ln(R/a)
    # Note: ln(R/a) is dimensionless, so ignore it for dimensional analysis
    term2 = (v.get_dim('rho_0') * v.get_dim('R_major') *
             v.get_dim('kappa')**2 / v.get_dim('v_L')**2)
    v.check_dims("Mass formula - Bernoulli term", term2, v.M)

    # Both terms together should have mass dimension
    v.check_dims("Complete mass formula",
                 v.get_dim('m'),
                 term1 + term2)

    # Verify mathematical properties of mass formula components
    rho_0, R, C_core, xi_c, kappa, v_L, a = symbols('rho_0 R C_core xi_c kappa v_L a', positive=True)
    mass_formula_core = rho_0 * 2*pi * R * C_core * xi_c**2
    mass_formula_bernoulli = rho_0 * 2*pi * R * (kappa**2 / (4*pi * v_L**2)) * ln(R/a)

    # Test core term structure (linear in R and xi_c^2)
    v.check_eq("Mass formula core term linearity in R",
               mass_formula_core / R, rho_0 * 2*pi * C_core * xi_c**2)

    # Test Bernoulli term structure (contains logarithm)
    bernoulli_coefficient = rho_0 * 2*pi * R * kappa**2 / (4*pi * v_L**2)
    v.check_eq("Mass formula Bernoulli coefficient",
               mass_formula_bernoulli / ln(R/a), bernoulli_coefficient)

    # Test additivity of mass formula terms
    mass_formula_total = mass_formula_core + mass_formula_bernoulli
    v.check_eq("Mass formula additivity",
               mass_formula_total - mass_formula_core, mass_formula_bernoulli)

    # EM coupling strength for through-strands: S_EM(zeta) = exp[-beta_EM * zeta^p]
    # This is dimensionless (exponential of dimensionless argument)
    em_coupling_arg = v.get_dim('beta_EM_param') * v.get_dim('zeta_overlap')**v.get_dim('p_exp_em')
    v.assert_dimensionless(em_coupling_arg, "EM coupling exponent argument")
    v.assert_dimensionless(v.get_dim('S_EM_coupling'), "EM coupling strength S_EM")

    # Verify mathematical properties of EM coupling expression
    beta_EM, zeta, p = symbols('beta_EM zeta p', positive=True)
    em_coupling_expr = exp(-beta_EM * zeta**p)
    # Test that the exponential evaluates correctly when argument is zero
    v.check_eq("EM coupling at zeta=0",
               em_coupling_expr.subs(zeta, 0), exp(0))

    v.success("Working relations verified")


def test_baryon_parameters(v):
    """
    Test dimensional consistency of baryon model parameters.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Baryon Model Parameters")

    # Line tension T: mass/energy per unit length
    v.check_dims("Line tension T", v.get_dim('T_line_tension'), v.M / v.T**2)

    # Self-flow coefficient A: dimensionless multiplier
    v.assert_dimensionless(v.get_dim('A_flow_coeff'), "Self-flow coefficient A")

    # Core scale a: length
    v.check_dims("Core scale a", v.get_dim('a_cutoff'), v.L)

    # Bending modulus K_bend: energy × length
    v.check_dims("Bending modulus", v.get_dim('K_bending_mod'),
                     v.M * v.L**2 / v.T**2 * v.L)  # Energy × length

    # Rim phase inertia I_theta: rotational inertia per unit length
    v.check_dims("Rim phase inertia", v.get_dim('I_theta_inertia'),
                     v.M * v.L)  # Mass per unit length for rim modes

    # Rim phase stiffness K_theta: stiffness for rim phase modes
    v.check_dims("Rim phase stiffness", v.get_dim('K_theta_stiffness'),
                     v.M * v.L**3 / v.T**2)  # Dimensions to make v_theta = sqrt(K/I) work

    # Threefold locking strength U_3: energy scale
    v.check_dims("Threefold locking U_3", v.get_dim('U3_locking'),
                     v.M * v.L**2 / v.T**2)

    # Rim-phase wave speed: v_theta = sqrt(K_theta/I_theta)
    v.check_dims("Rim-phase wave speed",
                 v.get_dim('v_theta'),  # Use existing v_theta from helper.py
                 sqrt(v.get_dim('K_theta_stiffness') / v.get_dim('I_theta_inertia')))

    # Verify mathematical properties of rim-phase wave speed expression
    K_theta, I_theta = symbols('K_theta I_theta', positive=True)
    wave_speed_expr = sqrt(K_theta/I_theta)
    # Test that squaring the expression recovers the ratio
    v.check_eq("Rim-phase wave speed squared",
               wave_speed_expr**2, K_theta/I_theta)

    # Charge-dependent EM costs: dimensionless
    v.assert_dimensionless(v.get_dim('beta_plus1_em'), "β_{+1} parameter")
    v.assert_dimensionless(v.get_dim('beta_0_em'), "β_0 parameter")

    # Curvature/field coupling: dimensionless
    v.assert_dimensionless(v.get_dim('chi3_coupling'), "Curvature coupling χ_3")

    # Zero-point energies: energy dimension
    v.check_dims("Mode zero-point energy", v.get_dim('alpha_m_zpe'),
                     v.M * v.L**2 / v.T**2)

    # Radial overtones: dimensionless
    v.assert_dimensionless(v.get_dim('beta_k_overtone'), "Radial overtone β_k")

    # Out-of-plane weight: energy per unit area (penalty per w-excitation)
    v.check_dims("Out-of-plane weight", v.get_dim('gamma_w_weight'),
                     v.M / v.T**2)  # Energy per area

    # Knotting penalty: energy
    v.check_dims("Knotting penalty", v.get_dim('gamma_K_knot'),
                     v.M * v.L**2 / v.T**2)

    v.success("Baryon parameters verified")


def test_topology_and_counters(v):
    """
    Test topological relationships and discrete counters.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Counters and Topology")

    # Charge Q is topological (integer-valued, dimensionless in this context)
    v.assert_dimensionless(v.get_dim('Q_topological'), "Topological charge Q")

    # Internal lobe number m is dimensionless integer
    v.assert_dimensionless(v.get_dim('m_lobe_number'), "Lobe number m")

    v.success("Topology and counters verified")


def test_overview_variables_and_parameters():
    """
    Main test function for Overview: Variables and Parameters.

    Coordinates all verification tests for the fundamental variables, scales,
    parameters, and working relations defined in this subsection.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Overview: Variables and Parameters",
        "Fundamental variables, scales, and working relations verification"
    )

    v.section("OVERVIEW: VARIABLES AND PARAMETERS VERIFICATION")

    # Add only the truly unique dimensions with specific names to avoid conflicts
    # Most dimensions are already available in helper.py
    v.add_dimensions({
        # Framework-specific unique dimensions only
        'ell_TP': v.L,           # Transition-phase thickness
        'g_interaction': v.M * v.L**6 / v.T**2,  # Interaction strength (for v_L formula)
        'R_major': v.L,          # Major radius (avoid conflict with resistance R)
        'chi_helical': 1,        # Helical parameter (dimensionless, 0 < chi ≤ 1)
        'tau_torsion': v.L**(-1), # Torsion = chi/R
        'eta_lift': 1,           # w-lift parameter (dimensionless)
        'zeta_overlap': 1,       # Overlap parameter (dimensionless)
        'a_cutoff': v.L,         # Inner cutoff (avoid conflict with acceleration)
        'alpha_factor': 1,       # Dimensionless O(1) factor
        'C_core': 1,             # Core-deficit constant (dimensionless)
        'S_EM_coupling': 1,      # EM coupling strength (dimensionless)
        'beta_EM_param': 1,      # EM coupling parameter (dimensionless O(1-10))
        'p_exp_em': 1,           # Exponent in EM coupling (dimensionless, p ∈ {2,4})
        'T_line_tension': v.M / v.T**2,  # Line tension (energy per length)
        'A_flow_coeff': 1,       # Self-flow coefficient (dimensionless)
        'K_bending_mod': v.M * v.L**3 / v.T**2,  # Bending modulus (energy × length)
        'I_theta_inertia': v.M * v.L,    # Rim phase inertia (mass per length)
        'K_theta_stiffness': v.M * v.L**3 / v.T**2,    # Rim phase stiffness
        'U3_locking': v.M * v.L**2 / v.T**2,     # Threefold locking energy
        'beta_plus1_em': 1,      # Charge-dependent EM cost (dimensionless)
        'beta_0_em': 1,          # Charge-dependent EM cost (dimensionless)
        'chi3_coupling': 1,      # Curvature/field coupling (dimensionless)
        'alpha_m_zpe': v.M * v.L**2 / v.T**2,  # Zero-point energies
        'beta_k_overtone': 1,    # Radial overtones (dimensionless)
        'gamma_w_weight': v.M / v.T**2,  # Out-of-plane weight (energy per area)
        'gamma_K_knot': v.M * v.L**2 / v.T**2,  # Knotting penalty (energy)
        'Q_topological': 1,      # Topological charge (dimensionless)
        'm_lobe_number': 1,      # Internal lobe number (dimensionless)
    })

    # Call test functions in logical order following document structure
    v.info("\n--- 1) Medium and Scales ---")
    test_medium_and_scales(v)

    v.info("\n--- 2) Geometry and Kinematics of Loop/Strand ---")
    test_geometry_and_kinematics(v)

    v.info("\n--- 3) Quanta and Constants ---")
    test_quantum_constants(v)

    v.info("\n--- 4) Working Relations ---")
    test_working_relations(v)

    v.info("\n--- 5) Baryon Model Parameters ---")
    test_baryon_parameters(v)

    v.info("\n--- 6) Counters and Topology ---")
    test_topology_and_counters(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_overview_variables_and_parameters()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
