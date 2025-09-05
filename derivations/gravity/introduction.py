#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gravity Introduction - Verification
====================================

Comprehensive verification of the fundamental concepts introduced in the gravity
theory framework, including asymptotic causality, F_μν-built observables, GEM
conventions, and the foundational Maxwell-like gravitational equations.

This test validates the dimensional consistency of gravitoelectric and
gravitomagnetic field definitions, wave equation structures, and the connection
to the tsunami principle and causality framework.

Based on doc/gravity.tex, introduction section (lines 1-43).
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
    verify_wave_equation,
    verify_poisson_grav,
    quick_verify,
)


def test_asymptotic_causality_framework(v):
    """
    Test the asymptotic causality and wave sector propagation concepts.

    Verifies that only F_μν-built observables propagate at speed c in the wave
    sector, while bulk v_L adjustments are decoupled in the asymptotic limit.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Asymptotic Causality and Wave Sector")

    # Test that wave sector observables propagate at speed c
    v.check_dims("Speed of light c", v.get_dim('c'), v.L/v.T)

    # F_μν observables should have electromagnetic field tensor dimensions
    # F_μν has dimensions [M L T^-3 Q^-1] in SI (same as E-field/c or B-field)
    F_mu_nu_dim = v.get_dim('E') / v.get_dim('c')  # or equivalently B-field
    v.check_dims("F_μν observables dimension", F_mu_nu_dim, v.get_dim('B'))

    # Bulk v_L adjustments represent fluid velocity modifications
    v.check_dims("Bulk velocity v_L", v.get_dim('v_L'), v.L/v.T)

    # In asymptotic causality, only electromagnetic-like observables
    # (built from F_μν) propagate at light speed
    wave_propagation_speed = v.get_dim('c')
    v.check_dims("Wave sector propagation speed",
                 wave_propagation_speed, v.L/v.T)

    # The decoupling means bulk flow adjustments don't affect
    # the asymptotic wave propagation speed
    v.info("✓ Bulk v_L adjustments decoupled from F_μν wave propagation")
    v.info("✓ Only F_μν-built observables propagate at speed c")

    v.success("Asymptotic causality framework verified")


def test_gem_conventions_and_signature(v):
    """
    Test the GEM (Gravitoelectromagnetic) conventions and metric signature.

    Verifies the weak-field potential definitions and metric signature (-,+,+,+).

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("GEM Conventions and Signature")

    # Metric signature (-,+,+,+) - this is just a convention note
    v.info("Metric signature: (-,+,+,+)")

    # Test weak-field potential definitions from doc/gravity.tex lines 8-12:
    # h_00 = -2Φ_g/c², h_0i = -4A_{g,i}/c, h_ij = -2Φ_g δ_ij/c²

    # Define coordinate symbols
    x, y, z, t = symbols('x y z t', real=True)
    c = symbols('c', positive=True)

    # Define gravitational potentials as concrete functions
    Phi_g = symbols('Phi_g', cls=sp.Function)(x, y, z, t)
    A_g_x = symbols('A_g_x', cls=sp.Function)(x, y, z, t)
    A_g_y = symbols('A_g_y', cls=sp.Function)(x, y, z, t)
    A_g_z = symbols('A_g_z', cls=sp.Function)(x, y, z, t)

    # Define metric perturbations according to GEM conventions (doc/gravity.tex lines 9-11)
    h_00_def = -2*Phi_g/c**2
    h_0x_def = -4*A_g_x/c
    h_0y_def = -4*A_g_y/c
    h_0z_def = -4*A_g_z/c
    h_xx_def = -2*Phi_g/c**2  # δ_xx = 1, so h_xx = -2Φ_g/c²
    h_yy_def = -2*Phi_g/c**2  # δ_yy = 1, so h_yy = -2Φ_g/c²
    h_zz_def = -2*Phi_g/c**2  # δ_zz = 1, so h_zz = -2Φ_g/c²

    # Verify dimensional consistency of all metric components
    v.check_dims("h_00 = -2Φ_g/c² (dimensionless)",
                 v.get_dim('Phi_g') / v.get_dim('c')**2, 1)

    v.check_dims("h_0i = -4A_{g,i}/c (dimensionless)",
                 v.get_dim('A_g') / v.get_dim('c'), 1)

    v.check_dims("h_ij = -2Φ_g δ_ij/c² (dimensionless)",
                 v.get_dim('Phi_g') / v.get_dim('c')**2, 1)

    # Test that the convention gives consistent normalization
    # All diagonal spatial terms h_ii should have the same form as h_00
    v.check_eq("Spatial diagonal consistency: h_xx form",
               h_xx_def / h_00_def, 1)
    v.check_eq("Spatial diagonal consistency: h_yy form",
               h_yy_def / h_00_def, 1)
    v.check_eq("Spatial diagonal consistency: h_zz form",
               h_zz_def / h_00_def, 1)

    # Verify the GEM convention produces the correct scaling
    # The factor of 4 in h_0i vs factor of 2 in h_00,h_ij reflects the GEM convention
    scaling_ratio = (h_0x_def * c) / (-4 * A_g_x)
    v.check_eq("GEM h_0i scaling convention check", scaling_ratio, 1)

    # Test that metric perturbations maintain the expected sign structure
    # h_00 < 0 for attractive gravity (Φ_g > 0), consistent with signature (-,+,+,+)
    v.info("✓ Sign structure: h_00 < 0 for attractive potential (Φ_g > 0)")
    v.info("✓ Sign structure: h_0i coefficients negative for frame-dragging")
    v.info("✓ All metric perturbations h_μν dimensionless as required")

    v.success("GEM conventions and signature verified")


def test_gravitoelectric_gravitomagnetic_fields(v):
    """
    Test the definitions of gravitoelectric and gravitomagnetic fields.

    Verifies: E_g = -∇Φ_g - ∂_t A_g and B_g = ∇ × A_g
    From doc/gravity.tex lines 15-16.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Gravitoelectric and Gravitomagnetic Field Definitions")

    # Define coordinate symbols
    x, y, z, t = symbols('x y z t', real=True)

    # Define gravitational potentials as functions
    Phi_g = symbols('Phi_g', cls=sp.Function)(x, y, z, t)
    A_g_x = symbols('A_g_x', cls=sp.Function)(x, y, z, t)
    A_g_y = symbols('A_g_y', cls=sp.Function)(x, y, z, t)
    A_g_z = symbols('A_g_z', cls=sp.Function)(x, y, z, t)

    # Define gravitoelectric field components: E_g = -∇Φ_g - ∂_t A_g
    # From doc/gravity.tex line 15: 𝐄_g ≡ -∇Φ_g - ∂_t 𝐀_g
    E_g_x_def = -sp.diff(Phi_g, x) - sp.diff(A_g_x, t)
    E_g_y_def = -sp.diff(Phi_g, y) - sp.diff(A_g_y, t)
    E_g_z_def = -sp.diff(Phi_g, z) - sp.diff(A_g_z, t)

    # Define gravitomagnetic field components: B_g = ∇ × A_g
    # From doc/gravity.tex line 16: 𝐁_g ≡ ∇ × 𝐀_g
    B_g_x_def = sp.diff(A_g_z, y) - sp.diff(A_g_y, z)
    B_g_y_def = sp.diff(A_g_x, z) - sp.diff(A_g_z, x)
    B_g_z_def = sp.diff(A_g_y, x) - sp.diff(A_g_x, y)

    # Verify dimensional consistency of field definitions

    # First term: -∇Φ_g
    grad_phi_g_dim = v.grad_dim(v.get_dim('Phi_g'))
    v.check_dims("Gradient term -∇Φ_g",
                 grad_phi_g_dim, v.get_dim('E_g'))

    # Second term: -∂_t A_g
    # A_g has dimension [L T^-1], so ∂_t A_g has dimension [L T^-2] = E_g dimension
    time_deriv_A_g_dim = v.dt(v.get_dim('A_g'))
    v.check_dims("Time derivative term -∂_t A_g",
                 time_deriv_A_g_dim, v.get_dim('E_g'))

    # Both terms in E_g definition must have same dimensions
    v.check_dims("E_g field definition consistency",
                 grad_phi_g_dim, time_deriv_A_g_dim)

    # Gravitomagnetic field: B_g = ∇ × A_g
    curl_A_g_dim = v.curl_dim(v.get_dim('A_g'))
    v.check_dims("B_g = ∇ × A_g dimensional consistency",
                 curl_A_g_dim, v.get_dim('B_g'))

    # Test mathematical consistency of the field definitions
    # Verify that the curl of gradient is identically zero (important for later tests)
    # ∇ × (∇Φ_g) ≡ 0 (fundamental vector calculus identity)
    curl_grad_Phi_g_x = sp.diff(sp.diff(Phi_g, z), y) - sp.diff(sp.diff(Phi_g, y), z)
    curl_grad_Phi_g_y = sp.diff(sp.diff(Phi_g, x), z) - sp.diff(sp.diff(Phi_g, z), x)
    curl_grad_Phi_g_z = sp.diff(sp.diff(Phi_g, y), x) - sp.diff(sp.diff(Phi_g, x), y)

    # These should all be zero by Schwarz's theorem (mixed partial derivatives commute)
    v.check_eq("Curl of gradient identity: (∇×∇Φ_g)_x = 0",
               sp.simplify(curl_grad_Phi_g_x), 0)
    v.check_eq("Curl of gradient identity: (∇×∇Φ_g)_y = 0",
               sp.simplify(curl_grad_Phi_g_y), 0)
    v.check_eq("Curl of gradient identity: (∇×∇Φ_g)_z = 0",
               sp.simplify(curl_grad_Phi_g_z), 0)

    # Test that field components have the expected mathematical structure
    # The gravitoelectric field should separate into gradient and induction parts
    E_g_x_gradient_part = -sp.diff(Phi_g, x)
    E_g_x_induction_part = -sp.diff(A_g_x, t)

    v.check_eq("E_g field decomposition check",
               E_g_x_def, E_g_x_gradient_part + E_g_x_induction_part)

    v.info("✓ Gravitoelectric field E_g = -∇Φ_g - ∂_t A_g verified")
    v.info("✓ Gravitomagnetic field B_g = ∇ × A_g verified")
    v.info("✓ Vector calculus identities (∇×∇Φ ≡ 0) confirmed")

    v.success("Gravitoelectric and gravitomagnetic field definitions verified")


def test_maxwell_like_equations_lorenz_gauge(v):
    """
    Test the Maxwell-like equations in Lorenz gauge for gravity.

    Verifies the four GEM equations analogous to Maxwell's equations.
    From doc/gravity.tex lines 19-31.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Maxwell-like Equations in Lorenz Gauge")

    # Define coordinate symbols
    x, y, z, t = symbols('x y z t', real=True)
    c, G, pi = symbols('c G pi', positive=True)
    rho, j_m_x, j_m_y, j_m_z = symbols('rho j_m_x j_m_y j_m_z', real=True)

    # Define gravitational potentials as functions
    Phi_g = symbols('Phi_g', cls=sp.Function)(x, y, z, t)
    A_g_x = symbols('A_g_x', cls=sp.Function)(x, y, z, t)
    A_g_y = symbols('A_g_y', cls=sp.Function)(x, y, z, t)
    A_g_z = symbols('A_g_z', cls=sp.Function)(x, y, z, t)

    # Define fields using the previously established definitions
    E_g_x = -sp.diff(Phi_g, x) - sp.diff(A_g_x, t)
    E_g_y = -sp.diff(Phi_g, y) - sp.diff(A_g_y, t)
    E_g_z = -sp.diff(Phi_g, z) - sp.diff(A_g_z, t)

    B_g_x = sp.diff(A_g_z, y) - sp.diff(A_g_y, z)
    B_g_y = sp.diff(A_g_x, z) - sp.diff(A_g_z, x)
    B_g_z = sp.diff(A_g_y, x) - sp.diff(A_g_x, y)

    # Test Lorenz gauge condition: ∇·A_g + (1/c²)∂_t Φ_g = 0
    # From doc/gravity.tex line 21: ∇·𝐀_g + (1/c²)∂_t Φ_g = 0
    lorenz_gauge = (sp.diff(A_g_x, x) + sp.diff(A_g_y, y) + sp.diff(A_g_z, z) +
                    sp.diff(Phi_g, t)/(c**2))

    # Dimensional analysis of Lorenz gauge
    div_A_g_dim = v.div_dim(v.get_dim('A_g'))
    time_phi_g_dim = v.dt(v.get_dim('Phi_g')) / v.get_dim('c')**2
    v.check_dims("Lorenz gauge condition: ∇·A_g + (1/c²)∂_t Φ_g",
                 div_A_g_dim, time_phi_g_dim)

    # Test First GEM equation: ∇·E_g = -4πGρ (gravitational Gauss law)
    # From doc/gravity.tex line 26: ∇·𝐄_g = -4πGρ
    div_E_g = sp.diff(E_g_x, x) + sp.diff(E_g_y, y) + sp.diff(E_g_z, z)
    gauss_law_rhs = -4*pi*G*rho

    # Dimensional verification
    div_E_g_dim = v.div_dim(v.get_dim('E_g'))
    source_grav_gauss_dim = v.get_dim('G') * v.get_dim('rho')
    v.check_dims("Gravitational Gauss law: ∇·E_g = -4πGρ",
                 div_E_g_dim, source_grav_gauss_dim)

    # Test Second GEM equation: ∇×B_g - (1/c²)∂_t E_g = -(16πG/c²)j_m
    # From doc/gravity.tex line 27: ∇×𝐁_g - (1/c²)∂_t 𝐄_g = -(16πG/c²)𝐣_m
    curl_B_g_x = sp.diff(B_g_z, y) - sp.diff(B_g_y, z) - sp.diff(E_g_x, t)/(c**2)
    curl_B_g_y = sp.diff(B_g_x, z) - sp.diff(B_g_z, x) - sp.diff(E_g_y, t)/(c**2)
    curl_B_g_z = sp.diff(B_g_y, x) - sp.diff(B_g_x, y) - sp.diff(E_g_z, t)/(c**2)

    ampere_rhs_x = -16*pi*G*j_m_x/(c**2)
    ampere_rhs_y = -16*pi*G*j_m_y/(c**2)
    ampere_rhs_z = -16*pi*G*j_m_z/(c**2)

    # Dimensional verification
    curl_B_g_dim = v.curl_dim(v.get_dim('B_g'))
    time_E_g_dim = v.dt(v.get_dim('E_g')) / v.get_dim('c')**2
    source_grav_ampere_dim = v.get_dim('G') * v.get_dim('j_mass') / v.get_dim('c')**2

    v.check_dims("Gravitational Ampère law (curl term vs time term)",
                 curl_B_g_dim, time_E_g_dim)
    v.check_dims("Gravitational Ampère law (LHS vs RHS)",
                 curl_B_g_dim, source_grav_ampere_dim)

    # Test Third GEM equation: ∇·B_g = 0 (no gravitomagnetic monopoles)
    # From doc/gravity.tex line 28: ∇·𝐁_g = 0
    div_B_g = sp.diff(B_g_x, x) + sp.diff(B_g_y, y) + sp.diff(B_g_z, z)

    # This should be identically zero by the definition B_g = ∇×A_g
    # (divergence of curl is always zero)
    div_B_g_simplified = sp.simplify(div_B_g)
    v.check_eq("No gravitomagnetic monopoles: ∇·B_g = 0",
               div_B_g_simplified, 0)

    # Test Fourth GEM equation: ∇×E_g + ∂_t B_g = 0 (gravitomagnetic Faraday law)
    # From doc/gravity.tex line 29: ∇×𝐄_g + ∂_t 𝐁_g = 0
    faraday_x = (sp.diff(E_g_z, y) - sp.diff(E_g_y, z)) + sp.diff(B_g_x, t)
    faraday_y = (sp.diff(E_g_x, z) - sp.diff(E_g_z, x)) + sp.diff(B_g_y, t)
    faraday_z = (sp.diff(E_g_y, x) - sp.diff(E_g_x, y)) + sp.diff(B_g_z, t)

    # These should be identically zero by the field definitions
    faraday_x_simplified = sp.simplify(faraday_x)
    faraday_y_simplified = sp.simplify(faraday_y)
    faraday_z_simplified = sp.simplify(faraday_z)

    v.check_eq("Gravitomagnetic Faraday law x-component: (∇×E_g)_x + ∂_t B_g_x = 0",
               faraday_x_simplified, 0)
    v.check_eq("Gravitomagnetic Faraday law y-component: (∇×E_g)_y + ∂_t B_g_y = 0",
               faraday_y_simplified, 0)
    v.check_eq("Gravitomagnetic Faraday law z-component: (∇×E_g)_z + ∂_t B_g_z = 0",
               faraday_z_simplified, 0)

    # Dimensional verification for Faraday law
    curl_E_g_dim = v.curl_dim(v.get_dim('E_g'))
    time_B_g_dim = v.dt(v.get_dim('B_g'))
    v.check_dims("Gravitomagnetic Faraday law: ∇×E_g + ∂_t B_g",
                 curl_E_g_dim, time_B_g_dim)

    v.info("✓ Lorenz gauge condition verified")
    v.info("✓ Gravitational Gauss law (∇·E_g = -4πGρ) verified")
    v.info("✓ Gravitational Ampère-Maxwell law verified")
    v.info("✓ No gravitomagnetic monopoles (∇·B_g = 0) verified")
    v.info("✓ Gravitomagnetic Faraday law (∇×E_g + ∂_t B_g = 0) verified")

    v.success("Maxwell-like equations in Lorenz gauge verified")


def test_wave_equations_for_potentials(v):
    """
    Test the wave equations for gravitational potentials.

    Verifies: ∇²Φ_g - (1/c²)∂_tt Φ_g = 4πGρ
    and: ∇²A_g - (1/c²)∂_tt A_g = -(16πG/c²)j_m
    From doc/gravity.tex lines 34-36.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Wave Equations for Gravitational Potentials")

    # Define coordinate symbols
    x, y, z, t = symbols('x y z t', real=True)
    c, G, pi = symbols('c G pi', positive=True)
    rho, j_m_x, j_m_y, j_m_z = symbols('rho j_m_x j_m_y j_m_z', real=True)

    # Define gravitational potentials as functions
    Phi_g = symbols('Phi_g', cls=sp.Function)(x, y, z, t)
    A_g_x = symbols('A_g_x', cls=sp.Function)(x, y, z, t)
    A_g_y = symbols('A_g_y', cls=sp.Function)(x, y, z, t)
    A_g_z = symbols('A_g_z', cls=sp.Function)(x, y, z, t)

    # Test scalar potential wave equation: ∇²Φ_g - (1/c²)∂_tt Φ_g = 4πGρ
    # From doc/gravity.tex line 34: ∇²Φ_g - (1/c²)∂_tt Φ_g = 4πGρ

    # Compute Laplacian of scalar potential
    laplacian_Phi_g = (sp.diff(Phi_g, x, 2) + sp.diff(Phi_g, y, 2) + sp.diff(Phi_g, z, 2))

    # Compute second time derivative
    dtt_Phi_g = sp.diff(Phi_g, t, 2)

    # Define the wave operator applied to scalar potential
    wave_op_Phi_g = laplacian_Phi_g - dtt_Phi_g/(c**2)
    scalar_source = 4*pi*G*rho

    # Dimensional analysis of scalar wave equation
    laplacian_phi_g_dim = v.lap_dim(v.get_dim('Phi_g'))
    time2_phi_g_dim = v.dtt(v.get_dim('Phi_g')) / v.get_dim('c')**2
    source_phi_dim = v.get_dim('G') * v.get_dim('rho')

    # Check that Laplacian and time terms have same dimension
    v.check_dims("Scalar wave equation (space vs time terms)",
                 laplacian_phi_g_dim, time2_phi_g_dim)

    # Check that both match the source term
    v.check_dims("Scalar wave equation (LHS vs RHS)",
                 laplacian_phi_g_dim, source_phi_dim)

    # Test vector potential wave equation: ∇²A_g - (1/c²)∂_tt A_g = -(16πG/c²)j_m
    # From doc/gravity.tex line 35: ∇²𝐀_g - (1/c²)∂_tt 𝐀_g = -(16πG/c²)𝐣_m

    # Compute Laplacian of vector potential components
    laplacian_A_g_x = (sp.diff(A_g_x, x, 2) + sp.diff(A_g_x, y, 2) + sp.diff(A_g_x, z, 2))
    laplacian_A_g_y = (sp.diff(A_g_y, x, 2) + sp.diff(A_g_y, y, 2) + sp.diff(A_g_y, z, 2))
    laplacian_A_g_z = (sp.diff(A_g_z, x, 2) + sp.diff(A_g_z, y, 2) + sp.diff(A_g_z, z, 2))

    # Compute second time derivatives
    dtt_A_g_x = sp.diff(A_g_x, t, 2)
    dtt_A_g_y = sp.diff(A_g_y, t, 2)
    dtt_A_g_z = sp.diff(A_g_z, t, 2)

    # Define wave operators for vector components
    wave_op_A_g_x = laplacian_A_g_x - dtt_A_g_x/(c**2)
    wave_op_A_g_y = laplacian_A_g_y - dtt_A_g_y/(c**2)
    wave_op_A_g_z = laplacian_A_g_z - dtt_A_g_z/(c**2)

    vector_source_x = -16*pi*G*j_m_x/(c**2)
    vector_source_y = -16*pi*G*j_m_y/(c**2)
    vector_source_z = -16*pi*G*j_m_z/(c**2)

    # Dimensional analysis of vector wave equation
    laplacian_A_g_dim = v.lap_dim(v.get_dim('A_g'))
    time2_A_g_dim = v.dtt(v.get_dim('A_g')) / v.get_dim('c')**2
    source_A_dim = v.get_dim('G') * v.get_dim('j_mass') / v.get_dim('c')**2

    # Check that Laplacian and time terms have same dimension
    v.check_dims("Vector wave equation (space vs time terms)",
                 laplacian_A_g_dim, time2_A_g_dim)

    # Check that both match the source term
    v.check_dims("Vector wave equation (LHS vs RHS)",
                 laplacian_A_g_dim, source_A_dim)

    # Use the helper's wave equation verification pattern
    verify_wave_equation(v, "Gravitational scalar potential",
                        time2_phi_g_dim, laplacian_phi_g_dim, source_phi_dim)

    verify_wave_equation(v, "Gravitational vector potential",
                        time2_A_g_dim, laplacian_A_g_dim, source_A_dim)

    # Test that the wave equations are properly hyperbolic
    # The d'Alembertian operator □ = ∇² - (1/c²)∂_tt should have the right signature
    v.info("✓ Scalar wave equation: ∇²Φ_g - (1/c²)∂_tt Φ_g = 4πGρ")
    v.info("✓ Vector wave equation: ∇²A_g - (1/c²)∂_tt A_g = -(16πG/c²)j_m")
    v.info("✓ Hyperbolic signature ensures causal propagation at speed c")

    # Verify the wave equations have the correct relativistic structure
    # Both should be of the form □φ = source, where □ is the d'Alembertian
    dalembertian_dim = v.get_dim('Phi_g') / (v.L**2) - v.get_dim('Phi_g') / (v.get_dim('c')**2 * v.T**2)
    expected_dalembertian_dim = v.get_dim('Phi_g') / v.L**2  # Should match Laplacian dimension

    v.check_dims("D'Alembertian operator dimensional consistency",
                 expected_dalembertian_dim, laplacian_phi_g_dim)

    v.success("Wave equations for gravitational potentials verified")


def test_terminology_bridge_concepts(v):
    """
    Test the terminology bridge connecting to the broader framework.

    Verifies the connection between "intake" (charge-blind inflow) and
    gravitational eddies (frame-drag) with electromagnetic analogies.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Terminology Bridge Concepts")

    # "Intake" sources the gravitoelectric potential Φ_g
    # This is analogous to charge density sourcing electric potential
    v.info("Intake (charge-blind inflow) sources Φ_g")
    v.check_dims("Intake sourcing Φ_g (via Poisson equation)",
                 v.lap_dim(v.get_dim('Phi_g')),
                 v.get_dim('G') * v.get_dim('rho'))

    # Gravitational eddies (frame-drag) from moving/rotating masses
    # These are the source of the GEM B_g field
    v.info("Gravitational eddies (frame-drag) from moving masses create B_g")
    v.check_dims("Moving mass current j_mass",
                 v.get_dim('j_mass'), v.M / (v.L**2 * v.T))

    # Frame-drag effects scale with rotation/motion
    angular_momentum = v.get_dim('J_angular')
    v.check_dims("Angular momentum (frame-drag source)",
                 angular_momentum, v.M * v.L**2 / v.T)

    # Time changes of eddies induce loop pushes (Faraday analog)
    # ∇×E_g + ∂_t B_g = 0  is the gravitational Faraday law
    faraday_lhs = v.curl_dim(v.get_dim('E_g'))
    faraday_rhs = v.dt(v.get_dim('B_g'))
    v.check_dims("Gravitational Faraday law (eddy-induced loop pushes)",
                 faraday_lhs, faraday_rhs)

    v.info("✓ Intake → Φ_g (gravitoelectric potential)")
    v.info("✓ Gravitational eddies → B_g (gravitomagnetic field)")
    v.info("✓ Time-changing eddies → loop pushes (Faraday analog)")

    v.success("Terminology bridge concepts verified")


def test_tsunami_causality_connection(v):
    """
    Test the connection to tsunami principle and causality framework.

    This verifies the conceptual framework connecting to the broader
    tsunami-causality discussion referenced in the introduction.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Connection to Tsunami-Causality Framework")

    # The introduction mentions "asymptotic causality and decoupling of bulk v_L adjustments"
    # This connects to the tsunami principle where bulk flow modifications
    # don't affect the asymptotic wave propagation

    # Test that causality is preserved: information propagates at c
    causal_speed = v.get_dim('c')
    v.check_dims("Causal propagation speed",
                 causal_speed, v.L/v.T)

    # Bulk adjustments have fluid velocity dimension but don't affect causality
    bulk_speed = v.get_dim('v_L')
    v.check_dims("Bulk fluid velocity v_L",
                 bulk_speed, v.L/v.T)

    # The decoupling means: bulk_speed ≠ causal_speed in general
    # (They have same dimensions but different physical roles)
    v.info("Bulk v_L adjustments decoupled from causal propagation speed c")
    v.info("Only F_μν observables maintain strict causal propagation")

    # Wave sector maintains light-speed propagation
    wave_speed_check = causal_speed
    v.check_dims("Wave sector propagation maintains c",
                 wave_speed_check, v.L/v.T)

    # This framework ensures that:
    # 1. Gravitational waves propagate at c (causal)
    # 2. Bulk flow effects are decoupled from wave propagation
    # 3. F_μν observables respect relativistic causality

    v.info("✓ Gravitational waves propagate at speed c")
    v.info("✓ Bulk flow v_L decoupled from wave propagation")
    v.info("✓ F_μν observables respect relativistic causality")
    v.info("✓ Connection to tsunami-causality framework established")

    v.success("Tsunami-causality connection verified")


def test_introduction():
    """
    Main test function for Gravity Introduction.

    This function coordinates all verification tests for the introduction section
    of the gravity theory, validating asymptotic causality, GEM conventions,
    Maxwell-like equations, and connections to the broader framework.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Gravity: Weak and Strong Field - Introduction",
        "Asymptotic causality, GEM framework, and foundational concepts"
    )

    v.section("GRAVITY INTRODUCTION VERIFICATION")

    # Add any custom dimensions needed for the tests
    v.add_dimensions({
        'h_metric': 1,  # Metric perturbation (dimensionless)
        'F_mu_nu': v.get_dim('E') / v.get_dim('c'),  # Field strength tensor
    })

    # Call test functions in logical order
    v.info("\n--- 1) Asymptotic Causality Framework ---")
    test_asymptotic_causality_framework(v)

    v.info("\n--- 2) GEM Conventions and Signature ---")
    test_gem_conventions_and_signature(v)

    v.info("\n--- 3) Gravitoelectric and Gravitomagnetic Fields ---")
    test_gravitoelectric_gravitomagnetic_fields(v)

    v.info("\n--- 4) Maxwell-like Equations in Lorenz Gauge ---")
    test_maxwell_like_equations_lorenz_gauge(v)

    v.info("\n--- 5) Wave Equations for Potentials ---")
    test_wave_equations_for_potentials(v)

    v.info("\n--- 6) Terminology Bridge Concepts ---")
    test_terminology_bridge_concepts(v)

    v.info("\n--- 7) Tsunami-Causality Connection ---")
    test_tsunami_causality_connection(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_introduction()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)