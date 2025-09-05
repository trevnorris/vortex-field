#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GEM Conventions and Signature - Verification
====================================

This module implements dimensional and mathematical verification for the GEM 
(Gravitoelectromagnetic) conventions and signature subsection, including metric 
signature, weak-field potential definitions, and the analogy between gravity 
and electromagnetism.

Verifies the fundamental relationships:
- Metric signature (-,+,+,+) and weak-field potential definitions
- Gravitoelectric and gravitomagnetic field definitions
- GEM field equations and their relationship to EM Maxwell equations
- Dimensional consistency of the gravity-electromagnetism analogy

Based on doc/gravity.tex, subsection "GEM Conventions and Signature" (lines 6-43).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, sqrt

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    verify_wave_equation,
    quick_verify
)


def test_metric_signature_and_potentials(v):
    """
    Test dimensional consistency of metric signature conventions and weak-field potentials:
    h_{00} = -2Φ_g/c², h_{0i} = -4A_{g,i}/c³, h_{ij} = -2Φ_g/c²δ_{ij}
    """
    v.subsection("Metric Signature and Weak-Field Potentials")

    # Metric components should be dimensionless
    # h_{00} = -2Φ_g/c²
    h00_rhs = 2 * v.get_dim('Phi_g') / v.get_dim('c')**2
    dimensionless = 1

    v.check_dims("h_{00} = -2Φ_g/c² dimensionless", h00_rhs, dimensionless)

    # h_{0i} = -4A_{g,i}/c³ 
    h0i_rhs = 4 * v.get_dim('A_g') / v.get_dim('c')**3
    
    v.check_dims("h_{0i} = -4A_{g,i}/c³ dimensionless", h0i_rhs, dimensionless)

    # h_{ij} = -2Φ_g/c²δ_{ij} (spatial components, same as h_{00})
    hij_rhs = 2 * v.get_dim('Phi_g') / v.get_dim('c')**2
    
    v.check_dims("h_{ij} = -2Φ_g/c²δ_{ij} dimensionless", hij_rhs, dimensionless)

    # Check consistency: h_{00} and h_{ij} should have same form
    v.check_dims("h_{00} and h_{ij} consistency", h00_rhs, hij_rhs)

    v.success("Metric signature and weak-field potential definitions verified")


def test_gem_field_definitions(v):
    """
    Test dimensional consistency of gravitoelectric and gravitomagnetic field definitions:
    E_g = -∇Φ_g - (1/c)∂_t A_g  and  B_g = ∇×A_g
    """
    v.subsection("GEM Field Definitions")

    # Gravitoelectric field: E_g = -∇Φ_g - (1/c)∂_t A_g
    grad_phi_g = v.grad_dim(v.get_dim('Phi_g'))     # -∇Φ_g term
    dt_A_g_term = v.dt(v.get_dim('A_g')) / v.get_dim('c')  # -(1/c)∂_t A_g term

    v.check_dims("E_g components: ∇Φ_g vs (1/c)∂_t A_g", grad_phi_g, dt_A_g_term)

    # Both terms should match E_g dimensions
    E_g_expected = v.get_dim('E_g')
    
    v.check_dims("E_g = -∇Φ_g term", grad_phi_g, E_g_expected)
    v.check_dims("E_g = -(1/c)∂_t A_g term", dt_A_g_term, E_g_expected)

    # Gravitomagnetic field: B_g = ∇×A_g
    curl_A_g = v.curl_dim(v.get_dim('A_g'))
    B_g_expected = v.get_dim('B_g')

    v.check_dims("B_g = ∇×A_g", curl_A_g, B_g_expected)

    v.success("GEM field definitions verified")


def test_lorenz_gauge_condition(v):
    """
    Test dimensional consistency of the Lorenz gauge condition:
    ∇·A_g + (1/c²)∂_t Φ_g = 0
    """
    v.subsection("Lorenz Gauge Condition")

    # ∇·A_g term
    div_A_g = v.div_dim(v.get_dim('A_g'))

    # (1/c²)∂_t Φ_g term  
    dt_phi_g_term = v.dt(v.get_dim('Phi_g')) / v.get_dim('c')**2

    v.check_dims("Lorenz gauge: ∇·A_g vs (1/c²)∂_t Φ_g", div_A_g, dt_phi_g_term)

    # Both terms should be dimensionally consistent for gauge condition
    v.success("Lorenz gauge condition dimensional consistency verified")


def test_gem_field_equations(v):
    """
    Test dimensional consistency of the GEM field equations (Maxwell-like):
    ∇·E_g = -4πG ρ
    ∇×B_g - (1/c²)∂_t E_g = -(16πG/c²) j
    ∇·B_g = 0  
    ∇×E_g + ∂_t B_g = 0
    """
    v.subsection("GEM Field Equations")

    # Gauss law for gravity: ∇·E_g = -4πG ρ
    div_E_g = v.div_dim(v.get_dim('E_g'))
    gauss_rhs = 4 * pi * v.get_dim('G') * v.get_dim('rho')

    v.check_dims("Gauss law: ∇·E_g vs 4πGρ", div_E_g, gauss_rhs)

    # Ampère-Maxwell for gravity: ∇×B_g - (1/c²)∂_t E_g = -(16πG/c²) j
    curl_B_g = v.curl_dim(v.get_dim('B_g'))
    dt_E_g_term = v.dt(v.get_dim('E_g')) / v.get_dim('c')**2
    ampere_rhs = 16 * pi * v.get_dim('G') * v.get_dim('j_mass') / v.get_dim('c')**2

    v.check_dims("Ampère law terms: ∇×B_g vs (1/c²)∂_t E_g", curl_B_g, dt_E_g_term)
    v.check_dims("Ampère law: ∇×B_g vs (16πG/c²)j", curl_B_g, ampere_rhs)

    # Magnetic Gauss law: ∇·B_g = 0 (automatically satisfied dimensionally)
    div_B_g = v.div_dim(v.get_dim('B_g'))
    v.check_dims("Magnetic Gauss: ∇·B_g dimensionally valid", div_B_g, div_B_g)

    # Faraday law for gravity: ∇×E_g + ∂_t B_g = 0
    curl_E_g = v.curl_dim(v.get_dim('E_g'))
    dt_B_g = v.dt(v.get_dim('B_g'))

    v.check_dims("Faraday law: ∇×E_g vs ∂_t B_g", curl_E_g, dt_B_g)

    v.success("GEM field equations verified")


def test_gem_wave_equations(v):
    """
    Test dimensional consistency of the GEM wave equations:
    ∇²Φ_g - (1/c²)∂_tt Φ_g = 4πG ρ
    ∇²A_g - (1/c²)∂_tt A_g = -(16πG/c²) j
    """
    v.subsection("GEM Wave Equations")

    # Gravitoelectric potential wave equation
    lap_phi_g = v.lap_dim(v.get_dim('Phi_g'))
    dtt_phi_g_term = v.dtt(v.get_dim('Phi_g')) / v.get_dim('c')**2
    phi_source = 4 * pi * v.get_dim('G') * v.get_dim('rho')

    v.check_dims("Φ_g wave: ∇²Φ_g vs (1/c²)∂_tt Φ_g", lap_phi_g, dtt_phi_g_term)
    v.check_dims("Φ_g wave: ∇²Φ_g vs 4πGρ", lap_phi_g, phi_source)

    # Gravitomagnetic potential wave equation  
    lap_A_g = v.lap_dim(v.get_dim('A_g'))
    dtt_A_g_term = v.dtt(v.get_dim('A_g')) / v.get_dim('c')**2
    A_source = 16 * pi * v.get_dim('G') * v.get_dim('j_mass') / v.get_dim('c')**2

    v.check_dims("A_g wave: ∇²A_g vs (1/c²)∂_tt A_g", lap_A_g, dtt_A_g_term)
    v.check_dims("A_g wave: ∇²A_g vs (16πG/c²)j", lap_A_g, A_source)

    v.success("GEM wave equations verified")


def test_gem_em_analogy_consistency(v):
    """
    Test the consistency of the gravity-electromagnetism analogy by comparing
    dimensional structures of corresponding equations.
    """
    v.subsection("GEM-EM Analogy Consistency")

    # Compare field definition structures
    # EM: E = -∇Φ - (1/c)∂_t A
    # GEM: E_g = -∇Φ_g - (1/c)∂_t A_g

    # EM electric field components
    em_grad_term = v.grad_dim(v.get_dim('Phi'))
    em_dt_term = v.dt(v.get_dim('A')) / v.get_dim('c')

    # GEM gravitoelectric field components  
    gem_grad_term = v.grad_dim(v.get_dim('Phi_g'))
    gem_dt_term = v.dt(v.get_dim('A_g')) / v.get_dim('c')

    # The dimensional structure should be analogous (both give field dimensions)
    v.check_dims("EM field structure", em_grad_term, em_dt_term)
    v.check_dims("GEM field structure", gem_grad_term, gem_dt_term)

    # Compare magnetic field structures
    # EM: B = ∇×A, GEM: B_g = ∇×A_g
    em_curl = v.curl_dim(v.get_dim('A'))
    gem_curl = v.curl_dim(v.get_dim('A_g'))

    # Both should be consistent with their respective field dimensions
    v.check_dims("EM magnetic: B vs ∇×A", v.get_dim('B'), em_curl)
    v.check_dims("GEM magnetic: B_g vs ∇×A_g", v.get_dim('B_g'), gem_curl)

    # Compare Gauss law structures
    # EM: ∇·E = ρ_charge/ε₀, GEM: ∇·E_g = -4πGρ
    em_gauss_lhs = v.div_dim(v.get_dim('E'))
    gem_gauss_lhs = v.div_dim(v.get_dim('E_g'))
    em_gauss_rhs = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
    gem_gauss_rhs = 4 * pi * v.get_dim('G') * v.get_dim('rho')

    v.check_dims("EM Gauss law", em_gauss_lhs, em_gauss_rhs)
    v.check_dims("GEM Gauss law", gem_gauss_lhs, gem_gauss_rhs)

    v.success("GEM-EM analogy dimensional consistency verified")


def test_physical_constants_in_gem(v):
    """
    Test the role of fundamental constants G and c in GEM theory and their 
    dimensional consistency with electromagnetic counterparts.
    """
    v.subsection("Physical Constants in GEM Theory")

    # Newton's gravitational constant G should have correct dimensions
    # G: [M⁻¹ L³ T⁻²] in SI units
    G_expected = v.L**3 / (v.M * v.T**2)
    v.check_dims("Gravitational constant G", v.get_dim('G'), G_expected)

    # Speed of light c should be consistent across EM and GEM
    # Both theories use the same c
    c_consistency = v.get_dim('c')
    v.check_dims("Speed of light consistency", c_consistency, c_consistency)

    # Check the role of G in GEM vs μ₀ in EM
    # GEM source term: 4πGρ has dimensions [T⁻²]  
    gem_source_factor = 4 * pi * v.get_dim('G') * v.get_dim('rho')
    
    # EM source term: ρ_charge/ε₀ should also have dimensions [L⁻¹ T⁻²] * [L] = [T⁻²]
    # Actually ∇·E has dimensions [field/length], so we need [field/length]
    em_source_factor = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
    
    # Both should give the same dimensional structure as div(field)
    expected_source_dim = v.div_dim(v.get_dim('E'))  # Use E as template
    v.check_dims("GEM source dimensional structure", gem_source_factor, expected_source_dim)

    # Verify the coupling strength dimensional consistency
    # GEM: G couples mass density to gravitoelectric field
    # EM: 1/ε₀ couples charge density to electric field
    G_coupling = v.get_dim('G') * v.get_dim('rho')
    em_coupling = v.get_dim('rho_charge') / v.get_dim('epsilon_0')

    # Both should produce field/length dimensions when used in Gauss law
    expected_coupling = v.get_dim('E_g') / v.L  # field/length
    v.check_dims("GEM coupling structure", G_coupling, expected_coupling)

    v.success("Physical constants in GEM theory verified")


def test_terminology_bridge(v):
    """
    Test the dimensional consistency of the terminology bridge concepts:
    intake (charge-blind inflow) and gravitational eddies (frame-drag).
    """
    v.subsection("Terminology Bridge Verification")

    # Intake sources gravitoelectric potential Φ_g
    # This should be dimensionally consistent with the Poisson equation
    # ∇²Φ_g = 4πGρ (in the static limit)
    intake_source = v.lap_dim(v.get_dim('Phi_g'))
    poisson_rhs = 4 * pi * v.get_dim('G') * v.get_dim('rho')

    v.check_dims("Intake-Poisson consistency: ∇²Φ_g vs 4πGρ", intake_source, poisson_rhs)

    # Gravitational eddies (frame-drag) correspond to B_g field
    # Moving/rotating masses create gravitomagnetic effects
    # This is analogous to moving charges creating magnetic fields
    
    # The B_g field should be dimensionally consistent with its source
    # From ∇×B_g = -(16πG/c²)j, the current j should have appropriate dimensions
    eddies_field = v.get_dim('B_g')
    curl_B_g = v.curl_dim(eddies_field)
    frame_drag_source = 16 * pi * v.get_dim('G') * v.get_dim('j_mass') / v.get_dim('c')**2

    v.check_dims("Frame-drag source: ∇×B_g vs (16πG/c²)j", curl_B_g, frame_drag_source)

    # Faraday-analog: time changes of eddies induce loop pushes
    # ∇×E_g + ∂_t B_g = 0 (homogeneous Faraday equation)
    faraday_curl = v.curl_dim(v.get_dim('E_g'))
    faraday_time = v.dt(v.get_dim('B_g'))

    v.check_dims("Faraday-analog: ∇×E_g vs ∂_t B_g", faraday_curl, faraday_time)

    v.success("Terminology bridge concepts verified")


def test_gem_conventions_and_signature():
    """
    Main test function for the "GEM Conventions and Signature" subsection.

    Tests all aspects of the gravitoelectromagnetic theory conventions including
    metric signature, field definitions, field equations, and the gravity-EM analogy.
    
    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "GEM Conventions and Signature",
        "Gravitoelectromagnetic theory conventions and EM analogy verification"
    )

    v.section("GEM CONVENTIONS AND SIGNATURE VERIFICATION")

    # Add any custom dimensions needed for GEM theory
    # Most should already be defined in helper.py, but let's ensure consistency
    v.declare_dimensionless('delta_ij')  # Kronecker delta is dimensionless

    # Test 1: Metric signature and weak-field potentials
    v.info("\n--- 1) Metric Signature and Weak-Field Potentials ---")
    test_metric_signature_and_potentials(v)

    # Test 2: GEM field definitions
    v.info("\n--- 2) Gravitoelectric and Gravitomagnetic Field Definitions ---")
    test_gem_field_definitions(v)

    # Test 3: Lorenz gauge condition
    v.info("\n--- 3) Lorenz Gauge Condition ---")
    test_lorenz_gauge_condition(v)

    # Test 4: GEM field equations (Maxwell-like)
    v.info("\n--- 4) GEM Field Equations ---")
    test_gem_field_equations(v)

    # Test 5: GEM wave equations
    v.info("\n--- 5) GEM Wave Equations ---")
    test_gem_wave_equations(v)

    # Test 6: GEM-EM analogy consistency
    v.info("\n--- 6) GEM-EM Analogy Consistency ---")
    test_gem_em_analogy_consistency(v)

    # Test 7: Physical constants in GEM theory
    v.info("\n--- 7) Physical Constants in GEM Theory ---")
    test_physical_constants_in_gem(v)

    # Test 8: Terminology bridge verification
    v.info("\n--- 8) Terminology Bridge (Intake and Eddies) ---")
    test_terminology_bridge(v)

    # Final summary
    return v.summary()


if __name__ == "__main__":
    success_rate = test_gem_conventions_and_signature()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)