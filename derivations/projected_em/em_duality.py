#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EM Duality - Verification
==========================================

Complete verification of electromagnetic duality properties in the projected
electromagnetic theory. This section tests the symmetry relationships between
electric and magnetic fields in Maxwell's equations and their implications
for field transformations.

Key concepts tested:
- Homogeneous Maxwell equations and their dimensional consistency
- Electric-magnetic field duality transformations
- Field tensor antisymmetry properties
- Dimensional analysis of electromagnetic field relationships
- Source-free field equation symmetries

Based on doc/projected_em.tex, homogeneous Maxwell equations and field definitions.
"""

import os
import sys
import sympy as sp
from sympy import symbols, simplify, diff, integrate, limit, oo, pi, sqrt

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    batch_check_dims,
    quick_verify,
    verify_conservation_law,
    verify_poisson_equation
)

# Initialize verification helper
v = PhysicsVerificationHelper(
    "EM Duality",
    "Verification of electromagnetic duality in projected field theory"
)


def test_homogeneous_maxwell_dimensions(v):
    """
    Test dimensional consistency of the homogeneous Maxwell equations.

    From doc/projected_em.tex:
    ∇·B = 0
    ∇×E + ∂_t B = 0

    These are exact identities in the projected theory.
    """
    v.section("Homogeneous Maxwell Equations - Dimensional Analysis")

    # Test ∇·B = 0 dimensions
    v.check_dims("Divergence of magnetic field",
                 v.div_dim(v.dims['B']),
                 v.dims['B'] / v.L)

    # Test Faraday's law: ∇×E + ∂_t B = 0 dimensions
    curl_E_dim = v.curl_dim(v.dims['E'])
    dB_dt_dim = v.dt(v.dims['B'])

    v.check_dims("Curl of electric field",
                 curl_E_dim,
                 v.dims['B'] / v.T)

    v.check_dims("Time derivative of magnetic field",
                 dB_dt_dim,
                 v.dims['B'] / v.T)

    # Test that both terms in Faraday's law have same dimensions
    v.check_dims("Faraday's law dimensional consistency",
                 curl_E_dim,
                 dB_dt_dim)


def test_field_definition_dimensions(v):
    """
    Test the fundamental field definitions from potentials:
    B = ∇×A
    E = -∂_t A - ∇Φ
    """
    v.section("Field Definitions - Dimensional Analysis")

    # Test B = ∇×A dimensions
    v.check_dims("B from curl of A",
                 v.curl_dim(v.dims['A']),
                 v.dims['B'])

    # Test E = -∂_t A - ∇Φ dimensions
    induction_term = v.dt(v.dims['A'])
    potential_term = v.grad_dim(v.dims['Phi'])

    v.check_dims("E induction term: ∂_t A",
                 induction_term,
                 v.dims['E'])

    v.check_dims("E potential term: ∇Φ",
                 potential_term,
                 v.dims['E'])

    # Test that both terms of E have consistent dimensions
    v.check_dims("E field components dimensional consistency",
                 induction_term,
                 potential_term)


def test_electromagnetic_duality_dimensions(v):
    """
    Test dimensional consistency of electromagnetic duality transformations.
    In vacuum, the duality E ↔ cB involves impedance relationships.
    """
    v.section("Electromagnetic Duality - Dimensional Analysis")

    # Test E ↔ cB duality dimensions
    E_to_B_transform = v.dims['E'] / v.dims['c']
    B_to_E_transform = v.dims['B'] * v.dims['c']

    v.check_dims("E/c has magnetic field dimensions",
                 E_to_B_transform,
                 v.dims['B'])

    v.check_dims("cB has electric field dimensions",
                 B_to_E_transform,
                 v.dims['E'])

    # Test vacuum impedance relationship: Z₀ = √(μ₀/ε₀) = μ₀c = 1/(ε₀c)
    impedance_from_mu = v.dims['mu_0'] * v.dims['c']
    impedance_from_eps = 1 / (v.dims['epsilon_0'] * v.dims['c'])

    v.check_dims("Vacuum impedance from μ₀c",
                 impedance_from_mu,
                 v.dims['Z_0'])

    v.check_dims("Vacuum impedance from 1/(ε₀c)",
                 impedance_from_eps,
                 v.dims['Z_0'])

    # Test E/H and B/D ratios for duality
    E_over_H = v.dims['E'] / v.dims['H']
    B_over_D = v.dims['B'] / v.dims['D']

    v.check_dims("E/H ratio has impedance dimensions",
                 E_over_H,
                 v.dims['Z_0'])

    v.check_dims("B/D ratio has impedance dimensions",
                 B_over_D,
                 v.dims['Z_0'])


def test_field_tensor_dimensions(v):
    """
    Test electromagnetic field tensor dimensional properties.
    F_μν components correspond to E and B field components.
    """
    v.section("Field Tensor - Dimensional Analysis")

    # Field tensor components in 4D
    # F₀ᵢ components are electric field related
    # Fᵢⱼ components are magnetic field related

    v.check_dims("Field tensor temporal components (E-like)",
                 v.dt(v.dims['A']),  # ∂₀Aᵢ ~ ∂_t A ~ E
                 v.dims['E'])

    v.check_dims("Field tensor spatial components (B-like)",
                 v.curl_dim(v.dims['A']),  # ∂ᵢAⱼ - ∂ⱼAᵢ ~ B
                 v.dims['B'])

    # Test 4-potential consistency
    v.check_dims("4-potential temporal component A₀",
                 v.dims['A0'],
                 v.dims['A'])

    v.check_dims("4-potential spatial components Aᵢ",
                 v.dims['A1'],
                 v.dims['A'])


def test_lorentz_invariant_combinations(v):
    """
    Test Lorentz invariant combinations of electromagnetic fields.
    These are preserved under electromagnetic duality transformations.
    """
    v.section("Lorentz Invariants - Dimensional Analysis")

    # First invariant: E·B (pseudoscalar)
    E_dot_B = v.dims['E'] * v.dims['B']
    v.check_dims("E·B invariant",
                 E_dot_B,
                 v.dims['E'] * v.dims['B'])

    # Second invariant: B² - E²/c² (scalar)
    # In natural units where c=1, this is B² - E²
    B_squared = v.dims['B']**2
    E_squared_over_c_squared = (v.dims['E']**2) / (v.dims['c']**2)

    v.check_dims("B² term in invariant",
                 B_squared,
                 v.dims['B']**2)

    v.check_dims("E²/c² term in invariant",
                 E_squared_over_c_squared,
                 v.dims['B']**2)  # Should match B² dimensions

    # Energy density related invariants
    electric_energy_density = v.dims['epsilon_0'] * (v.dims['E']**2)
    magnetic_energy_density = (v.dims['B']**2) / v.dims['mu_0']

    v.check_dims("Electric energy density",
                 electric_energy_density,
                 v.dims['u_EM'])

    v.check_dims("Magnetic energy density",
                 magnetic_energy_density,
                 v.dims['u_EM'])

    # Test that both energy densities have same dimensions (required for duality)
    v.check_dims("EM energy density consistency",
                 electric_energy_density,
                 magnetic_energy_density)


def test_poynting_vector_duality(v):
    """
    Test Poynting vector and its role in electromagnetic duality.
    S = (1/μ₀) E × B represents energy flow.
    """
    v.section("Poynting Vector - Dimensional Analysis")

    # Poynting vector: S = (1/μ₀) E × B
    poynting_dim = (v.dims['E'] * v.dims['B']) / v.dims['mu_0']

    v.check_dims("Poynting vector S = E×B/μ₀",
                 poynting_dim,
                 v.dims['S_poynting'])

    # Alternative form using impedance: S = E²/Z₀ = cε₀E²
    poynting_from_E = v.dims['c'] * v.dims['epsilon_0'] * (v.dims['E']**2)
    poynting_from_B = v.dims['c'] * (v.dims['B']**2) / v.dims['mu_0']

    v.check_dims("Poynting vector from E field",
                 poynting_from_E,
                 v.dims['S_poynting'])

    v.check_dims("Poynting vector from B field",
                 poynting_from_B,
                 v.dims['S_poynting'])

    # Test consistency of all Poynting vector forms
    v.check_dims("Poynting vector consistency (E×B vs E² forms)",
                 poynting_dim,
                 poynting_from_E)

    v.check_dims("Poynting vector consistency (E×B vs B² forms)",
                 poynting_dim,
                 poynting_from_B)


def test_maxwell_wave_equations(v):
    """
    Test wave equations derived from Maxwell's equations.
    From taking curl of Ampère-Maxwell and using Faraday's law.
    """
    v.section("Maxwell Wave Equations")

    # Define wave equation operator dimensions
    t, x = define_symbols_batch(['t', 'x'], real=True)

    # Wave operator: ∇² - (1/c²)∂²/∂t²
    laplacian_dim = v.dims['laplacian']  # ∇²
    time_second_deriv_coeff = 1 / (v.dims['c']**2)  # 1/c²

    v.check_dims("Laplacian operator",
                 laplacian_dim,
                 v.L**(-2))

    v.check_dims("Time derivative coefficient 1/c²",
                 time_second_deriv_coeff,
                 v.T**2 / v.L**2)

    # Test that both terms in wave equation have same dimensions when applied to fields
    wave_spatial_term = laplacian_dim * v.dims['E']
    wave_temporal_term = (time_second_deriv_coeff * v.dims['E']) / (v.T**2)

    v.check_dims("Spatial term ∇²E in wave equation",
                 wave_spatial_term,
                 v.dims['E'] / v.L**2)

    v.check_dims("Temporal term (1/c²)∂²E/∂t² in wave equation",
                 wave_temporal_term,
                 v.dims['E'] / v.L**2)

    # Test dimensional consistency of wave equation
    v.check_dims("Wave equation dimensional consistency",
                 wave_spatial_term,
                 wave_temporal_term)


def test_electromagnetic_field_relationships(v):
    """
    Test fundamental relationships between E and B fields in different contexts.
    """
    v.section("Electromagnetic Field Relationships")

    # Test plane wave relationships in vacuum
    # For plane waves: |E| = c|B| (in SI units)
    plane_wave_E_from_B = v.dims['c'] * v.dims['B']

    v.check_dims("Plane wave E from B: |E| = c|B|",
                 plane_wave_E_from_B,
                 v.dims['E'])

    # Test electromagnetic energy equipartition in vacuum
    # Electric and magnetic energy densities are equal for plane waves
    electric_energy = v.dims['epsilon_0'] * v.dims['E']**2
    magnetic_energy = v.dims['B']**2 / v.dims['mu_0']

    v.check_dims("Electromagnetic energy equipartition",
                 electric_energy,
                 magnetic_energy)

    # Test momentum density of electromagnetic field
    # p_EM = S/c² = ε₀ E × B
    momentum_density_from_poynting = v.dims['S_poynting'] / (v.dims['c']**2)
    momentum_density_from_fields = v.dims['epsilon_0'] * v.dims['E'] * v.dims['B']

    v.check_dims("EM momentum density from Poynting",
                 momentum_density_from_poynting,
                 v.M / (v.L**2 * v.T))

    v.check_dims("EM momentum density from fields",
                 momentum_density_from_fields,
                 v.M / (v.L**2 * v.T))

    v.check_dims("EM momentum density consistency",
                 momentum_density_from_poynting,
                 momentum_density_from_fields)


def test_em_duality():
    """Test runner compatible function for EM duality."""
    return run_all_tests()


def run_all_tests():
    """Run all electromagnetic duality tests."""
    test_homogeneous_maxwell_dimensions(v)
    test_field_definition_dimensions(v)
    test_electromagnetic_duality_dimensions(v)
    test_field_tensor_dimensions(v)
    test_lorentz_invariant_combinations(v)
    test_poynting_vector_duality(v)
    test_maxwell_wave_equations(v)
    test_electromagnetic_field_relationships(v)

    # Print summary and return success rate
    return v.summary()


if __name__ == "__main__":
    run_all_tests()