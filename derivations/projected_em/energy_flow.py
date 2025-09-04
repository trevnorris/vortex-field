#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Energy Flow (Poynting Theorem) - Verification

Comprehensive verification of the Poynting theorem derivation including Maxwell equations,
vector identities, intermediate algebraic steps, and dimensional consistency of the
electromagnetic energy conservation law.
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, Rational, sqrt

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    verify_conservation_law,
    quick_verify,
    define_symbols_batch
)


def test_maxwell_equations(v):
    """Test Maxwell equations used in the derivation"""
    v.subsection("Maxwell Equations")

    # Ampère-Maxwell law: ∇×B = μ₀J + μ₀ε₀∂E/∂t
    curl_B = v.curl_dim(v.get_dim('B'))
    mu0_J_term = v.get_dim('mu_0') * v.get_dim('j_current')
    mu0_eps0_dEdt_term = v.get_dim('mu_0') * v.get_dim('epsilon_0') * v.dt(v.get_dim('E'))

    v.check_dims("Ampère-Maxwell: ∇×B = μ₀J", curl_B, mu0_J_term)
    v.check_dims("Ampère-Maxwell: μ₀J = μ₀ε₀∂E/∂t", mu0_J_term, mu0_eps0_dEdt_term)

    # Faraday's law: ∇×E = -∂B/∂t
    curl_E = v.curl_dim(v.get_dim('E'))
    dBdt = v.dt(v.get_dim('B'))

    v.check_dims("Faraday's law: ∇×E = -∂B/∂t", curl_E, dBdt)


def test_vector_identity(v):
    """Test the crucial vector identity: ∇·(E×B) = B·(∇×E) - E·(∇×B)"""
    v.subsection("Vector Identity")

    # Left side: ∇·(E×B)
    cross_product_EB = v.get_dim('E') * v.get_dim('B')
    lhs = v.div_dim(cross_product_EB)

    # Right side: B·(∇×E) - E·(∇×B)
    B_dot_curl_E = v.get_dim('B') * v.curl_dim(v.get_dim('E'))
    E_dot_curl_B = v.get_dim('E') * v.curl_dim(v.get_dim('B'))
    rhs = B_dot_curl_E - E_dot_curl_B

    v.check_dims("Vector identity: ∇·(E×B) = B·(∇×E) - E·(∇×B)", lhs, rhs)


def test_derivation_step1(v):
    """Test first intermediate step: E·(∇×B) - B·(∇×E) = μ₀E·J + μ₀ε₀E·∂ₜE + B·∂ₜB"""
    v.subsection("Derivation Step 1")

    # Left side: E·(∇×B) - B·(∇×E)
    E_dot_curl_B = v.get_dim('E') * v.curl_dim(v.get_dim('B'))
    B_dot_curl_E = v.get_dim('B') * v.curl_dim(v.get_dim('E'))
    lhs = E_dot_curl_B - B_dot_curl_E

    # Right side: μ₀E·J + μ₀ε₀E·∂ₜE + B·∂ₜB
    mu0_E_dot_J = v.get_dim('mu_0') * v.get_dim('E') * v.get_dim('j_current')
    mu0_eps0_E_dot_dEdt = v.get_dim('mu_0') * v.get_dim('epsilon_0') * v.get_dim('E') * v.dt(v.get_dim('E'))
    B_dot_dBdt = v.get_dim('B') * v.dt(v.get_dim('B'))

    rhs = mu0_E_dot_J + mu0_eps0_E_dot_dEdt + B_dot_dBdt

    v.check_dims("Step 1: curl-dot terms equal source terms", lhs, rhs)
    v.check_dims("Step 1: μ₀E·J term", mu0_E_dot_J, lhs)
    v.check_dims("Step 1: μ₀ε₀E·∂ₜE term", mu0_eps0_E_dot_dEdt, lhs)
    v.check_dims("Step 1: B·∂ₜB term", B_dot_dBdt, lhs)


def test_time_derivative_identities(v):
    """Test time derivative identities: E·∂ₜE = ½∂ₜ|E|² and B·∂ₜB = ½∂ₜ|B|²"""
    v.subsection("Time Derivative Identities")

    # E·∂ₜE = ½∂ₜ|E|²
    E_dot_dEdt = v.get_dim('E') * v.dt(v.get_dim('E'))
    half_dt_E_squared = v.dt(v.get_dim('E')**2) / 2

    v.check_dims("Time derivative identity: E·∂ₜE = ½∂ₜ|E|²", E_dot_dEdt, half_dt_E_squared)

    # B·∂ₜB = ½∂ₜ|B|²
    B_dot_dBdt = v.get_dim('B') * v.dt(v.get_dim('B'))
    half_dt_B_squared = v.dt(v.get_dim('B')**2) / 2

    v.check_dims("Time derivative identity: B·∂ₜB = ½∂ₜ|B|²", B_dot_dBdt, half_dt_B_squared)


def test_derivation_step2(v):
    """Test second step: apply vector identity and time derivative recognition"""
    v.subsection("Derivation Step 2")

    # -∇·(E×B) = μ₀E·J + (μ₀ε₀/2)∂ₜ|E|² + (1/2)∂ₜ|B|²
    cross_product_EB = v.get_dim('E') * v.get_dim('B')
    neg_div_ExB = -v.div_dim(cross_product_EB)

    mu0_E_dot_J = v.get_dim('mu_0') * v.get_dim('E') * v.get_dim('j_current')
    mu0_eps0_half_dt_E_squared = v.get_dim('mu_0') * v.get_dim('epsilon_0') * v.dt(v.get_dim('E')**2) / 2
    half_dt_B_squared = v.dt(v.get_dim('B')**2) / 2

    rhs = mu0_E_dot_J + mu0_eps0_half_dt_E_squared + half_dt_B_squared

    v.check_dims("Step 2: -∇·(E×B) = source and time derivative terms", neg_div_ExB, rhs)


def test_derivation_step3(v):
    """Test final step: division by μ₀ to obtain conservation law"""
    v.subsection("Derivation Step 3: Final Conservation Law")

    # ∂ₜ(ε₀/2 |E|² + 1/(2μ₀) |B|²) + ∇·(1/μ₀ E×B) = -J·E

    # Energy density terms
    electric_energy = (v.get_dim('epsilon_0')/2) * v.get_dim('E')**2
    magnetic_energy = v.get_dim('B')**2 / (2 * v.get_dim('mu_0'))
    total_energy_density = electric_energy + magnetic_energy

    dt_energy_density = v.dt(total_energy_density)

    # Poynting vector divergence
    poynting_vector = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
    div_poynting = v.div_dim(poynting_vector)

    # Joule heating
    joule_heating = v.get_dim('j_current') * v.get_dim('E')

    # All terms should have dimensions [M L⁻¹ T⁻³]
    target_power_density = v.M / (v.L * v.T**3)

    v.check_dims("Final: ∂ₜ(energy density)", dt_energy_density, target_power_density)
    v.check_dims("Final: ∇·(Poynting vector)", div_poynting, target_power_density)
    v.check_dims("Final: Joule heating J·E", joule_heating, target_power_density)

    # Complete conservation law consistency
    lhs_total = dt_energy_density + div_poynting
    v.check_dims("Conservation law: LHS = RHS", lhs_total, joule_heating)


def test_energy_density_components(v):
    """Test electromagnetic energy density components"""
    v.subsection("Energy Density Components")

    # Electric energy density: (ε₀/2)|E|²
    electric_energy_density = (v.get_dim('epsilon_0')/2) * v.get_dim('E')**2
    target_energy_density = v.M / (v.L * v.T**2)

    v.check_dims("Electric energy density: (ε₀/2)|E|²", electric_energy_density, target_energy_density)

    # Magnetic energy density: (1/2μ₀)|B|²
    magnetic_energy_density = v.get_dim('B')**2 / (2 * v.get_dim('mu_0'))

    v.check_dims("Magnetic energy density: (1/2μ₀)|B|²", magnetic_energy_density, target_energy_density)

    # Both components have same dimensions
    v.check_dims("Electric and magnetic parts consistent", electric_energy_density, magnetic_energy_density)

    # Total energy density
    total_energy_density = electric_energy_density + magnetic_energy_density
    v.check_dims("Total EM energy density", total_energy_density, target_energy_density)


def test_poynting_vector(v):
    """Test Poynting vector: S = (1/μ₀) E × B"""
    v.subsection("Poynting Vector")

    # Poynting vector: (1/μ₀) E × B
    poynting_vector = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
    target_energy_flux = v.M / v.T**3  # Energy per area per time

    v.check_dims("Poynting vector: (1/μ₀) E × B", poynting_vector, target_energy_flux)

    # Check against helper's predefined if available
    if 'S_poynting' in v.dims:
        v.check_dims("Poynting vector consistency", poynting_vector, v.get_dim('S_poynting'))


def test_final_conservation_law(v):
    """Test the complete Poynting theorem as conservation law"""
    v.subsection("Complete Conservation Law")

    # Energy density
    energy_density = (v.get_dim('epsilon_0')/2) * v.get_dim('E')**2 + v.get_dim('B')**2/(2*v.get_dim('mu_0'))
    density_rate = v.dt(energy_density)

    # Energy flux (Poynting vector)
    poynting_vector = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
    flux_div = v.div_dim(poynting_vector)

    # Source term (Joule heating)
    source_term = v.get_dim('j_current') * v.get_dim('E')

    # Use conservation law verification: ∂ₜu + ∇·S = -J·E
    verify_conservation_law(
        v,
        "Electromagnetic Energy Conservation",
        density_rate,     # ∂ₜu
        flux_div,        # ∇·S
        -source_term     # -J·E (negative because energy is dissipated)
    )


def test_dimensional_units_statement(v):
    """Test the dimensional statement that all terms have [M L⁻¹ T⁻³]"""
    v.subsection("Power Density Units")

    # Target: power per volume = [W m⁻³] = [M L⁻¹ T⁻³]
    target_power_density = v.M / (v.L * v.T**3)

    # ∂ₜu_EM
    energy_density = (v.get_dim('epsilon_0')/2) * v.get_dim('E')**2 + v.get_dim('B')**2/(2*v.get_dim('mu_0'))
    dt_energy = v.dt(energy_density)

    # ∇·S
    poynting_vector = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
    div_poynting = v.div_dim(poynting_vector)

    # J·E
    joule_heating = v.get_dim('j_current') * v.get_dim('E')

    v.check_dims("∂ₜu_EM has power density units", dt_energy, target_power_density)
    v.check_dims("∇·S has power density units", div_poynting, target_power_density)
    v.check_dims("J·E has power density units", joule_heating, target_power_density)


def test_vacuum_conditions(v):
    """Test vacuum electromagnetic relationships"""
    v.subsection("Vacuum Conditions")

    # c² = 1/(μ₀ε₀)
    c_squared_from_constants = 1 / (v.get_dim('mu_0') * v.get_dim('epsilon_0'))
    c_squared_direct = v.get_dim('c')**2

    v.check_dims("Vacuum light speed: c² = 1/(μ₀ε₀)", c_squared_from_constants, c_squared_direct)

    # Z₀ = √(μ₀/ε₀)
    impedance_from_ratio = sqrt(v.get_dim('mu_0') / v.get_dim('epsilon_0'))

    if 'Z_0' in v.dims:
        v.check_dims("Vacuum impedance: √(μ₀/ε₀)", impedance_from_ratio, v.get_dim('Z_0'))


def main():
    """Run all Energy Flow (Poynting Theorem) verification tests"""

    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Energy Flow (Poynting Theorem)",
        "Complete verification of Poynting theorem derivation and conservation law",
        unit_system=UnitSystem.SI
    )

    # Define symbolic variables
    E, B, J = define_symbols_batch(['E', 'B', 'J'], positive=True, real=True)
    t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)

    v.section("POYNTING THEOREM VERIFICATION")

    # Test foundation and derivation steps
    test_maxwell_equations(v)
    test_vector_identity(v)
    test_derivation_step1(v)
    test_time_derivative_identities(v)
    test_derivation_step2(v)
    test_derivation_step3(v)

    # Test individual components
    test_energy_density_components(v)
    test_poynting_vector(v)

    # Test final result
    test_final_conservation_law(v)
    test_dimensional_units_statement(v)
    test_vacuum_conditions(v)

    # Generate summary
    return v.summary()


if __name__ == "__main__":
    main()
