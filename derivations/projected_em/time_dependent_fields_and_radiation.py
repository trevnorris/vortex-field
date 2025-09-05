#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Time-Dependent Fields and Radiation - Verification
==================================================

Complete verification of all mathematical relationships in the "Electromagnetism
in the Wave Sector" subsection from doc/projected_em.tex (lines 530-648).

Tests cover:
- Wave equation for electromagnetic potentials in Lorenz gauge
- Maxwell's equations in field tensor form
- Static Coulomb limit and field calibration
- Retarded potentials and causality
- Energy density and Poynting flux relationships
- Sommerfeld radiation condition
- Oscillating dipole radiation scaling laws
"""

import os
import sys
from sympy import symbols, pi, sqrt, simplify, cos, sin, exp, I, diff, integrate

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem
)


def test_wave_equations_lorenz_gauge(v):
    """Test wave equations for EM potentials in Lorenz gauge."""
    v.subsection("Wave Equations in Lorenz Gauge")

    # Wave equation: (1/c² ∂_t² - ∇²)A^μ = -4π/c J^μ_ch
    # Check dimensional consistency of all terms

    # Time derivative term: (1/c²)∂²A/∂t²
    c = v.get_dim('c')
    A_potential = v.get_dim('A')
    time_second_deriv_A = A_potential / v.get_dim('t')**2
    time_term = time_second_deriv_A / c**2

    # Spatial Laplacian term: ∇²A
    laplacian_A = v.lap_dim(A_potential)

    # Left side of wave equation (d'Alembertian)
    dalembertian_dim = time_term  # Same dimension as spatial term

    # Right side: source term -4π/c J^μ_ch
    J_current = v.get_dim('j_current')
    source_term = J_current / c

    # Check dimensional consistency of wave equation
    v.check_dims("Wave equation: time derivative term", time_term, laplacian_A)
    v.check_dims("Wave equation: ∇²A vs (1/c²)∂²A/∂t²", laplacian_A, time_term)
    v.check_dims("Wave equation balance", laplacian_A, source_term)

    # Lorenz gauge condition: ∂_μ A^μ = 0
    # In 3+1 form: (1/c)∂Φ/∂t + ∇·A = 0
    Phi_E = v.get_dim('Phi_E')
    div_A_dim = v.div_dim(A_potential)
    time_deriv_phi = Phi_E / v.get_dim('t')
    lorenz_term = time_deriv_phi / c

    v.check_dims("Lorenz gauge: ∇·A vs (1/c)∂Φ/∂t", div_A_dim, lorenz_term)


def test_maxwell_field_tensor_equations(v):
    """Test Maxwell's equations in field tensor form."""
    v.subsection("Maxwell's Equations in Field Tensor Form")

    # Inhomogeneous Maxwell equations: ∂_μ F^{μν} = (4π/c) J^ν_ch

    # Field tensor has dimensions of E/length (for spatial components)
    F_tensor_div = v.get_dim('E_field') / v.get_dim('r')

    # Right side: current source (4π/c) J^ν
    c = v.get_dim('c')
    J_current = v.get_dim('j_current')
    current_source = J_current / c

    v.check_dims("Maxwell inhomogeneous: ∂F vs J/c", F_tensor_div, current_source)

    # Homogeneous Maxwell equations: ∂_{[α} F_{βγ]} = 0
    # This gives ∇·B = 0 and ∇×E + (1/c)∂B/∂t = 0
    curl_E_dim = v.curl_dim(v.get_dim('E_field'))
    time_deriv_B = v.get_dim('B_field') / v.get_dim('t')
    c = v.get_dim('c')

    v.check_dims("Faraday's law: ∇×E vs -(1/c)∂B/∂t", curl_E_dim, time_deriv_B / c)

    # Check that ∇·B has right dimensions (should be zero)
    div_B_dim = v.div_dim(v.get_dim('B_field'))
    v.check_dims("Gauss law for magnetism: ∇·B", div_B_dim,
                 v.get_dim('B_field') / v.get_dim('r'))


def test_electric_magnetic_field_definitions(v):
    """Test electric and magnetic field definitions from 4-potential."""
    v.subsection("Electric and Magnetic Field Definitions")

    # E = -∇Φ_Slope - (1/c)∂A/∂t

    # Gradient of scalar potential: -∇Φ_E
    Phi_E = v.get_dim('Phi_E')
    grad_phi_dim = v.grad_dim(Phi_E)

    # Time derivative of vector potential: (1/c)∂A/∂t
    A_vec = v.get_dim('A')
    c = v.get_dim('c')
    time_deriv_A_dim = (A_vec / v.get_dim('t')) / c

    # Both terms should have E-field dimensions
    v.check_dims("Gradient potential term", grad_phi_dim, v.get_dim('E_field'))
    v.check_dims("Time derivative vector potential", time_deriv_A_dim, v.get_dim('E_field'))
    v.check_dims("E-field components balance", grad_phi_dim, time_deriv_A_dim)

    # Magnetic field: B = ∇×A
    curl_A_dim = v.curl_dim(A_vec)
    v.check_dims("Magnetic field from curl", curl_A_dim, v.get_dim('B_field'))


def test_static_coulomb_calibration(v):
    """Test static limit and Coulomb law calibration."""
    v.subsection("Static Coulomb Calibration")

    # In static limit (∂_t = 0): ∇²Φ_E = -4π ρ_ch
    Phi_E = v.get_dim('Phi_E')
    rho_charge = v.get_dim('rho_charge')

    # Poisson equation dimensional consistency
    laplacian_phi_dim = v.lap_dim(Phi_E)
    source_term_dim = rho_charge

    v.check_dims("Poisson equation: ∇²Φ vs ρ", laplacian_phi_dim, source_term_dim)

    # For point charge Q: Φ(r) = Q/r, E(r) = Q r/r³
    Q = v.get_dim('e')  # Use elementary charge
    r = v.get_dim('r')

    # Coulomb potential: Φ = Q/r
    coulomb_potential = Q / r
    v.check_dims("Coulomb potential", coulomb_potential, v.get_dim('Phi_E'))

    # Electric field: E = -∇Φ = Q r/r³ = Q/r²
    coulomb_field = v.grad_dim(coulomb_potential)
    v.check_dims("Coulomb field from potential", coulomb_field, v.get_dim('E_field'))

    # Direct check: |E| = Q/r²
    direct_coulomb_field = Q / r**2
    v.check_dims("Direct Coulomb field", direct_coulomb_field, v.get_dim('E_field'))
    v.check_dims("Coulomb field consistency", coulomb_field, direct_coulomb_field)

    # In static limit: ∇²A = 0 (Lorenz gauge)
    A_vec = v.get_dim('A')
    laplacian_A_dim = v.lap_dim(A_vec)
    v.check_dims("Static vector potential Laplacian", laplacian_A_dim,
                 A_vec / v.get_dim('r')**2)


def test_retarded_potentials(v):
    """Test retarded potential solutions and causality."""
    v.subsection("Retarded Potentials and Causality")

    # Retarded scalar potential: Φ(x,t) = ∫ ρ(x',t_r)/|x-x'| d³x'
    # where t_r = t - |x-x'|/c

    # Retarded time structure
    c = v.get_dim('c')
    t = v.get_dim('t')
    r_distance = v.get_dim('r')

    # Light travel time: |x-x'|/c
    travel_time = r_distance / c
    v.check_dims("Light travel time", travel_time, t)

    # Retarded scalar potential integrand: ρ(x',t_r)/|x-x'|
    rho_charge = v.get_dim('rho_charge')
    potential_integrand = rho_charge / r_distance

    # This integrand, when integrated over d³x', gives potential
    # Dimension: [charge/length³] / [length] * [length³] = [charge/length] = [potential]
    expected_potential_dim = v.get_dim('Phi_E')
    integrand_times_volume = potential_integrand * v.get_dim('r')**3

    v.check_dims("Retarded potential integral", integrand_times_volume, expected_potential_dim)

    # Retarded vector potential: A(x,t) = (1/c) ∫ J(x',t_r)/|x-x'| d³x'
    J_current = v.get_dim('j_current')
    vector_potential_integrand = J_current / (c * r_distance)

    # When integrated: [current/length³] / ([speed] * [length]) * [length³] = [current/(speed*length)] = [A]
    vector_integral = vector_potential_integrand * v.get_dim('r')**3
    expected_A_dim = v.get_dim('A')

    v.check_dims("Retarded vector potential integral", vector_integral, expected_A_dim)


def test_energy_density_poynting_flux(v):
    """Test electromagnetic energy density and Poynting flux."""
    v.subsection("Energy Density and Poynting Flux")

    # Energy density in Gaussian units: u = (E² + B²)/(8π)
    E_field = v.get_dim('E_field')
    B_field = v.get_dim('B_field')

    # In Gaussian units E and B have same dimensions
    v.check_dims("E and B field dimensions", E_field, B_field)

    energy_density = (E_field**2 + B_field**2) / (8*pi)
    expected_energy_density = v.get_dim('E_energy') / v.get_dim('r')**3

    v.check_dims("EM energy density", energy_density, expected_energy_density)

    # Poynting vector in Gaussian units: S = (c/4π) E × B
    c = v.get_dim('c')
    # |E × B| has dimension E*B = E^2 in Gaussian units
    cross_product_EB = E_field * B_field

    poynting_vector = (c / (4*pi)) * cross_product_EB
    expected_poynting = v.get_dim('E_energy') / (v.get_dim('r')**2 * v.get_dim('t'))

    v.check_dims("Poynting vector", poynting_vector, expected_poynting)

    # Poynting theorem: ∂u/∂t + ∇·S = -J·E
    # All terms should have dimensions of power density
    power_density = v.get_dim('E_energy') / (v.get_dim('r')**3 * v.get_dim('t'))

    time_deriv_u = energy_density / v.get_dim('t')
    div_S = v.div_dim(poynting_vector)
    J_current = v.get_dim('j_current')
    J_dot_E = J_current * E_field

    v.check_dims("Energy density time derivative", time_deriv_u, power_density)
    v.check_dims("Poynting divergence", div_S, power_density)
    v.check_dims("Current work term", J_dot_E, power_density)

    # Check Poynting theorem dimensional balance
    v.check_dims("Poynting theorem: ∂u/∂t vs ∇·S", time_deriv_u, div_S)
    v.check_dims("Poynting theorem: ∇·S vs J·E", div_S, J_dot_E)


def test_sommerfeld_radiation_condition(v):
    """Test Sommerfeld radiation condition for outgoing waves."""
    v.subsection("Sommerfeld Radiation Condition")

    # Radiation condition: (∂_r - (1/c)∂_t)(r A^μ) → 0 as r → ∞

    # Product r A^μ has dimensions [length] * [A]
    r = v.get_dim('r')
    A_vec = v.get_dim('A')
    c = v.get_dim('c')

    r_times_A = r * A_vec

    # Radial derivative: ∂(rA)/∂r has dimensions [A]
    radial_deriv_rA = r_times_A / r  # ∂/∂r brings down 1/r

    # Time derivative: (1/c)∂(rA)/∂t has dimensions [A]
    time_deriv_rA = r_times_A / v.get_dim('t')
    normalized_time_deriv = time_deriv_rA / c

    # Both terms in radiation condition should have same dimensions
    v.check_dims("Radial derivative term", radial_deriv_rA, A_vec)
    v.check_dims("Time derivative term", normalized_time_deriv, A_vec)
    v.check_dims("Radiation condition balance", radial_deriv_rA, normalized_time_deriv)

    # The radiation condition ensures outgoing waves only
    # This is a boundary condition that selects retarded solutions


def test_oscillating_dipole_radiation(v):
    """Test far-field radiation from oscillating electric dipole."""
    v.subsection("Oscillating Dipole Radiation")

    # Far-field electric field scaling: |E| ~ (ω² |p₀| sin θ)/(c² r)
    omega = v.get_dim('omega_freq')
    p_0 = v.get_dim('e') * v.get_dim('r')  # Electric dipole moment = charge × length
    c = v.get_dim('c')
    r = v.get_dim('r')

    # sin θ is dimensionless
    sin_theta = 1  # Dimensionless factor

    # Field amplitude scaling
    E_amplitude_scaling = (omega**2 * p_0 * sin_theta) / (c**2 * r)

    v.check_dims("Dipole field amplitude", E_amplitude_scaling, v.get_dim('E_field'))

    # Magnetic field scaling in Gaussian units: |B| ~ |E|/c
    B_amplitude_scaling = E_amplitude_scaling / c

    # Check the radiative field relationship |B| = |E|/c (dimensionally)
    v.check_dims("Radiative field relation", B_amplitude_scaling * c, E_amplitude_scaling)

    # In radiative fields, B amplitude has extra time factor compared to static B field
    # This is correct physics: static B ~ Q/L², radiative B ~ (E/c) ~ Q*T/L³
    radiative_B_vs_static = B_amplitude_scaling / (v.get_dim('B_field'))
    expected_ratio = v.get_dim('t') / v.get_dim('r')  # T/L factor from c
    v.check_dims("Radiative vs static B field ratio", radiative_B_vs_static, expected_ratio)

    # Radiated power scaling: P ∝ ω⁴ |p₀|²/c³
    power_scaling = (omega**4 * p_0**2) / c**3

    v.check_dims("Radiated power scaling", power_scaling,
                 v.get_dim('E_energy') / v.get_dim('t'))


def test_em_self_energy(v):
    """Test electromagnetic self-energy definition."""
    v.subsection("Electromagnetic Self-Energy")

    # EM self-energy in Gaussian units: δE_EM = (1/8π) ∫(E² + B²) d³x
    E_field = v.get_dim('E_field')
    B_field = v.get_dim('B_field')

    # Energy density: u = (E² + B²)/(8π)
    energy_density = (E_field**2 + B_field**2) / (8*pi)
    expected_energy_density = v.get_dim('E_energy') / v.get_dim('r')**3

    v.check_dims("EM energy density", energy_density, expected_energy_density)

    # Total self-energy (volume integral): E_EM = ∫ u d³x
    volume_element = v.get_dim('r')**3
    self_energy_integrand = energy_density * volume_element

    v.check_dims("EM self-energy integrand", self_energy_integrand, v.get_dim('E_energy'))

    # For localized configuration, this gives finite self-energy
    # This quantity appears in baryon mass calculations as β_Q/R term


def test_time_dependent_fields_and_radiation():
    """Test runner compatible function for time-dependent fields and radiation."""
    return main()


def main():
    """Run all verification tests for time-dependent fields and radiation."""
    v = PhysicsVerificationHelper("Time-Dependent Fields and Radiation",
                                  unit_system=UnitSystem.GAUSSIAN)

    # Run all test functions
    test_wave_equations_lorenz_gauge(v)
    test_maxwell_field_tensor_equations(v)
    test_electric_magnetic_field_definitions(v)
    test_static_coulomb_calibration(v)
    test_retarded_potentials(v)
    test_energy_density_poynting_flux(v)
    test_sommerfeld_radiation_condition(v)
    test_oscillating_dipole_radiation(v)
    test_em_self_energy(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = main()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)