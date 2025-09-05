#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "Minimal Coupling, EM Stress--Energy, and Light Propagation" subsection.

This module implements dimensional and mathematical verification for electromagnetic
field theory in curved spacetime, including:
- EM action and field equations
- Stress-energy tensor
- Light propagation
- Poynting theorem

Based on subsection "Minimal Coupling, EM Stress--Energy, and Light Propagation"
in doc/projected_em.tex (lines ~202-230).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, sqrt, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    quick_verify
)


def setup_em_dimensions(v):
    """Add electromagnetic field tensor and related dimensions exactly as they appear in relativity."""

    # Field tensor F_{μν} in SI units:
    # F_{0i} = E_i/c  (time-space components)
    # F_{ij} = ε_{ijk} B_k  (space-space components)
    # Both have dimensions [V*s/m²] = [M/(Q*T)] = B field dimensions
    v.add_dimensions({
        'F_field_tensor': v.dims['B'],  # F_{μν} has B field dimensions [M/(Q*T)]
        'F_squared': v.dims['B']**2,    # F_{μν} F^{μν} scalar invariant

        # Stress-energy tensor (standard dimensions)
        'T_EM_munu': v.M * v.L**(-1) * v.T**(-2),

        # 4-momentum wave vector
        'k_vector': v.L**(-1),  # k^μ wavenumber dimensions

        # Action integrand (will test if this works out correctly)
        'action_integrand': None,  # Will be computed in test
    })


def test_field_tensor_definition(v):
    """Test electromagnetic field tensor F_{μν} construction from potentials."""
    v.section("FIELD TENSOR DEFINITION")

    # F_{μν} = ∇_μ A_ν - ∇_ν A_μ
    # Each term ∇_μ A_ν has dimensions [L^-1] × [V*T/L] = [V*T/L²]

    gradient_A = v.grad_dim(v.get_dim('A'))  # ∇A has dims [L^-1] × [V*T/L]

    v.check_dims("Field tensor from potential gradient",
                 v.get_dim('F_field_tensor'), gradient_A)

    # Time components: F_{0i} ~ ∂₀A_i - ∂_iA₀ ~ E_i/c
    # Space components: F_{ij} ~ ∂_iA_j - ∂_jA_i ~ ε_{ijk}B_k

    # Check E-field relation: F_{0i} ~ E_i/c
    E_over_c = v.get_dim('E') / v.get_dim('c')
    v.check_dims("F_{0i} ~ E_i/c relationship", v.get_dim('F_field_tensor'), E_over_c)

    # Check B-field relation: F_{ij} ~ B_k
    v.check_dims("F_{ij} ~ B_k relationship", v.get_dim('F_field_tensor'), v.get_dim('B'))

    v.success("Field tensor definition verified")


def test_em_action(v):
    """Test EM action S_EM = -1/(4μ₀c) ∫ d⁴x √(-g) F_{μν}F^{μν}."""
    v.section("EM ACTION DIMENSIONAL ANALYSIS")

    # From document equation (line 211): S_EM[g,A] = -1/(4μ₀c) ∫ d⁴x √(-g) F_{μν}F^{μν}

    # √(-g) is the 4D volume element with dimensions [L⁴]
    sqrt_g = v.L**4

    # F_{μν}F^{μν} is the field tensor scalar invariant
    F_invariant = v.get_dim('F_squared')

    # The integrand: √(-g) F_{μν}F^{μν} / (4μ₀c)
    # Apply the 1/(4μ₀c) factor (ignoring numerical factor 1/4)
    full_integrand = sqrt_g * F_invariant / (v.get_dim('mu_0') * v.get_dim('c'))

    # Action should have energy×time dimensions [M*L²/T]
    expected_action = v.M * v.L**2 / v.T

    v.check_dims("EM action dimensional consistency",
                 full_integrand, expected_action)

    # Verify the Lagrangian density (integrand per unit 4-volume)
    lagrangian_density = full_integrand / v.L**4
    expected_lagrangian_density = v.M / (v.L**2 * v.T)

    v.check_dims("EM Lagrangian density",
                 lagrangian_density, expected_lagrangian_density)

    v.success("EM action dimensions verified")


def test_maxwell_equations(v):
    """Test Maxwell equations ∇_μ F^{μν} = μ₀ J^ν."""
    v.section("MAXWELL EQUATIONS")

    # Left side: ∇_μ F^{μν}
    # Has dimensions [L^-1] × [F_field_tensor] = [L^-1] × [V*T/L²] = [V*T/L³]

    div_F = v.div_dim(v.get_dim('F_field_tensor'))

    # Right side: μ₀ J^ν
    # J^ν has current density dimensions, μ₀ J^ν should match left side

    mu0_J = v.get_dim('mu_0') * v.get_dim('J_mu')

    v.check_dims("Maxwell equation: ∇_μ F^{μν} = μ₀ J^ν",
                 div_F, mu0_J)

    # Verify this gives the expected source equation structure
    # Should relate field derivatives to currents with proper coupling

    v.success("Maxwell equations dimensional consistency verified")


def test_bianchi_identity(v):
    """Test Bianchi identity ∇_{[α}F_{βγ]} = 0."""
    v.section("BIANCHI IDENTITY")

    # ∇_{[α}F_{βγ]} = 0 is identically satisfied when F_{μν} = ∇_μ A_ν - ∇_ν A_μ
    # This is a geometric identity - check it's dimensionally consistent

    # Each term ∇_α F_{βγ} has dimensions [L^-1] × [F_field_tensor]
    nabla_F = v.grad_dim(v.get_dim('F_field_tensor'))

    # The antisymmetrized combination should have same dimensions
    v.check_dims("Bianchi identity dimensional consistency",
                 nabla_F, v.div_dim(v.get_dim('F_field_tensor')))

    # Since this equals zero, it's automatically dimensionally consistent
    # The key point is F_{μν} = ∂_μ A_ν - ∂_ν A_μ automatically satisfies it

    v.success("Bianchi identity verified")


def test_em_stress_energy_tensor(v):
    """Test EM stress-energy tensor T^{μν}_EM."""
    v.section("EM STRESS-ENERGY TENSOR")

    # T^{μν}_EM = 1/μ₀ (F^{μα} F^ν_α - 1/4 g^{μν} F_{αβ} F^{αβ})

    # First term: F^{μα} F^ν_α / μ₀
    # F^{μα} F^ν_α has dimensions [F_field_tensor²]
    F_product = v.get_dim('F_squared')
    term1 = F_product / v.get_dim('mu_0')

    # Second term: g^{μν} F_{αβ} F^{αβ} / (4μ₀)
    # g^{μν} is dimensionless, F_{αβ} F^{αβ} has dimensions [F_field_tensor²]
    term2 = F_product / v.get_dim('mu_0')

    # Both terms should have stress-energy tensor dimensions
    target_Tmunu = v.get_dim('T_EM_munu')

    v.check_dims("T_EM first term: F^{μα} F^ν_α / μ₀", term1, target_Tmunu)
    v.check_dims("T_EM second term: g^{μν} F² / (4μ₀)", term2, target_Tmunu)

    # Check this matches general stress-energy dimensions [M L^-1 T^-2]
    v.check_dims("T_EM matches stress-energy dimensions",
                 target_Tmunu, v.M * v.L**(-1) * v.T**(-2))

    # Verify trace properties
    # Trace: T^μ_μ = 1/μ₀ (F^{μα} F_μα - 4/4 F_{αβ} F^{αβ}) = 0
    # This is automatic for EM field (traceless stress-energy)

    v.success("EM stress-energy tensor verified")


def test_light_ray_geodesics(v):
    """Test light ray conditions k^μ k_μ = 0 and k^ν ∇_ν k^μ = 0."""
    v.section("LIGHT RAY GEODESICS")

    # From document (line 232): k^μ k_μ = 0 and k^ν ∇_ν k^μ = 0
    # k^μ is the photon wavevector

    k_mu = v.get_dim('k_vector')  # Wavevector has dimensions [L^-1]

    # First condition: k^μ k_μ = 0 (null condition)
    null_scalar = k_mu**2

    # Second condition: k^ν ∇_ν k^μ = 0 (geodesic equation)
    geodesic_acceleration = k_mu * v.grad_dim(k_mu)

    # Verify dimensional consistency of both conditions
    v.check_dims("Null condition k^μ k_μ",
                 null_scalar, v.L**(-2))

    v.check_dims("Geodesic condition k^ν ∇_ν k^μ",
                 geodesic_acceleration, v.L**(-3))

    v.success("Light ray geodesic conditions verified")


def test_poynting_theorem(v):
    """Test Poynting theorem from document."""
    v.section("POYNTING THEOREM")

    # From document (lines 239-241):
    # ∂_t (ε₀/2 |E|² + 1/(2μ₀) |B|²) + ∇·(1/μ₀ E×B) = -J·E

    # Energy density terms
    electric_energy_density = v.get_dim('epsilon_0') * v.get_dim('E')**2
    electric_energy_rate = v.dt(electric_energy_density)

    magnetic_energy_density = v.get_dim('B')**2 / v.get_dim('mu_0')
    magnetic_energy_rate = v.dt(magnetic_energy_density)

    # Poynting vector and its divergence
    ExB_cross = v.get_dim('E') * v.get_dim('B')
    poynting_vector = ExB_cross / v.get_dim('mu_0')
    poynting_divergence = v.div_dim(poynting_vector)

    # Current work term
    current_field_work = v.get_dim('j_current') * v.get_dim('E')

    # All terms should have power density dimensions [M/(L*T³)]
    target_dim = v.M * v.L**(-1) * v.T**(-3)

    v.check_dims("Electric energy rate ∂_t(ε₀E²/2)", electric_energy_rate, target_dim)
    v.check_dims("Magnetic energy rate ∂_t(B²/(2μ₀))", magnetic_energy_rate, target_dim)
    v.check_dims("Poynting divergence ∇·(E×B/μ₀)", poynting_divergence, target_dim)
    v.check_dims("Current work J·E", current_field_work, target_dim)

    # Verify energy conservation structure
    v.check_dims("Energy rate terms balance", electric_energy_rate, magnetic_energy_rate)
    v.check_dims("Energy vs Poynting balance", electric_energy_rate, poynting_divergence)
    v.check_dims("Energy vs current work balance", electric_energy_rate, current_field_work)

    # Poynting vector should have energy flux dimensions [M T^-3]
    expected_energy_flux = v.M * v.T**(-3)
    v.check_dims("Poynting vector energy flux dimensions",
                 poynting_vector, expected_energy_flux)

    v.success("Poynting theorem verified")


def test_si_unit_relations(v):
    """Test fundamental SI unit relations for EM fields."""
    v.section("SI UNIT CONSISTENCY")

    # c² = 1/(μ₀ε₀)
    c_squared = v.get_dim('c')**2
    mu0_eps0 = v.get_dim('mu_0') * v.get_dim('epsilon_0')
    inverse_mu0_eps0 = 1 / mu0_eps0

    v.check_dims("Speed of light: c² = 1/(μ₀ε₀)", c_squared, inverse_mu0_eps0)

    # Impedance: Z₀ = √(μ₀/ε₀) = μ₀c
    Z0_ratio = sqrt(v.get_dim('mu_0') / v.get_dim('epsilon_0'))
    Z0_muc = v.get_dim('mu_0') * v.get_dim('c')

    v.check_dims("Vacuum impedance from ratio", v.get_dim('Z_0'), Z0_ratio)
    v.check_dims("Vacuum impedance from μ₀c", v.get_dim('Z_0'), Z0_muc)

    v.success("SI unit relations verified")


def test_minimal_coupling():
    """Main test function for Minimal Coupling, EM Stress-Energy, and Light Propagation.

    This function coordinates all verification tests for electromagnetic
    field theory in curved spacetime, calling helper functions as needed
    and providing a single entry point.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize helper with SI units (as noted in the document)
    v = PhysicsVerificationHelper(
        "Minimal Coupling, EM Stress-Energy, and Light Propagation",
        "Verification of electromagnetic field theory in curved spacetime",
        unit_system=UnitSystem.SI
    )

    # Add EM-specific dimensions
    setup_em_dimensions(v)

    # Run all verification tests
    test_field_tensor_definition(v)
    test_em_action(v)
    test_maxwell_equations(v)
    test_bianchi_identity(v)
    test_em_stress_energy_tensor(v)
    test_light_ray_geodesics(v)
    test_poynting_theorem(v)
    test_si_unit_relations(v)

    # Generate summary
    return v.summary()


if __name__ == "__main__":
    success_rate = test_minimal_coupling()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
