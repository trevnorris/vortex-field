#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "Conservation Laws and Aether Drainage" subsection.

This module implements dimensional and mathematical verification for the subsection
covering mass continuity, momentum balance, and energy conservation with drainage
terms, as outlined in TEST.md.

Based on the mathematical framework document, Section covering conservation laws.
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    verify_conservation_law,
    quick_verify
)


def test_mass_continuity_local_form(v):
    """Test local mass continuity with drainage: dt rho + div(rho v) = -sum M_dot_i delta3."""
    # partial_t rho + nabla·(rho v) = -sum_i M_dot_i delta^(3)(r-r_i)

    lhs_rate = v.dt(v.get_dim('rho'))                           # [M L^-3 T^-1]
    lhs_flux = v.div_dim(v.get_dim('rho') * v.get_dim('v'))     # [M L^-3 T^-1]
    rhs_sink = v.get_dim('M_dot_i') * v.get_dim('delta3')       # [M L^-3 T^-1]

    target_mass_rate = v.M / (v.L**3 * v.T)

    v.check_dims("Mass continuity: partial_t rho", lhs_rate, target_mass_rate)
    v.check_dims("Mass continuity: div(rho v)", lhs_flux, target_mass_rate)
    v.check_dims("Mass continuity: sink term", rhs_sink, target_mass_rate)

    # Verify conservation law structure
    verify_conservation_law(v, "Mass continuity with drainage", lhs_rate, lhs_flux, rhs_sink)

    v.success("Local mass continuity verified")


def test_mass_continuity_linearized_form(v):
    """Test linearized mass continuity: dt delta_rho + rho_0 div delta_v = -sum M_dot_i delta3."""
    # partial_t delta_rho + rho_0 nabla·delta_v = -sum_i M_dot_i delta^(3)

    lhs_rate_lin = v.dt(v.get_dim('delta_rho'))                           # [M L^-3 T^-1]
    lhs_flux_lin = v.get_dim('rho_0') * v.div_dim(v.get_dim('v'))         # rho_0 div delta_v
    rhs_sink_lin = v.get_dim('M_dot_i') * v.get_dim('delta3')             # [M L^-3 T^-1]

    target_mass_rate = v.M / (v.L**3 * v.T)

    v.check_dims("Linearized mass: partial_t delta_rho", lhs_rate_lin, target_mass_rate)
    v.check_dims("Linearized mass: rho_0 div delta_v", lhs_flux_lin, target_mass_rate)
    v.check_dims("Linearized mass: sink term", rhs_sink_lin, target_mass_rate)

    v.success("Linearized mass continuity verified")


def test_global_mass_balance(v):
    """Test global mass balance dimensional consistency."""
    # d/dt ∫_V rho dV + ∮_∂V rho v·dA = -sum_{i∈V} M_dot_i

    # Volume integral: d/dt ∫ ρ dV → [M T^-1]
    mass_rate_V = v.get_dim('rho') * v.L**3 / v.T

    # Surface integral: ∮ ρ v·dA → [M T^-1]
    # rho [M L^-3] * v [L T^-1] * dA [L^2] = [M T^-1]
    mass_flux_S = v.get_dim('rho') * v.get_dim('v') * v.get_dim('dA')

    # Sinks: sum M_dot_i → [M T^-1]
    total_sinks = v.get_dim('M_dot_i')

    target_mass_flow = v.M / v.T

    v.check_dims("Global mass: d/dt ∫ρ dV", mass_rate_V, target_mass_flow)
    v.check_dims("Global mass: surface flux", mass_flux_S, target_mass_flow)
    v.check_dims("Global mass: total sinks", total_sinks, target_mass_flow)

    v.success("Global mass balance verified")


def test_momentum_balance_all_terms(v):
    """Test all terms in momentum balance equation."""
    # partial_t(rho v) + nabla·(rho v⊗v + p I) = -rho nabla Phi_g - rho nabla Q - sum M_dot_i v_*i delta3 + f_ext

    # LHS terms
    mom_rate = v.dt(v.get_dim('rho') * v.get_dim('v'))                    # [M L^-2 T^-2]
    mom_flux = v.div_dim(v.get_dim('rho') * v.get_dim('v') * v.get_dim('v'))  # div(rho v⊗v)
    press_grad = v.grad_dim(v.get_dim('p'))                               # nabla p

    # RHS force terms
    grav_force = v.get_dim('rho') * v.grad_dim(v.get_dim('Phi_g'))        # rho nabla Phi_g
    quant_force = v.get_dim('rho') * v.grad_dim(v.get_dim('Q'))           # rho nabla Q
    drain_mom = v.get_dim('M_dot_i') * v.get_dim('v_sink') * v.get_dim('delta3')  # drainage momentum
    f_ext = v.get_dim('f_ext')                                            # external force density

    target_mom = v.M / (v.L**2 * v.T**2)

    # Test all momentum terms
    momentum_terms = [
        ("partial_t(rho v)", mom_rate),
        ("div(rho v tensor v)", mom_flux),
        ("nabla p", press_grad),
        ("rho nabla Phi_g", grav_force),
        ("rho nabla Q", quant_force),
        ("drain momentum", drain_mom),
        ("f_ext", f_ext)
    ]

    for name, term in momentum_terms:
        v.check_dims(f"Momentum term: {name}", term, target_mom)

    v.success("Momentum balance all terms verified")


def test_global_momentum_balance(v):
    """Test global momentum balance dimensional consistency."""
    # d/dt ∫_V ρv dV + ∮_∂V (ρv⊗v + pI)·dA = -∫_V ρ∇Φg dV - ∫_V ρ∇Q dV - sum_{i∈V} M_dot_i v_*i + ∫_V f_ext dV

    # Volume momentum rate: d/dt ∫ ρv dV → [M L T^-2]
    mom_rate_V = v.get_dim('rho') * v.get_dim('v') * v.L**3 / v.T

    # Surface momentum flux: ∮ (ρv⊗v + pI)·dA → [M L T^-2]
    # This is trickier - we need the tensor contraction with area normal
    mom_flux_surface = (v.get_dim('rho') * v.get_dim('v') * v.get_dim('v') + v.get_dim('p')) * v.get_dim('dA')

    # Force terms integrated over volume → [M L T^-2]
    target_force = v.M * v.L / v.T**2

    v.check_dims("Global momentum: d/dt integral rho v dV", mom_rate_V, target_force)
    v.check_dims("Global momentum: surface flux", mom_flux_surface, target_force)

    v.success("Global momentum balance verified")


def test_energy_density_components(v):
    """Test all components of energy density."""
    # e = (1/2)ρv² + u(ρ) + ρΦg + e_Q

    # Energy density components
    kin_energy = v.get_dim('rho') * v.get_dim('v')**2                     # (1/2)ρv²
    int_energy = v.get_dim('u')                                           # internal energy density u(ρ)
    grav_energy = v.get_dim('rho') * v.get_dim('Phi_g')                   # ρΦg
    quant_energy = v.get_dim('e_Q')                                       # quantum energy density

    target_edens = v.M / (v.L * v.T**2)

    energy_components = [
        ("kinetic", kin_energy),
        ("internal", int_energy),
        ("gravitational", grav_energy),
        ("quantum", quant_energy)
    ]

    for name, term in energy_components:
        v.check_dims(f"Energy density: {name}", term, target_edens)

    v.success("Energy density components verified")


def test_energy_flux_terms(v):
    """Test energy flux (Poynting-like) terms."""
    # S = (e + p)v + S_Q

    # Convective flux: (e + p)v
    S_conv = (v.get_dim('e') + v.get_dim('p')) * v.get_dim('v')

    # Quantum/dispersion flux (if present)
    S_Q = v.get_dim('S_flux')

    target_eflux = v.M / v.T**3

    v.check_dims("Energy flux: convective (e+p)v", S_conv, target_eflux)
    v.check_dims("Energy flux: quantum S_Q", S_Q, target_eflux)

    v.success("Energy flux terms verified")


def test_local_energy_balance(v):
    """Test local energy balance equation."""
    # partial_t e + nabla·S = -ρv·∇Φg - sum M_dot_i ε_i delta3 + Π_Q

    # LHS terms
    e_rate = v.dt(v.get_dim('e'))                                         # [M L^-1 T^-3]
    div_S = v.div_dim(v.get_dim('S_flux'))                                # [M L^-1 T^-3]

    # RHS terms
    grav_work = v.get_dim('rho') * v.get_dim('v') * v.grad_dim(v.get_dim('Phi_g'))  # ρv·∇Φg
    drain_power = v.get_dim('M_dot_i') * v.get_dim('eps_spec') * v.get_dim('delta3')  # drainage power
    Pi_Q = v.get_dim('Pi_Q')                                              # quantum work term

    target_power_dens = v.M / (v.L * v.T**3)

    energy_balance_terms = [
        ("partial_t e", e_rate),
        ("nabla dot S", div_S),
        ("rho v dot nabla Phi_g", grav_work),
        ("drain power", drain_power),
        ("Pi_Q", Pi_Q)
    ]

    for name, term in energy_balance_terms:
        v.check_dims(f"Energy balance: {name}", term, target_power_dens)

    v.success("Local energy balance verified")


def test_global_energy_balance(v):
    """Test global energy balance dimensional consistency."""
    # d/dt ∫_V e dV + ∮_∂V S·dA = -∫_V ρv·∇Φg dV - sum_{i∈V} M_dot_i ε_i + ∫_V Π_Q dV

    # Volume energy rate: d/dt ∫ e dV → [M L^2 T^-3] (power)
    E_rate_V = v.get_dim('e') * v.L**3 / v.T

    # Surface energy flux: ∮ S·dA → [M L^2 T^-3]
    Surf_energy_flux = v.get_dim('S_flux') * v.get_dim('dA')

    # Integrated power terms → [M L^2 T^-3]
    target_power = v.M * v.L**2 / v.T**3

    v.check_dims("Global energy: d/dt integral e dV", E_rate_V, target_power)
    v.check_dims("Global energy: surface S dot dA", Surf_energy_flux, target_power)

    v.success("Global energy balance verified")


def test_4d_3d_consistency(v):
    """Test 4D→3D consistency for drainage terms."""
    # ∫ M_dot_i δ^(4)(r_4-r_4,i) dw = M_dot_i δ^(3)(r-r_i)

    # 4D form integrated over w: [M T^-1][L^-4][L] = [M L^-3 T^-1]
    delta4_integrated = v.get_dim('M_dot_i') * v.get_dim('delta4') * v.get_dim('w')

    # 3D form: [M T^-1][L^-3] = [M L^-3 T^-1]
    delta3_form = v.get_dim('M_dot_i') * v.get_dim('delta3')

    v.check_dims("4D to 3D consistency: integrated delta4", delta4_integrated, delta3_form)

    v.success("4D to 3D consistency verified")


def test_noether_current_naming(v):
    """Test Noether-style current dimensional consistency (optional)."""
    # Mass current: J_m = ρv
    # Momentum flux tensor: Π_ij = ρv_i v_j + pδ_ij
    # Energy current: S

    # Mass current J_m = ρv → [M L^-2 T^-1]
    J_m = v.get_dim('rho') * v.get_dim('v')
    v.check_dims("Mass current J_m", J_m, v.M / (v.L**2 * v.T))

    # Momentum flux tensor → [M L^-1 T^-2]
    Pi_momentum = v.get_dim('rho') * v.get_dim('v') * v.get_dim('v') + v.get_dim('p')
    v.check_dims("Momentum flux tensor", Pi_momentum, v.M / (v.L * v.T**2))

    # Energy current S → [M T^-3]
    S_energy = v.get_dim('S_flux')
    v.check_dims("Energy current S", S_energy, v.M / v.T**3)

    v.success("Noether current naming verified")


def test_sanity_reductions_and_diagnostics(v):
    """Test sanity reductions and diagnostic checks."""

    # 1) No drains -> standard conservation (diagnostic: should detect mismatch when comparing sink to 0)
    rhs_sink = v.get_dim('M_dot_i') * v.get_dim('delta3')
    v.info("DIAGNOSTIC: Testing no-drain condition (should detect non-zero sink)...")
    ok_nodrain = v.check_dims("No-drain diagnostic", rhs_sink, 0, record=False, verbose=False)
    quick_verify("Diagnostic: mass sink vanishes when M_dot=0", not ok_nodrain, helper=v, expected_failure=True)

    # 2) Energy drain must include specific energy factor
    e_rate = v.dt(v.get_dim('e'))
    bad_drain = v.get_dim('M_dot_i') * v.get_dim('delta3')  # missing eps_spec
    v.info("DIAGNOSTIC: Testing energy drain without specific energy (should be dimensionally wrong)...")
    ok_bad_drain = v.check_dims("Energy drain missing eps — diagnostic",
                                e_rate, bad_drain, record=False, verbose=False)
    quick_verify("Caught: energy drain must include specific energy factor", not ok_bad_drain, helper=v, expected_failure=True)

    v.success("Sanity reductions and diagnostics verified")


def test_drainage_rate_formula(v):
    """Test specific drainage rate formula from framework: Ṁ_i ≈ ρ_{4D}^0 Γ ξ_c²."""
    # Framework line 769-772: M_dot_i ~ rho_4D_0 * Gamma * xi_c^2

    # Individual components
    rho_4D_0 = v.get_dim('rho_4')             # Background 4D density [M L^-4]
    Gamma = v.get_dim('Gamma')                # Circulation quantum [L^2 T^-1] (from helper.py)
    xi_c_sq = v.get_dim('xi_c')**2            # Healing length squared [L^2]

    # Formula from framework
    M_dot_formula = rho_4D_0 * Gamma * xi_c_sq

    # Generic drainage rate
    M_dot_generic = v.get_dim('M_dot_i')      # [M T^-1]

    v.check_dims("Drainage rate formula", M_dot_formula, M_dot_generic)

    # Verify dimensional breakdown matches framework
    target_M_per_T = v.M / v.T
    v.check_dims("Drainage formula breakdown", M_dot_formula, target_M_per_T)

    v.success("Drainage rate formula verified")


def test_energy_barrier_formula(v):
    """Test energy barrier formula: ΔE ≈ (ρ_{4D}^0 Γ² ξ_c²)/(4π) ln(L/ξ_c)."""
    # Framework line 774-776: reconnection energy barrier

    # Components
    rho_4D_0 = v.get_dim('rho_4')             # [M L^-4]
    Gamma_sq = v.get_dim('Gamma')**2          # [L^4 T^-2] (from helper.py)
    xi_c_sq = v.get_dim('xi_c')**2            # [L^2]
    L_cutoff = v.L                            # Outer cutoff length [L] (dimension)
    pi = 1                                    # Dimensionless (sympy.pi is dimensionless)

    # Energy barrier formula: ΔE ~ (rho_4D_0 * Gamma^2 * xi_c^2) / (4π) * ln(L/xi_c)
    # Note: ln(L/xi_c) is dimensionless
    Delta_E = (rho_4D_0 * Gamma_sq * xi_c_sq) / (4 * pi)  # Ignoring ln factor for dimensions

    target_energy = v.M * v.L**2 / v.T**2     # [M L^2 T^-2]

    v.check_dims("Energy barrier formula", Delta_E, target_energy)

    v.success("Energy barrier formula verified")


def test_quantum_power_identity(v):
    """Test quantum power identity: -ρv·∇Q = -∇·S_Q + Π_Q."""
    # Framework line 832-833: quantum work-rate identity

    # Left hand side: -ρv·∇Q
    lhs = v.get_dim('rho') * v.get_dim('v') * v.grad_dim(v.get_dim('Q'))  # [M L^-1 T^-3]

    # Right hand side terms: -∇·S_Q + Π_Q
    div_S_Q = v.div_dim(v.get_dim('S_flux'))     # ∇·S_Q [M L^-1 T^-3]
    Pi_Q = v.get_dim('Pi_Q')                     # Quantum work power [M L^-1 T^-3]

    target_power_dens = v.M / (v.L * v.T**3)     # Power density [M L^-1 T^-3]

    # Check all terms have same dimensions
    v.check_dims("Quantum identity LHS: -ρv·∇Q", lhs, target_power_dens)
    v.check_dims("Quantum identity RHS: ∇·S_Q", div_S_Q, target_power_dens)
    v.check_dims("Quantum identity RHS: Π_Q", Pi_Q, target_power_dens)

    # Verify identity structure (all terms match)
    v.check_dims("Quantum power identity consistency", lhs, div_S_Q)

    v.success("Quantum power identity verified")


def test_enhanced_4d_3d_matching(v):
    """Test enhanced 4D-3D matching conditions with slab integration."""
    # Framework line 857-859: ∫_{-ξ_c}^{+ξ_c} ∂_w(ρ_bulk v_w) dw = -∑ M_dot_i δ³
    # This integrates flux jump across slab thickness using dimensionless window

    # Left side: ∫ ∂_w(ρ_bulk v_w) dw over slab thickness
    rho_bulk = v.get_dim('rho_4')                    # Bulk density [M L^-4]
    v_w = v.get_dim('v')                             # Bulk velocity in w direction [L T^-1]

    # ∂_w(ρ_bulk v_w): derivative of mass flux
    # [M L^-4] * [L T^-1] / [L] = [M L^-4 T^-1]
    bulk_flux_derivative = (rho_bulk * v_w) / v.get_dim('w')  # [M L^-4 T^-1]

    # Integration over slab thickness: ∫_{-ξ_c}^{+ξ_c} dw gives [L]
    # So: [M L^-4 T^-1] * [L] = [M L^-3 T^-1]
    xi_c_thickness = v.get_dim('xi_c')               # Slab thickness [L]
    lhs_integrated = bulk_flux_derivative * xi_c_thickness  # [M L^-3 T^-1]

    # Right side: -∑ M_dot_i δ³
    rhs_sinks = v.get_dim('M_dot_i') * v.get_dim('delta3')  # [M L^-3 T^-1]

    v.check_dims("4D-3D matching: slab integration", lhs_integrated, rhs_sinks)

    v.success("Enhanced 4D-3D matching verified")


def test_conservation_laws_and_aether_drainage():
    """
    Main test function implementing all verification categories from TEST.md.

    Tests 5 main categories:
    A) Mass continuity (4 tests)
    B) Momentum balance (7 tests)
    C) Energy conservation (8 tests)
    D) 4D<->3D consistency (1 test)
    E) Sanity/diagnostics (3 tests)
    """
    v = PhysicsVerificationHelper(
        "Conservation Laws & Aether Drainage",
        "Mass, momentum, energy balances with sinks; local and global forms",
        unit_system=UnitSystem.SI
    )

    v.section_header("Testing Conservation Laws and Aether Drainage")

    # Register additional symbols needed for conservation laws
    v.add_dimensions({
        # Pressure, internal energy density, enthalpy density
        'p': v.M/(v.L*v.T**2),           # [M L^-1 T^-2]
        'u': v.M/(v.L*v.T**2),           # internal energy density [M L^-1 T^-2]
        'w': v.M/(v.L*v.T**2),           # enthalpy density (if used)
        'e': v.M/(v.L*v.T**2),           # total energy density [M L^-1 T^-2]
        'e_Q': v.M/(v.L*v.T**2),         # quantum energy density [M L^-1 T^-2]
        'S_flux': v.M/(v.T**3),          # energy flux [M T^-3]

        # Surface/volume elements for global checks
        'dA': v.L**2,                    # area element [L^2]
        'dV': v.L**3,                    # volume element [L^3]

        # Specific energy in drainage terms
        'eps_spec': v.L**2/(v.T**2),     # specific energy [L^2 T^-2]

        # Force and momentum terms
        'f_ext': v.M/(v.L**2*v.T**2),    # external force density [M L^-2 T^-2]
        'v_sink': v.L/v.T,               # sink velocity [L T^-1]

        # Quantum potential (if not already present)
        'Q': v.L**2/(v.T**2),            # quantum potential [L^2 T^-2]

        # Quantum work term
        'Pi_Q': v.M/(v.L*v.T**3),        # quantum work power density [M L^-1 T^-3]

        # Perturbation density
        'delta_rho': v.M/(v.L**3),       # density perturbation [M L^-3]

        # 4D coordinate
        'w': v.L,                        # 4th spatial dimension [L]

        # Additional symbols for new tests (only add if not in helper.py)
        'xi_c': v.L,                     # Healing length [L] (alias for xi)
        'chi': 1/v.L                     # Window function [L^-1]
    }, allow_overwrite=True)

    # A) Mass continuity with drainage (3D slice)
    v.info("\n--- A) Mass continuity with drainage (3D slice) ---")
    v.section("Mass continuity with drainage (3D)")
    test_mass_continuity_local_form(v)
    test_mass_continuity_linearized_form(v)
    test_global_mass_balance(v)

    # B) Momentum balance with drainage and body forces
    v.info("\n--- B) Momentum balance with drainage and body forces ---")
    v.section("Momentum balance with drainage and body forces")
    test_momentum_balance_all_terms(v)
    test_global_momentum_balance(v)

    # C) Energy conservation with drainage
    v.info("\n--- C) Energy conservation with drainage ---")
    v.section("Energy conservation with drainage")
    test_energy_density_components(v)
    test_energy_flux_terms(v)
    test_local_energy_balance(v)
    test_global_energy_balance(v)

    # D) 4D<->3D consistency
    v.info("\n--- D) 4D<->3D consistency ---")
    v.section("4D->3D drainage consistency")
    test_4d_3d_consistency(v)

    # E) Noether-style current naming and sanity checks
    v.info("\n--- E) Noether currents and sanity checks ---")
    v.section("Noether currents and diagnostics")
    test_noether_current_naming(v)
    test_sanity_reductions_and_diagnostics(v)

    # F) New framework-specific tests (added to match updated documentation)
    v.info("\n--- F) Framework-specific formulas and identities ---")
    v.section("Specific formulas from updated framework")
    test_drainage_rate_formula(v)
    test_energy_barrier_formula(v)
    test_quantum_power_identity(v)
    test_enhanced_4d_3d_matching(v)

    # Final summary
    v.summary()


if __name__ == "__main__":
    test_conservation_laws_and_aether_drainage()
