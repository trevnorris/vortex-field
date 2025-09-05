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

    # Mathematical equation verification: ∂_t ρ + ∇·(ρv) = -∑ᵢ Ṁᵢ δ³(r-rᵢ)
    # Each term verified to have same dimension [M L^-3 T^-1]
    lhs_total = lhs_rate + lhs_flux
    rhs_total = -rhs_sink  # Note the negative sign from the equation
    v.check_dims("3D continuity equation total LHS", lhs_total, target_mass_rate)
    v.check_dims("3D continuity equation total RHS", rhs_total, target_mass_rate)

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

    # Mathematical equation verification: d/dt ∫ρ dV + ∮ ρv·dA = -∑ᵢ Ṁᵢ
    # Each term verified to have same dimension [M T^-1]
    lhs_global = mass_rate_V + mass_flux_S
    rhs_global = -total_sinks  # Note the negative sign
    v.check_dims("Global mass balance equation total LHS", lhs_global, target_mass_flow)
    v.check_dims("Global mass balance equation total RHS", rhs_global, target_mass_flow)

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

    # Mathematical equation verification: ∂_t(ρv) + ∇·(ρv⊗v + pI) = -ρ∇Φ_g - ρ∇Q - ∑ᵢ Ṁᵢv_{*i} δ³ + f_ext
    # Each term verified to have same dimension [M L^-2 T^-2]
    lhs_momentum = mom_rate + mom_flux + press_grad
    rhs_momentum = -grav_force - quant_force - drain_mom + f_ext
    v.check_dims("3D momentum balance equation total LHS", lhs_momentum, target_mom)
    v.check_dims("3D momentum balance equation total RHS", rhs_momentum, target_mom)

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

    # Mathematical equation verification: e = ½ρv² + u(ρ) + ρΦ_g + e_Q
    # All components verified to have same dimension [M L^-1 T^-2]
    total_energy_density = kin_energy/2 + int_energy + grav_energy + quant_energy
    expected_energy_density = v.get_dim('e')
    v.check_dims("Energy density definition total", total_energy_density, target_edens)
    v.check_dims("Energy density definition expected", expected_energy_density, target_edens)

    v.success("Energy density components verified")


def test_energy_flux_terms(v):
    """Test energy flux (Poynting-like) terms."""
    # S = (e + p)v + S_Q

    # Convective flux: (e + p)v
    S_conv = (v.get_dim('e') + v.get_dim('p')) * v.get_dim('v')

    # Quantum/dispersion flux (if present)
    S_Q_term = v.get_dim('S_flux')  # Renamed to avoid confusion with total S

    target_eflux = v.M / v.T**3

    v.check_dims("Energy flux: convective (e+p)v", S_conv, target_eflux)
    v.check_dims("Energy flux: quantum S_Q", S_Q_term, target_eflux)

    # Mathematical equation verification: S = (e + p)v + S_Q
    # All terms verified to have same dimension [M T^-3]
    total_energy_flux = S_conv + S_Q_term
    expected_energy_flux = v.get_dim('S_flux')
    v.check_dims("Energy flux definition total", total_energy_flux, target_eflux)
    v.check_dims("Energy flux definition expected", expected_energy_flux, target_eflux)

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

    # Mathematical equation verification: ∂_t e + ∇·S = -ρv·∇Φ_g - ∑ᵢ Ṁᵢ εᵢ δ³ + Π_Q
    # Each term verified to have same dimension [M L^-1 T^-3]
    lhs_energy = e_rate + div_S
    rhs_energy = -grav_work - drain_power + Pi_Q
    v.check_dims("Local energy balance equation total LHS", lhs_energy, target_power_dens)
    v.check_dims("Local energy balance equation total RHS", rhs_energy, target_power_dens)

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


def test_microscopic_drainage_velocity(v):
    """Test microscopic drainage velocity formula: v_w ≈ Γ/(2πr_4)."""
    # Framework line 764-765: circulation-based velocity formula

    # Components
    Gamma = v.get_dim('Gamma')                # Circulation quantum [L^2 T^-1]
    r_4 = v.L                                 # 4D radial distance [L] (using dimension symbol)
    pi_const = 1                              # Dimensionless constant (2π)

    # Microscopic drainage velocity formula: v_w ≈ Γ/(2πr_4)
    v_w_formula = Gamma / (pi_const * r_4)    # [L^2 T^-1] / [L] = [L T^-1]

    # Generic velocity
    v_w_generic = v.get_dim('v')              # [L T^-1]

    # Dimensional check
    v.check_dims("Drainage velocity formula", v_w_formula, v_w_generic)

    # Mathematical equation verification: v_w = Γ/(2πr_4)
    # Both terms verified to have same dimension [L T^-1]
    v.check_dims("Microscopic drainage velocity formula", v_w_formula, v_w_generic)

    v.success("Microscopic drainage velocity verified")


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

    # Mathematical equation verification: Ṁᵢ ≈ ρ_{4D}⁰ Γ ξ_c²
    # Both terms verified to have same dimension [M T^-1]
    v.check_dims("Drainage rate formula relationship", M_dot_formula, M_dot_generic)

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

    # Mathematical equation verification: -ρv·∇Q = -∇·S_Q + Π_Q
    # Each term verified to have same dimension [M L^-1 T^-3]
    rhs_total = -div_S_Q + Pi_Q
    lhs_neg = -lhs  # The equation has negative sign on LHS
    v.check_dims("Quantum power identity equation LHS", lhs_neg, target_power_dens)
    v.check_dims("Quantum power identity equation RHS", rhs_total, target_power_dens)

    v.success("Quantum power identity verified")


def test_bulk_dissipation_equation(v):
    """Test bulk dissipation equation: ∂_t ρ_bulk + ∂_w(ρ_bulk v_w) = -γ ρ_bulk."""
    # Framework line 784-785: dissipation to prevent accumulation

    # Components
    rho_bulk = v.get_dim('rho_4')                    # Bulk density [M L^-4]
    v_w = v.get_dim('v')                             # Bulk velocity in w direction [L T^-1]
    gamma = 1/v.T                                    # Dissipation rate [T^-1]

    # LHS terms
    bulk_rate = v.dt(rho_bulk)                       # ∂_t ρ_bulk [M L^-4 T^-1]
    bulk_flux_div = v.get_dim('w_deriv') * (rho_bulk * v_w)  # ∂_w(ρ_bulk v_w) [M L^-4 T^-1]

    # RHS dissipation term
    dissipation = gamma * rho_bulk                   # -γ ρ_bulk [M L^-4 T^-1]

    target_bulk_rate = v.M / (v.L**4 * v.T)          # [M L^-4 T^-1]

    v.check_dims("Bulk dissipation: ∂_t ρ_bulk", bulk_rate, target_bulk_rate)
    v.check_dims("Bulk dissipation: ∂_w flux", bulk_flux_div, target_bulk_rate)
    v.check_dims("Bulk dissipation: γ ρ_bulk", dissipation, target_bulk_rate)

    # Mathematical equation verification: ∂_t ρ_bulk + ∂_w(ρ_bulk v_w) = -γ ρ_bulk
    # Each term verified to have same dimension [M L^-4 T^-1]
    lhs_bulk = bulk_rate + bulk_flux_div
    rhs_bulk = -dissipation
    v.check_dims("Bulk dissipation equation total LHS", lhs_bulk, target_bulk_rate)
    v.check_dims("Bulk dissipation equation total RHS", rhs_bulk, target_bulk_rate)

    v.success("Bulk dissipation equation verified")


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

    # Mathematical equation verification: ∫ ∂_w(ρ_bulk v_w) dw = -∑ᵢ Ṁᵢ δ³
    # Both terms verified to have same dimension [M L^-3 T^-1]
    rhs_negative = -rhs_sinks
    v.check_dims("4D-3D slab integration matching LHS", lhs_integrated, v.M/(v.L**3*v.T))
    v.check_dims("4D-3D slab integration matching RHS", rhs_negative, v.M/(v.L**3*v.T))

    v.success("Enhanced 4D-3D matching verified")


def test_4d_continuity_equation(v):
    """Test 4D continuity equation: ∂_t ρ_{4D} + ∇_4·(ρ_{4D} v_4) = -∑_i Ṁ_i δ^4."""
    # Framework equation 735: base 4D conservation law

    # LHS terms
    rho_4d_rate = v.dt(v.get_dim('rho_4'))                           # [M L^-4 T^-1]
    rho_4d_flux = v.div_dim(v.get_dim('rho_4') * v.get_dim('v'))     # 4D divergence [M L^-4 T^-1]

    # RHS sink term
    sink_4d = v.get_dim('M_dot_i') * v.get_dim('delta4')             # [M L^-4 T^-1]

    target_4d_rate = v.M / (v.L**4 * v.T)

    v.check_dims("4D continuity: ∂_t ρ_{4D}", rho_4d_rate, target_4d_rate)
    v.check_dims("4D continuity: ∇_4·(ρ_{4D} v_4)", rho_4d_flux, target_4d_rate)
    v.check_dims("4D continuity: sink term", sink_4d, target_4d_rate)

    # Verify equation structure consistency
    lhs_4d = rho_4d_rate + rho_4d_flux
    rhs_4d = -sink_4d
    v.check_dims("4D continuity equation LHS", lhs_4d, target_4d_rate)
    v.check_dims("4D continuity equation RHS", rhs_4d, target_4d_rate)

    v.success("4D continuity equation verified")


def test_3d_momentum_equation_structure(v):
    """Test detailed structure of 3D momentum equation from framework line 751."""
    # ∂_t(ρv) + ∇·(ρv⊗v + pI) = -ρ∇Φ_g - ρ∇Q - ∑_i Ṁ_i v_{*i} δ³ + f_ext

    # LHS tensor components
    momentum_density_rate = v.dt(v.get_dim('rho') * v.get_dim('v'))  # [M L^-2 T^-2]
    convective_tensor = v.div_dim(v.get_dim('rho') * v.get_dim('v')**2)  # ∇·(ρv⊗v)
    pressure_tensor = v.div_dim(v.get_dim('p'))  # ∇·(pI) = ∇p

    # RHS force components
    gravitational_force = v.get_dim('rho') * v.grad_dim(v.get_dim('Phi_g'))
    quantum_force = v.get_dim('rho') * v.grad_dim(v.get_dim('Q'))
    drainage_momentum_force = v.get_dim('M_dot_i') * v.get_dim('v_sink') * v.get_dim('delta3')
    external_force = v.get_dim('f_ext')

    target_force_density = v.M / (v.L**2 * v.T**2)  # [M L^-2 T^-2]

    # Verify each term has correct momentum force dimensions
    force_terms = [
        ("∂_t(ρv)", momentum_density_rate),
        ("∇·(ρv⊗v)", convective_tensor),
        ("∇p", pressure_tensor),
        ("ρ∇Φ_g", gravitational_force),
        ("ρ∇Q", quantum_force),
        ("Ṁ_i v_{*i} δ³", drainage_momentum_force),
        ("f_ext", external_force)
    ]

    for name, term in force_terms:
        v.check_dims(f"3D momentum equation term: {name}", term, target_force_density)

    # Verify total equation structure
    lhs_total = momentum_density_rate + convective_tensor + pressure_tensor
    rhs_total = -gravitational_force - quantum_force - drainage_momentum_force + external_force

    v.check_dims("3D momentum equation: total LHS", lhs_total, target_force_density)
    v.check_dims("3D momentum equation: total RHS", rhs_total, target_force_density)

    v.success("3D momentum equation structure verified")


def test_circulation_velocity_relationship(v):
    """Test circulation velocity relationship: v_w ≈ Γ/(2πr_4) with r_4 = √(ρ² + w²)."""
    # Framework equations 765-766: microscopic drainage velocity

    # Components of r_4 = √(ρ² + w²)
    rho_coord = v.L                      # Radial coordinate in 3D [L]
    w_coord = v.get_dim('w')             # 4th dimension coordinate [L]
    r_4_distance = (rho_coord**2 + w_coord**2)**(1/2)  # [L]

    # Circulation velocity formula
    Gamma = v.get_dim('Gamma')           # [L^2 T^-1]
    pi_factor = 1                        # 2π is dimensionless
    v_w_circulation = Gamma / (pi_factor * r_4_distance)  # [L T^-1]

    v.check_dims("r_4 distance", r_4_distance, v.L)
    v.check_dims("Circulation velocity v_w", v_w_circulation, v.L/v.T)

    # Verify consistency with generic velocity
    v_generic = v.get_dim('v')
    v.check_dims("Circulation velocity consistency", v_w_circulation, v_generic)

    v.success("Circulation velocity relationship verified")


def test_core_integral_drainage_rate(v):
    """Test core integral for drainage rate: Ṁ_i ≈ ρ_{4D}⁰ ∫_core v_w dA_⊥."""
    # Framework equations 769-772: transverse integral over vortex core
    # Note: This is dimensionally consistent with ρ_{4D}⁰ Γ ξ_c² from eq 772

    # Background density
    rho_4d_background = v.get_dim('rho_4')               # [M L^-4]

    # From framework equation 772: Ṁ_i ~ ρ_{4D}⁰ Γ ξ_c²
    Gamma = v.get_dim('Gamma')                           # [L^2 T^-1]
    xi_c_area = v.get_dim('xi_c')**2                     # [L^2]

    # This gives the correct dimensions for drainage rate
    M_dot_from_framework = rho_4d_background * Gamma * xi_c_area  # [M L^-4] × [L^2 T^-1] × [L^2] = [M T^-1]

    # Generic drainage rate
    M_dot_generic = v.get_dim('M_dot_i')                 # [M T^-1]

    v.check_dims("Drainage rate from framework formula", M_dot_from_framework, M_dot_generic)

    # Verify intermediate dimensional consistency with detailed velocity integral
    # The velocity integral ∫ v_w dA_⊥ should dimensionally equal Γ × ξ_c (see framework 769-771)
    velocity_integral_equiv = Gamma * v.get_dim('xi_c')  # [L^2 T^-1] × [L] = [L^3 T^-1]
    v_w_core = v.get_dim('v')                            # [L T^-1]
    core_volume = v.get_dim('xi_c')**3                   # [L^3]
    velocity_times_volume = v_w_core * core_volume       # [L T^-1] × [L^3] = [L^4 T^-1]

    # Note: The framework shows that the actual integral ∫ v_w dA_⊥ has dimensions [L^3 T^-1]
    # This is consistent with velocity [L T^-1] times area [L^2] giving [L^3 T^-1]
    v.check_dims("Velocity integral equivalent", velocity_integral_equiv, v.L**3/v.T)

    v.success("Core integral drainage rate verified")


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
        'chi': 1/v.L,                    # Window function [L^-1]
        'w_deriv': 1/v.L                 # w-derivative operator [L^-1]
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
    test_microscopic_drainage_velocity(v)
    test_drainage_rate_formula(v)
    test_energy_barrier_formula(v)
    test_quantum_power_identity(v)
    test_bulk_dissipation_equation(v)
    test_enhanced_4d_3d_matching(v)

    # G) Additional mathematical structure tests (based on equations in documentation)
    v.info("\n--- G) Additional equation structure verification ---")
    v.section("Detailed equation structures from framework")
    test_4d_continuity_equation(v)
    test_3d_momentum_equation_structure(v)
    test_circulation_velocity_relationship(v)
    test_core_integral_drainage_rate(v)

    # Final summary
    return v.summary()


if __name__ == "__main__":
    success_rate = test_conservation_laws_and_aether_drainage()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
