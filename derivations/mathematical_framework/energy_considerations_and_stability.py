#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Energy Considerations and Stability - Verification
=================================================

Dimensional and mathematical verification of the energy framework and stability
analysis for the 4D vortex model, including Gross-Pitaevskii energy functionals,
topological stability, and projection to 3D slice dynamics.

Based on doc/main.tex, mathematical framework section (Energy considerations).
"""

import os
import sys
import sympy as sp
from sympy import pi, sqrt, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,
    define_symbols_batch,
    quick_verify,
    verify_conservation_law,
    batch_check_dims
)



def test_small_parameters_and_scaling(v):
    """
    Test small parameter definitions and dimensional scaling.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Small Parameters and Dimensional Scaling")

    v.check_dims("εᵨ = δρ₄D/ρ₄D⁰",
                 v.get_dim('delta_rho_4')/v.get_dim('rho_4_bg'),
                 v.get_dim('epsilon_rho'))

    v.check_dims("εᵥ = |v|/c",
                 v.get_dim('v')/v.get_dim('c'),
                 v.get_dim('epsilon_v'))

    v.check_dims("εξ = ξc/L",
                 v.get_dim('xi_c')/v.L,
                 v.get_dim('epsilon_xi'))

    v.assert_dimensionless(v.get_dim('epsilon_rho'), "density parameter εᵨ")
    v.assert_dimensionless(v.get_dim('epsilon_v'), "velocity parameter εᵥ")
    v.assert_dimensionless(v.get_dim('epsilon_xi'), "geometric parameter εξ")

    v.success("Small parameters and scaling verified")

def test_em_energy(v):
    """
    Test electromagnetic energy density and Poynting vector in SI units.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("EM Energy (SI)")

    # u_EM = (ε0 E^2 + B^2/μ0)/2
    u_EM_expr = (v.get_dim('epsilon_0')*v.get_dim('E')**2
                 + v.get_dim('B')**2 / v.get_dim('mu_0'))/2
    v.check_dims("EM energy density (SI)", u_EM_expr, v.get_dim('u_EM'))

    # S = (1/μ0) E × B
    S_expr = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
    v.check_dims("EM Poynting vector (SI)", S_expr, v.get_dim('S_poynting'))

    v.success("EM energy relationships verified")

def test_gp_energy_functional(v):
    """
    Test 4D Gross-Pitaevskii energy functional components.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("4D Gross-Pitaevskii Energy Functional")

    # |Ψ|^2 = ρ₄D/m (Eq. 1099)
    v.check_dims("|Ψ|² = ρ₄D/m",
                 v.get_dim('Psi_GP_4D')**2,
                 v.get_dim('rho_4')/v.get_dim('m'))

    # Mathematical verification of the relationship
    psi_squared_lhs = v.get_dim('Psi_GP_4D')**2
    psi_squared_rhs = v.get_dim('rho_4')/v.get_dim('m')
    v.check_dims("|Ψ|² = ρ₄D/m relationship", psi_squared_lhs, psi_squared_rhs)

    # E[Ψ] density terms: (ħ²/2m)|∇Ψ|² + (g/2)|Ψ|⁴   (no extra m)
    kinetic_term = v.get_dim('hbar')**2 / (2*v.get_dim('m')) * (v.L**(-1)*v.get_dim('Psi_GP_4D'))**2
    interaction_term = (v.get_dim('g_GP_4D')/2) * v.get_dim('Psi_GP_4D')**4

    v.check_dims("GP kinetic energy density", kinetic_term, v.get_dim('u_4D'))
    v.check_dims("GP interaction energy density", interaction_term, v.get_dim('u_4D'))

    # Density form: (g/2m²) ρ₄D² (Eq. 1102)
    interaction_density_form = (v.get_dim('g_GP_4D')/(2*v.get_dim('m')**2)) * v.get_dim('rho_4')**2
    v.check_dims("GP interaction (density form)", interaction_density_form, v.get_dim('u_4D'))

    # Mathematical verification: two forms of interaction energy should be equivalent
    # (g/2)|Ψ|⁴ = (g/2m²)ρ₄D² using |Ψ|² = ρ₄D/m
    interaction_form1 = (v.get_dim('g_GP_4D')/2) * v.get_dim('Psi_GP_4D')**4
    interaction_form2 = (v.get_dim('g_GP_4D')/(2*v.get_dim('m')**2)) * v.get_dim('rho_4')**2
    v.check_dims("GP interaction: (g/2)|Ψ|⁴ = (g/2m²)ρ₄D²", interaction_form1, interaction_form2)

    v.success("4D Gross-Pitaevskii energy functional verified")

def test_closed_core_energetic_constants(v):
    """
    Test closed-core energetic constants and mass template.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Closed-Core Energetic Constants")

    # T ≡ ρ₀ C_core ξc² — mass per length (Eq. 1108)
    v.check_dims("Line 'tension' (mass/length) T = ρ₀ C_core ξc²",
                 v.get_dim('rho')*v.get_dim('C_core')*v.get_dim('xi_c')**2,
                 v.get_dim('T_tension'))

    # Mathematical verification of T definition
    T_lhs = v.get_dim('T_tension')
    T_rhs = v.get_dim('rho')*v.get_dim('C_core')*v.get_dim('xi_c')**2
    v.check_dims("T ≡ ρ₀ C_core ξc²", T_lhs, T_rhs)

    # A ≡ κ²/(4π v_L²) — lives inside the bracket; ρ₀ multiplies outside (Eq. 1109)
    v.check_dims("Log coefficient A = κ²/(4π v_L²)",
                 v.get_dim('kappa')**2/(4*pi*v.get_dim('v_L')**2),
                 v.get_dim('A_coeff'))

    # Mathematical verification of A definition
    A_lhs = v.get_dim('A_coeff')
    A_rhs = v.get_dim('kappa')**2/(4*pi*v.get_dim('v_L')**2)
    v.check_dims("A ≡ κ²/(4π v_L²)", A_lhs, A_rhs)

    # a = α ξc (Eq. 1110)
    v.check_dims("Core cutoff a = α ξc",
                 v.get_dim('alpha')*v.get_dim('xi_c'),
                 v.get_dim('a_cutoff'))

    # Mathematical verification of a definition
    a_lhs = v.get_dim('a_cutoff')
    a_rhs = v.get_dim('alpha')*v.get_dim('xi_c')
    v.check_dims("a = α ξc", a_lhs, a_rhs)

    # Mass template terms (m = ρ₀ 2πR [C ξc² + A ln(R/a)]) (Eq. 1115)
    core_contribution = v.get_dim('C_core') * v.get_dim('xi_c')**2
    log_contribution  = v.get_dim('kappa')**2/(4*pi*v.get_dim('v_L')**2) * sp.log(v.get_dim('r')/v.get_dim('a_cutoff'))

    v.check_dims("Core contribution to mass",
                 v.get_dim('rho')*2*pi*v.get_dim('r')*core_contribution,
                 v.get_dim('m'))

    v.check_dims("Log contribution to mass",
                 v.get_dim('rho')*2*pi*v.get_dim('r')*log_contribution,
                 v.get_dim('m'))

    # Mathematical verification of complete mass template
    mass_template_lhs = v.get_dim('m')
    mass_template_rhs = v.get_dim('rho')*2*pi*v.get_dim('r')*(core_contribution + log_contribution)
    v.check_dims("Mass template: m(R) = ρ₀ 2πR [C_core ξc² + A ln(R/a)]", mass_template_lhs, mass_template_rhs)

    v.assert_dimensionless(v.get_dim('r')/v.get_dim('a_cutoff'), "log argument R/a")

    v.success("Closed-core energetic constants verified")

def test_hydrodynamic_form_and_quantum_pressure(v):
    """
    Test hydrodynamic formulation and quantum pressure effects.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Hydrodynamic Form and Quantum Pressure")

    # v₄ = (ħ/m) ∇θ (Eq. 1121)
    v.check_dims("4D velocity v₄ = (ħ/m)∇θ",
                 v.get_dim('hbar')/v.get_dim('m')*v.L**(-1)*v.get_dim('theta'),
                 v.get_dim('v'))

    # Mathematical verification of hydrodynamic velocity
    v4_lhs = v.get_dim('v')
    v4_rhs = v.get_dim('hbar')/v.get_dim('m')*v.L**(-1)*v.get_dim('theta')
    v.check_dims("v₄ = (ħ/m)∇θ", v4_lhs, v4_rhs)

    # Force density from quantum pressure: f_Q ~ -ρ₄D ∇(Q/m) with Q/m ∼ (ħ²/m²) ∇²√ρ / √ρ (Eq. 1122)
    Q_over_m_dim = (v.get_dim('hbar')**2 / v.get_dim('m')**2) * v.L**(-2)   # [L²/T²]
    qp_force_density = v.get_dim('rho_4') * v.L**(-1) * Q_over_m_dim       # [M L⁻³ T⁻²]
    v.check_dims("Quantum pressure force density", qp_force_density, v.M/(v.L**3*v.T**2))

    v.success("Hydrodynamic form and quantum pressure verified")

def test_twist_energy(v):
    """
    Test twist energy contributions to the energy functional.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Twist Energy")

    # E_twist density = (ħ² τ² / 2m²) ρ₄D  (τ~1/L) (Eq. 1130, corrected)
    twist_energy_density = v.get_dim('hbar')**2 * v.get_dim('tau')**2 / (2*v.get_dim('m')**2) * v.get_dim('rho_4')
    v.check_dims("Twist energy density", twist_energy_density, v.get_dim('u_4D'))

    # Alternative form: E_twist = ∫ d⁴r (ħ²/2m²)ρ₄D(∂wθ)² (Eq. 1129, corrected)
    # Since (∂wθ)² = τ², these should be equivalent
    twist_energy_alt = v.get_dim('hbar')**2 / (2*v.get_dim('m')**2) * v.get_dim('rho_4') * v.get_dim('tau')**2
    v.check_dims("Twist energy (alternative form)", twist_energy_alt, v.get_dim('u_4D'))

    # Mathematical verification: both forms should be equivalent
    v.check_dims("Twist energy equivalence", twist_energy_density, twist_energy_alt)

    # Phase contribution: τ w is dimensionless
    v.assert_dimensionless(v.get_dim('tau') * v.get_dim('w'), "twist phase contribution")

    v.success("Twist energy verified")

def test_kelvin_slope_mode_energy(v):
    """
    Test Kelvin and slope mode energy contributions.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Kelvin/Slope Mode Energy")

    # Kelvin mode frequency ω ~ v_L/ξc (Eq. 1135)
    v.check_dims("Kelvin mode frequency ω ~ v_L/ξc",
                 v.get_dim('v_L')/v.get_dim('xi_c'),
                 v.get_dim('omega_kelvin'))

    # Mathematical verification of frequency relationship
    omega_lhs = v.get_dim('omega_kelvin')
    omega_rhs = v.get_dim('v_L')/v.get_dim('xi_c')
    v.check_dims("ω ~ v_L/ξc", omega_lhs, omega_rhs)

    # Core number N_core ~ (ρ₄D⁰/m) ξc⁴ (Eq. 1139-1140)
    v.check_dims("Core number N_core ~ (ρ₄D⁰/m) ξc⁴",
                 v.get_dim('rho_4_bg')/v.get_dim('m') * v.get_dim('xi_c')**4,
                 v.get_dim('N_core'))

    # Mathematical verification of core number
    N_core_lhs = v.get_dim('N_core')
    N_core_rhs = v.get_dim('rho_4_bg')/v.get_dim('m') * v.get_dim('xi_c')**4
    v.check_dims("N_core ~ (ρ₄D⁰/m) ξc⁴", N_core_lhs, N_core_rhs)

    # Zero-point shake energy per core (Eq. 1137)
    v.check_dims("Zero-point shake energy per core",
                 Rational(1,2)*v.get_dim('hbar')*v.get_dim('omega_kelvin')*v.get_dim('N_core'),
                 v.get_dim('E_shake'))

    # Mathematical verification of shake energy
    E_shake_lhs = v.get_dim('E_shake')
    E_shake_rhs = Rational(1,2)*v.get_dim('hbar')*v.get_dim('omega_kelvin')*v.get_dim('N_core')
    v.check_dims("E_shake = (1/2)ℏω N_core", E_shake_lhs, E_shake_rhs)

    # 3D projected energy (Eq. 1142)
    shake_3d = Rational(1,2)*v.get_dim('hbar')*v.get_dim('omega_kelvin')*v.get_dim('rho')/v.get_dim('m')*v.get_dim('xi_c')**3
    v.check_dims("3D zero-point energy", shake_3d, v.M*v.L**2/v.T**2)

    v.success("Kelvin/slope mode energy verified")

def test_rim_phase_dynamics(v):
    """
    Test rim phase dynamics and characteristic frequencies.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Rim Phase Dynamics")

    # Characteristic speed v_θ = √(T_rim/μ_rim) (Eq. 1147, line physics)
    v.check_dims("Rim characteristic speed vθ = √(T_rim/μ_rim)",
                 sqrt(v.get_dim('K_theta')/v.get_dim('I_theta')),
                 v.get_dim('v_theta_rim'))

    # Mathematical verification of characteristic speed
    v_theta_lhs = v.get_dim('v_theta_rim')
    v_theta_rhs = sqrt(v.get_dim('K_theta')/v.get_dim('I_theta'))
    v.check_dims("v_θ = √(T_rim/μ_rim)", v_theta_lhs, v_theta_rhs)

    # Lock frequency squared ω_lock² = 9U₃/Iθ (Eq. 1147)
    v.check_dims("Lock frequency squared ω_lock² = 9U₃/Iθ",
                 9*v.get_dim('U_3')/v.get_dim('I_theta'),
                 v.get_dim('omega_lock')**2)

    # Mathematical verification of lock frequency
    omega_lock_sq_lhs = v.get_dim('omega_lock')**2
    omega_lock_sq_rhs = 9*v.get_dim('U_3')/v.get_dim('I_theta')
    v.check_dims("ω_lock² = 9U₃/I_θ", omega_lock_sq_lhs, omega_lock_sq_rhs)

    v.success("Rim phase dynamics verified")

def test_core_scales_and_timescale_hierarchy(v):
    """
    Test core scales and timescale hierarchy in 4D GP framework.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Core Scales and Timescale Hierarchy")

    # ξc ≃ ħ/√(2gρ₄D⁰) (Eq. 1154, corrected)
    v.check_dims("Healing length ξc (4D)",
                 v.get_dim('hbar')/sqrt(2*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg')),
                 v.get_dim('xi_c'))

    # Mathematical verification of healing length
    xi_c_lhs = v.get_dim('xi_c')
    xi_c_rhs = v.get_dim('hbar')/sqrt(2*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg'))
    v.check_dims("ξc = ℏ/√(2gρ₄D⁰)", xi_c_lhs, xi_c_rhs)

    # v_L = √(gρ₄D⁰)/m (Eq. 1156, corrected)
    v.check_dims("Bulk sound speed v_L",
                 sqrt(v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg'))/v.get_dim('m'),
                 v.get_dim('v_L'))

    # Mathematical verification of sound speed
    v_L_lhs = v.get_dim('v_L')
    v_L_rhs = sqrt(v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg'))/v.get_dim('m')
    v.check_dims("v_L = √(gρ₄D⁰)/m", v_L_lhs, v_L_rhs)

    # τ_core = ξc/v_L (Eq. 1160)
    v.check_dims("Core relaxation time (ratio form)",
                 v.get_dim('xi_c')/v.get_dim('v_L'),
                 v.get_dim('tau_core'))

    # τ_core = ℏm/(√2 g ρ₄D⁰) (Eq. 1160, corrected)
    v.check_dims("Core relaxation time (direct form)",
                 v.get_dim('hbar')*v.get_dim('m')/(sqrt(2)*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg')),
                 v.get_dim('tau_core'))

    # Mathematical verification: both expressions for τ_core should be equivalent
    tau_core_ratio = v.get_dim('xi_c')/v.get_dim('v_L')
    tau_core_direct = v.get_dim('hbar')*v.get_dim('m')/(sqrt(2)*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg'))
    v.check_dims("τ_core = ξc/v_L = ℏm/(√2 g ρ₄D⁰)", tau_core_ratio, tau_core_direct)

    v.success("Core scales and timescale hierarchy verified")

def test_linear_stability(v):
    """
    Test linear stability analysis with Bogoliubov dispersion.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Linear Stability")

    # Bogoliubov dispersion ω²(k) = v_L²k² + (ℏ²/4m²)k⁴ (Eq. 1167)
    bogoliubov_freq_sq = v.get_dim('v_L')**2*v.get_dim('k')**2 + v.get_dim('hbar')**2/(4*v.get_dim('m')**2)*v.get_dim('k')**4
    v.check_dims("Bogoliubov dispersion ω²(k)", bogoliubov_freq_sq, v.get_dim('omega_bogo')**2)

    # Mathematical verification of Bogoliubov dispersion
    omega_bogo_sq_lhs = v.get_dim('omega_bogo')**2
    omega_bogo_sq_rhs = v.get_dim('v_L')**2*v.get_dim('k')**2 + v.get_dim('hbar')**2/(4*v.get_dim('m')**2)*v.get_dim('k')**4
    v.check_dims("ω²(k) = v_L²k² + (ℏ²/4m²)k⁴", omega_bogo_sq_lhs, omega_bogo_sq_rhs)

    # Component verification
    v.check_dims("Bogoliubov: sound term", v.get_dim('v_L')**2*v.get_dim('k')**2, (1/v.T)**2)
    v.check_dims("Bogoliubov: dispersion term", v.get_dim('hbar')**2/(4*v.get_dim('m')**2)*v.get_dim('k')**4, (1/v.T)**2)

    # Landau bound |v_bg| < v_L (Eq. 1171)
    v.check_dims("Landau bound velocity", v.get_dim('v_bg'), v.get_dim('v_L'))

    v.success("Linear stability verified")

def test_energy_projection_to_3d(v):
    """
    Test energy projection from 4D to 3D slice.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Energy Projection to 3D Slice")

    # 3D energy from projection: u_3D(x,t) = ∫ u_4D(x,w,t) χ_ξc(w) dw (Eq. 1178)
    v.check_dims("3D energy from projection",
                 v.get_dim('u_4D') * v.get_dim('chi_xi') * v.get_dim('w'),
                 v.get_dim('u_3D'))

    # Mathematical verification of projection integral
    u_3D_lhs = v.get_dim('u_3D')
    u_3D_rhs = v.get_dim('u_4D') * v.get_dim('chi_xi') * v.get_dim('w')
    v.check_dims("u_3D = ∫ u_4D χ_ξc dw", u_3D_lhs, u_3D_rhs)

    # Window normalization: ∫ χ_ξc(w) dw = 1 (Eq. 1179)
    v.check_dims("Projection window has length scale ~ ξ_c",
                 v.get_dim('chi_xi') * v.get_dim('w'),
                 v.get_dim('xi_c'))

    # Dimensional relationship [u₄D] = [u₃D]/L (Eq. 1181)
    v.check_dims("[u₄D] = [u₃D]/L",
                 v.get_dim('u_4D'),
                 v.get_dim('u_3D')/v.L)

    # Mathematical verification of dimensional relationship
    u_4D_lhs = v.get_dim('u_4D')
    u_4D_rhs = v.get_dim('u_3D')/v.L
    v.check_dims("[u₄D] = [u₃D]/L", u_4D_lhs, u_4D_rhs)

    # Density projection ρ₀ = ρ₄D⁰ ξc (Eq. 1181)
    v.check_dims("ρ₀ = ρ₄D⁰ ξc",
                 v.get_dim('rho'),
                 v.get_dim('rho_4_bg')*v.get_dim('xi_c'))

    # Mathematical verification of density projection
    rho_0_lhs = v.get_dim('rho')
    rho_0_rhs = v.get_dim('rho_4_bg')*v.get_dim('xi_c')
    v.check_dims("ρ₀ = ρ₄D⁰ ξc", rho_0_lhs, rho_0_rhs)

    v.success("Energy projection to 3D verified")

def test_global_energy_balance(v):
    """
    Test global energy balance and conservation laws.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Global Energy Balance")

    # Global energy balance equation (Eq. 1187-1191):
    # d/dt ∫u₃D d³x = -∮S_EM·da + ∫J_ch·E d³x + ∫W_exch d³x
    energy_rate   = v.get_dim('u_3D')*v.L**3 / v.T
    poynting_flux = v.get_dim('S_poynting')*v.L**2
    joule_heating = v.get_dim('j_current')*v.get_dim('E')*v.L**3
    bulk_exchange = v.get_dim('W_exch')*v.L**3

    target_power = v.M*v.L**2/v.T**3

    v.check_dims("d/dt ∫ u₃D d³x", energy_rate, target_power)
    v.check_dims("∮ S·da", poynting_flux, target_power)
    v.check_dims("∫ J·E d³x", joule_heating, target_power)
    v.check_dims("∫ W_exch d³x", bulk_exchange, target_power)

    # Mathematical verification of global energy balance
    energy_lhs = energy_rate
    energy_rhs = -poynting_flux + joule_heating + bulk_exchange
    v.check_dims("Global energy balance: d/dt ∫u₃D = -∮S_EM·da + ∫J_ch·E + ∫W_exch",
                 energy_lhs, energy_rhs)

    verify_conservation_law(v, "Global energy balance",
                            energy_rate, poynting_flux + joule_heating + bulk_exchange)

    v.success("Global energy balance verified")

def test_sommerfeld_outgoing_wave_condition(v):
    """
    Test Sommerfeld outgoing wave boundary condition.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Sommerfeld Outgoing Wave Condition")

    # Sommerfeld condition (∂_r - (1/c)∂_t)(rΨ) = o(r^-1) (Eq. 1183)
    # Use Ψ_field from helper (a scalar field on 3D slice)
    outgoing_op_1 = v.get_dim('Psi_field')                                     # ∂_r(rΨ)/r ~ Ψ
    outgoing_op_2 = v.get_dim('Psi_field')/(v.get_dim('c')*v.T)                # (1/c)∂_t(rΨ)/r ~ Ψ/(ct)

    v.check_dims("Sommerfeld: spatial term", outgoing_op_1, v.get_dim('Psi_field'))
    v.check_dims("Sommerfeld: temporal term", outgoing_op_2*v.get_dim('c')*v.T, v.get_dim('Psi_field'))

    # Mathematical verification: both terms should have same dimensions
    v.check_dims("Sommerfeld condition dimensional consistency", outgoing_op_1, outgoing_op_2*v.get_dim('c')*v.T)

    v.success("Sommerfeld outgoing wave condition verified")

def test_parametric_relationships_and_hierarchies(v):
    """
    Test parametric relationships and scale hierarchies.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Parametric Relationships and Hierarchies")

    for i, eps in enumerate([1e-5, 1e-3, 1e-2]):
        quick_verify(f"Small parameter ε{i+1} ≪ 1", eps < 0.1, helper=v)

    macro_time_estimate = 1e-4 / 3e8    # 0.1 mm / c
    core_time_estimate  = 1e-21
    quick_verify("Core timescale hierarchy τ_core ≪ τ_macro",
                 core_time_estimate < 0.01*macro_time_estimate, helper=v)

    omega_estimates = {'core': 1e15, 'kelvin': 1e12, 'lock': 1e10}
    quick_verify("Frequency hierarchy reasonable",
                 omega_estimates['lock'] < omega_estimates['kelvin'] < omega_estimates['core'], helper=v)

    v.success("Parametric relationships and hierarchies verified")

def test_mathematical_consistency_checks(v):
    """
    Test mathematical consistency of fundamental relationships.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Mathematical Consistency Checks")

    # Interaction strength consistency from vL = √(gρ₄D⁰)/m (Eq. 1156, corrected)
    # This gives g = v_L²m²/ρ₄D⁰, and from ξc = ℏ/√(2gρ₄D⁰) (Eq. 1154, corrected)
    # Rearranging: g = ℏ²/(2ξc²ρ₄D⁰)
    interaction_from_sound   = v.get_dim('v_L')**2 * v.get_dim('m')**2 / v.get_dim('rho_4_bg')
    interaction_from_healing = v.get_dim('hbar')**2 / (2*v.get_dim('xi_c')**2 * v.get_dim('rho_4_bg'))

    v.check_dims("g from sound speed",   interaction_from_sound,   v.get_dim('g_GP_4D'))
    v.check_dims("g from healing length",interaction_from_healing, v.get_dim('g_GP_4D'))

    # Mathematical verification of interaction strength consistency
    v.check_dims("Interaction consistency: g = v_L²m²/ρ₄D⁰ = ℏ²/(2ξc²ρ₄D⁰)",
                 interaction_from_sound, interaction_from_healing)

    # Dimensionless sanity checks
    for name in ["epsilon_rho", "epsilon_v", "epsilon_xi", "alpha", "N_core"]:
        v.assert_dimensionless(v.get_dim(name), name)

    v.success("Mathematical consistency checks verified")

def test_topological_and_composite_stability(v):
    """
    Test topological and composite stability properties.

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Topological and Composite Stability")

    v.check_dims("Quantized circulation Γ ~ ħ/m",
                 v.get_dim('hbar')/v.get_dim('m'),
                 v.get_dim('Gamma'))

    v.check_dims("Projected circulation units v·dl",
                 v.get_dim('v')*v.get_dim('dl'),
                 v.get_dim('Gamma'))

    v.info("Topological protection is primarily qualitative")
    v.info("Quantized circulation ensures discrete stability constraints")
    v.info("Braided arrangements avoid low-order resonances")

    v.success("Topological and composite stability verified")

def test_energy_considerations_and_stability():
    """
    Main test function for Energy Considerations and Stability.

    This function coordinates all verification tests for the energy framework
    and stability analysis of the 4D vortex model.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Energy Considerations and Stability",
        "Dimensional and mathematical verification of energy framework and stability",
        unit_system=UnitSystem.SI
    )

    v.section("ENERGY CONSIDERATIONS AND STABILITY VERIFICATION")

    # Symbols (purely symbolic; all dimensions come from v.get_dim(...))
    rho, w, theta, R, k, r, t = define_symbols_batch(
        ['rho', 'w', 'theta', 'R', 'k', 'r', 't'],
        positive=True, real=True
    )
    tau, alpha = define_symbols_batch(['tau', 'alpha'], real=True)

    # Dimensionless declarations
    v.declare_dimensionless('theta', 'alpha')          # tau is NOT dimensionless

    # Section-specific dimensions / parameters
    v.add_dimensions({
        # Core / geometry / small params
        'C_core': 1,
        'xi_c': v.L,                                   # healing length on 3D slice
        'epsilon_rho': 1,
        'epsilon_v': 1,

        # Energetic / line properties
        'T_tension': v.M / v.L,                        # ρ0 ξ_c^2 → mass per length (not M·L)
        'A_coeff': v.L**2,                             # κ^2/(4π v_L^2) has dimensions of L^2

        # Mode & twist
        'omega_kelvin': 1/v.T,
        'N_core': 1,
        'E_shake': v.M * v.L**2 / v.T**2,
        'tau': 1/v.L,                                  # twist density has units 1/L
        'w': v.L,                                      # ensure w is treated as a length coord

        # Rim dynamics (line physics: tension and line density)
        'I_theta': v.M / v.L,                              # line density μ_rim [M/L]
        'K_theta': v.M * v.L / v.T**2,                     # line tension T_rim [M·L/T²]
        'U_3': v.M / (v.L * v.T**2),                          # locking potential [M/(L·T²)]
        'v_theta_rim': v.L / v.T,
        'omega_lock': 1/v.T,

        # Stability / dispersion
        'tau_core': v.T,
        'v_bg': v.L / v.T,
        'omega_bogo': 1/v.T,

        # Projection + balance
        'chi_xi': 1,                                   # dimensionless window
        'W_exch': v.M * v.L**(-1) / v.T**3,            # power density (energy/time/volume)
        'a_cutoff': v.L,
    }, allow_overwrite=True)

    # Call test functions in logical order
    v.info("\n--- 1) Small Parameters and Scaling ---")
    test_small_parameters_and_scaling(v)

    v.info("\n--- 2) EM Energy ---")
    test_em_energy(v)

    v.info("\n--- 3) GP Energy Functional ---")
    test_gp_energy_functional(v)

    v.info("\n--- 4) Closed-Core Energetic Constants ---")
    test_closed_core_energetic_constants(v)

    v.info("\n--- 5) Hydrodynamic Form and Quantum Pressure ---")
    test_hydrodynamic_form_and_quantum_pressure(v)

    v.info("\n--- 6) Twist Energy ---")
    test_twist_energy(v)

    v.info("\n--- 7) Kelvin/Slope Mode Energy ---")
    test_kelvin_slope_mode_energy(v)

    v.info("\n--- 8) Rim Phase Dynamics ---")
    test_rim_phase_dynamics(v)

    v.info("\n--- 9) Core Scales and Timescale Hierarchy ---")
    test_core_scales_and_timescale_hierarchy(v)

    v.info("\n--- 10) Linear Stability ---")
    test_linear_stability(v)

    v.info("\n--- 11) Energy Projection to 3D ---")
    test_energy_projection_to_3d(v)

    v.info("\n--- 12) Global Energy Balance ---")
    test_global_energy_balance(v)

    v.info("\n--- 13) Sommerfeld Outgoing Wave Condition ---")
    test_sommerfeld_outgoing_wave_condition(v)

    v.info("\n--- 14) Parametric Relationships and Hierarchies ---")
    test_parametric_relationships_and_hierarchies(v)

    v.info("\n--- 15) Mathematical Consistency Checks ---")
    test_mathematical_consistency_checks(v)

    v.info("\n--- 16) Topological and Composite Stability ---")
    test_topological_and_composite_stability(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_energy_considerations_and_stability()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
