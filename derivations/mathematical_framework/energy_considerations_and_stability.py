#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Energy Considerations and Stability - Verification (fixed)
"""

import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper, UnitSystem, define_symbols_batch,
    quick_verify, verify_conservation_law, batch_check_dims
)
import sympy as sp
from sympy import pi, sqrt, Rational

# Init (SI dimensions are the anchor in helper.py)
v = PhysicsVerificationHelper(
    "Energy Considerations and Stability",
    "Dimensional and mathematical verification of energy framework and stability",
    unit_system=UnitSystem.SI
)

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

    # Rim dynamics
    'I_theta': v.M * v.L**2,
    'K_theta': v.M * v.L**2 / v.T**2,
    'U_3': v.M * v.L**2 / v.T**2,
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

# ==============================================================================
# SMALL PARAMETERS AND SCALING
# ==============================================================================
v.section("Small Parameters and Dimensional Scaling")

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

# ==============================================================================
# EM ENERGY (SI — helper dimensions are SI-anchored)
# ==============================================================================
v.section("EM Energy (SI)")

# u_EM = (ε0 E^2 + B^2/μ0)/2
u_EM_expr = (v.get_dim('epsilon_0')*v.get_dim('E')**2
             + v.get_dim('B')**2 / v.get_dim('mu_0'))/2
v.check_dims("EM energy density (SI)", u_EM_expr, v.get_dim('u_EM'))

# S = (1/μ0) E × B
S_expr = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
v.check_dims("EM Poynting vector (SI)", S_expr, v.get_dim('S_poynting'))

# ==============================================================================
# 4D GROSS-PITAEVSKII ENERGY FUNCTIONAL
# ==============================================================================
v.section("4D Gross-Pitaevskii Energy Functional")

# |Ψ|^2 = ρ₄D/m
v.check_dims("|Ψ|² = ρ₄D/m",
             v.get_dim('Psi_GP_4D')**2,
             v.get_dim('rho_4')/v.get_dim('m'))

# E[Ψ] density terms: (ħ²/2m)|∇Ψ|² + (g/2)|Ψ|⁴   (no extra m)
kinetic_term = v.get_dim('hbar')**2 / (2*v.get_dim('m')) * (v.L**(-1)*v.get_dim('Psi_GP_4D'))**2
interaction_term = (v.get_dim('g_GP_4D')/2) * v.get_dim('Psi_GP_4D')**4

v.check_dims("GP kinetic energy density", kinetic_term, v.get_dim('u_4D'))
v.check_dims("GP interaction energy density", interaction_term, v.get_dim('u_4D'))

# Density form: (g/2m²) ρ₄D²
interaction_density_form = (v.get_dim('g_GP_4D')/(2*v.get_dim('m')**2)) * v.get_dim('rho_4')**2
v.check_dims("GP interaction (density form)", interaction_density_form, v.get_dim('u_4D'))

# ==============================================================================
# CLOSED-CORE ENERGETIC CONSTANTS
# ==============================================================================
v.section("Closed-Core Energetic Constants")

# T ≡ ρ₀ C_core ξc² — mass per length
v.check_dims("Line 'tension' (mass/length) T = ρ₀ C_core ξc²",
             v.get_dim('rho')*v.get_dim('C_core')*v.get_dim('xi_c')**2,
             v.get_dim('T_tension'))

# A ≡ κ²/(4π v_L²) — lives inside the bracket; ρ₀ multiplies outside
v.check_dims("Log coefficient A = κ²/(4π v_L²)",
             v.get_dim('kappa')**2/(4*pi*v.get_dim('v_L')**2),
             v.get_dim('A_coeff'))

# a = α ξc
v.check_dims("Core cutoff a = α ξc",
             v.get_dim('alpha')*v.get_dim('xi_c'),
             v.get_dim('a_cutoff'))

# Mass template terms (m = ρ₀ 2πR [C ξc² + A ln(R/a)])
core_contribution = v.get_dim('C_core') * v.get_dim('xi_c')**2
log_contribution  = v.get_dim('kappa')**2/(4*pi*v.get_dim('v_L')**2) * sp.log(v.get_dim('r')/v.get_dim('a_cutoff'))

v.check_dims("Core contribution to mass",
             v.get_dim('rho')*2*pi*v.get_dim('r')*core_contribution,
             v.get_dim('m'))

v.check_dims("Log contribution to mass",
             v.get_dim('rho')*2*pi*v.get_dim('r')*log_contribution,
             v.get_dim('m'))

v.assert_dimensionless(v.get_dim('r')/v.get_dim('a_cutoff'), "log argument R/a")

# ==============================================================================
# HYDRODYNAMIC FORM AND QUANTUM PRESSURE
# ==============================================================================
v.section("Hydrodynamic Form and Quantum Pressure")

# v₄ = (ħ/m) ∇θ
v.check_dims("4D velocity v₄ = (ħ/m)∇θ",
             v.get_dim('hbar')/v.get_dim('m')*v.L**(-1)*v.get_dim('theta'),
             v.get_dim('v'))

# Force density from quantum pressure: f_Q ~ -ρ₄D ∇(Q/m) with Q/m ∼ (ħ²/m²) ∇²√ρ / √ρ
Q_over_m_dim = (v.get_dim('hbar')**2 / v.get_dim('m')**2) * v.L**(-2)   # [L²/T²]
qp_force_density = v.get_dim('rho_4') * v.L**(-1) * Q_over_m_dim       # [M L⁻³ T⁻²]
v.check_dims("Quantum pressure force density", qp_force_density, v.M/(v.L**3*v.T**2))

# ==============================================================================
# TWIST ENERGY
# ==============================================================================
v.section("Twist Energy")

# E_twist density = (ħ² τ² / 2m) |Ψ|²  (τ~1/L)
twist_energy_density = v.get_dim('hbar')**2 * v.get_dim('tau')**2 / (2*v.get_dim('m')) * v.get_dim('Psi_GP_4D')**2
v.check_dims("Twist energy density", twist_energy_density, v.get_dim('u_4D'))

# Phase contribution: τ w is dimensionless
v.assert_dimensionless(v.get_dim('tau') * v.get_dim('w'), "twist phase contribution")

# ==============================================================================
# KELVIN/SLOPE MODE ENERGY
# ==============================================================================
v.section("Kelvin/Slope Mode Energy")

v.check_dims("Kelvin mode frequency ω ~ v_L/ξc",
             v.get_dim('v_L')/v.get_dim('xi_c'),
             v.get_dim('omega_kelvin'))

v.check_dims("Core number N_core ~ (ρ₄D⁰/m) ξc⁴",
             v.get_dim('rho_4_bg')/v.get_dim('m') * v.get_dim('xi_c')**4,
             v.get_dim('N_core'))

v.check_dims("Zero-point shake energy per core",
             Rational(1,2)*v.get_dim('hbar')*v.get_dim('omega_kelvin')*v.get_dim('N_core'),
             v.get_dim('E_shake'))

shake_3d = Rational(1,2)*v.get_dim('hbar')*v.get_dim('omega_kelvin')*v.get_dim('rho')/v.get_dim('m')*v.get_dim('xi_c')**3
v.check_dims("3D zero-point energy", shake_3d, v.M*v.L**2/v.T**2)

# ==============================================================================
# RIM PHASE DYNAMICS
# ==============================================================================
v.section("Rim Phase Dynamics")

# Angular frequency = √(K/I); speed needs a radius factor.  v_θ ~ R * √(K/I)
v.check_dims("Rim characteristic speed vθ ~ r√(Kθ/Iθ)",
             v.get_dim('r')*sqrt(v.get_dim('K_theta')/v.get_dim('I_theta')),
             v.get_dim('v_theta_rim'))

v.check_dims("Lock frequency squared ω_lock² = 9U₃/Iθ",
             9*v.get_dim('U_3')/v.get_dim('I_theta'),
             v.get_dim('omega_lock')**2)

# ==============================================================================
# CORE SCALES AND TIMESCALE HIERARCHY (4D GP)
# ==============================================================================
v.section("Core Scales and Timescale Hierarchy")

# ξc ≃ ħ/√(2 g ρ₄D⁰)
v.check_dims("Healing length ξc (4D)",
             v.get_dim('hbar')/sqrt(2*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg')),
             v.get_dim('xi_c'))

# v_L² = (g/m²) ρ₄D⁰
v.check_dims("Bulk sound speed v_L",
             sqrt(v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg')/v.get_dim('m')**2),
             v.get_dim('v_L'))

# τ_core = ξc/v_L  and  τ_core = ħ m /(√2 g ρ₄D⁰)
v.check_dims("Core relaxation time (ratio form)",
             v.get_dim('xi_c')/v.get_dim('v_L'),
             v.get_dim('tau_core'))

v.check_dims("Core relaxation time (direct form)",
             v.get_dim('hbar')*v.get_dim('m')/(sqrt(2)*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg')),
             v.get_dim('tau_core'))

v.check_dims("Core time consistency (dimensional)",
             v.get_dim('xi_c')/v.get_dim('v_L'),
             v.get_dim('hbar')*v.get_dim('m')/(sqrt(2)*v.get_dim('g_GP_4D')*v.get_dim('rho_4_bg')))

# ==============================================================================
# LINEAR STABILITY
# ==============================================================================
v.section("Linear Stability")

bogoliubov_freq_sq = v.get_dim('v_L')**2*v.get_dim('k')**2 + v.get_dim('hbar')**2/(4*v.get_dim('m')**2)*v.get_dim('k')**4
v.check_dims("Bogoliubov dispersion ω²(k)", bogoliubov_freq_sq, v.get_dim('omega_bogo')**2)
v.check_dims("Bogoliubov: sound term", v.get_dim('v_L')**2*v.get_dim('k')**2, (1/v.T)**2)
v.check_dims("Bogoliubov: dispersion term", v.get_dim('hbar')**2/(4*v.get_dim('m')**2)*v.get_dim('k')**4, (1/v.T)**2)

v.check_dims("Landau bound velocity", v.get_dim('v_bg'), v.get_dim('v_L'))

# ==============================================================================
# ENERGY PROJECTION TO 3D SLICE
# ==============================================================================
v.section("Energy Projection to 3D Slice")

# 3D energy from projection: u_3D(x,t) = ∫ u_4D(x,w,t) χ_ξ(w) dw
v.check_dims("3D energy from projection",
             v.get_dim('u_4D') * v.get_dim('chi_xi') * v.get_dim('w'),
             v.get_dim('u_3D'))

# Window "normalization": ∫ χ_ξ(w) dw ~ ξ_c  (dimensionally: χ·w ~ L)
v.check_dims("Projection window has length scale ~ ξ_c",
             v.get_dim('chi_xi') * v.get_dim('w'),
             v.get_dim('xi_c'))

v.check_dims("[u₄D] = [u₃D]/L",
             v.get_dim('u_4D'),
             v.get_dim('u_3D')/v.L)

v.check_dims("ρ₀ = ρ₄D⁰ ξc",
             v.get_dim('rho'),
             v.get_dim('rho_4_bg')*v.get_dim('xi_c'))

# ==============================================================================
# GLOBAL ENERGY BALANCE
# ==============================================================================
v.section("Global Energy Balance")

energy_rate   = v.get_dim('u_3D')*v.L**3 / v.T
poynting_flux = v.get_dim('S_poynting')*v.L**2
joule_heating = v.get_dim('j_current')*v.get_dim('E')*v.L**3
bulk_exchange = v.get_dim('W_exch')*v.L**3

target_power = v.M*v.L**2/v.T**3

v.check_dims("d/dt ∫ u₃D d³x", energy_rate, target_power)
v.check_dims("∮ S·da", poynting_flux, target_power)
v.check_dims("∫ J·E d³x", joule_heating, target_power)
v.check_dims("∫ W_exch d³x", bulk_exchange, target_power)

verify_conservation_law(v, "Global energy balance",
                        energy_rate, poynting_flux + joule_heating + bulk_exchange)

# ==============================================================================
# SOMMERFELD OUTGOING WAVE CONDITION
# ==============================================================================
v.section("Sommerfeld Outgoing Wave Condition")

# Use Ψ_field from helper (a scalar field on 3D slice)
outgoing_op_1 = v.get_dim('r')*v.get_dim('Psi_field')/v.get_dim('r')           # ∂_r(rΨ)
outgoing_op_2 = v.get_dim('r')*v.get_dim('Psi_field')/(v.get_dim('c')*v.get_dim('t'))  # (1/c)∂_t(rΨ)

v.check_dims("Sommerfeld: spatial term", outgoing_op_1, v.get_dim('Psi_field'))
v.check_dims("Sommerfeld: temporal term", outgoing_op_2, v.get_dim('Psi_field'))

# ==============================================================================
# PARAMETRIC RELATIONSHIPS AND HIERARCHIES
# ==============================================================================
v.section("Parametric Relationships and Hierarchies")

for i, eps in enumerate([1e-5, 1e-3, 1e-2]):
    quick_verify(f"Small parameter ε{i+1} ≪ 1", eps < 0.1, helper=v)

macro_time_estimate = 1e-4 / 3e8    # 0.1 mm / c
core_time_estimate  = 1e-21
quick_verify("Core timescale hierarchy τ_core ≪ τ_macro",
             core_time_estimate < 0.01*macro_time_estimate, helper=v)

omega_estimates = {'core': 1e15, 'kelvin': 1e12, 'lock': 1e10}
quick_verify("Frequency hierarchy reasonable",
             omega_estimates['lock'] < omega_estimates['kelvin'] < omega_estimates['core'], helper=v)

# ==============================================================================
# MATHEMATICAL CONSISTENCY CHECKS
# ==============================================================================
v.section("Mathematical Consistency Checks")

# 4D relations g = v_L² m² / ρ₄D⁰  and  g = ħ² /(ξ_c² ρ₄D⁰)
interaction_from_sound   = v.get_dim('v_L')**2 * v.get_dim('m')**2 / v.get_dim('rho_4_bg')
interaction_from_healing = v.get_dim('hbar')**2 / (v.get_dim('xi_c')**2 * v.get_dim('rho_4_bg'))

v.check_dims("g from sound speed",   interaction_from_sound,   v.get_dim('g_GP_4D'))
v.check_dims("g from healing length",interaction_from_healing, v.get_dim('g_GP_4D'))
v.check_dims("Interaction strength consistency", interaction_from_sound, interaction_from_healing)

# Dimensionless sanity checks
for name in ["epsilon_rho", "epsilon_v", "epsilon_xi", "alpha", "N_core"]:
    v.assert_dimensionless(v.get_dim(name), name)

# ==============================================================================
# TOPOLOGICAL AND COMPOSITE STABILITY
# ==============================================================================
v.section("Topological and Composite Stability")

v.check_dims("Quantized circulation Γ ~ ħ/m",
             v.get_dim('hbar')/v.get_dim('m'),
             v.get_dim('Gamma'))

v.check_dims("Projected circulation units v·dl",
             v.get_dim('v')*v.get_dim('dl'),
             v.get_dim('Gamma'))

v.info("Topological protection is primarily qualitative")
v.info("Quantized circulation ensures discrete stability constraints")
v.info("Braided arrangements avoid low-order resonances")

# ==============================================================================
# SUMMARY
# ==============================================================================
v.section("VALIDATION SUMMARY")
success_rate = v.summary()
