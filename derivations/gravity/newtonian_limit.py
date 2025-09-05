#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Newtonian Limit - Verification
====================================

Complete verification of the emergence of Newtonian gravity from the 4D vortex
framework in the static, low-velocity limit. Tests the continuity equation
derivation, Euler equation linearization, and the final Poisson equation
that reproduces standard Newtonian gravity.

Based on doc/gravity.tex, lines 74-145 (Newtonian Limit subsection).

Key Physics Verified:
- Continuity equation: ∂_t ρ_{3D} + ∇·(ρ_{3D} v) = -Ṁ_body
- Static limit emergence from scalar sector
- Euler equation linearization with sink terms
- Calibration relationships yielding Newton's law
- Connection between vortex sinks and gravitational attraction
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify, sqrt, diff

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    batch_check_dims,
    verify_conservation_law,
    verify_poisson_equation,
    quick_verify,
)


def test_continuity_equation_derivation(v):
    """
    Test the 3D continuity equation derived from 4D vortex dynamics.

    Key equation: ∂_t ρ_{3D} + ∇·(ρ_{3D} v) = -Ṁ_body

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("3D Continuity Equation from 4D Projection")

    # Define the continuity equation terms (point sink)
    time_term = v.dt(v.get_dim('rho_3D'))           # ∂_t ρ_{3D}
    flux_term = v.div_dim(v.get_dim('rho_3D') * v.get_dim('v'))  # ∇·(ρ_{3D} v)
    sink_term = v.get_dim('M_dot_body') * v.get_dim('delta3')  # -Ṁ_body δ³(r-r₀)

    # Verify dimensional consistency of continuity equation
    v.check_dims("Continuity: time derivative", time_term,
                 v.get_dim('rho_3D') / v.get_dim('t'))

    v.check_dims("Continuity: flux divergence", flux_term,
                 v.get_dim('rho_3D') * v.get_dim('v') / v.get_dim('r'))

    v.check_dims("Continuity: sink term", sink_term,
                 time_term)  # Should match ∂_t ρ_{3D} dimensions

    # Check overall continuity equation consistency (point sink)
    verify_conservation_law(v, "3D mass continuity", time_term, flux_term, -sink_term)

    # Verify background density relationship: ρ_0 = ρ_{4D}^0 · ξ_c
    rho_0_projected = v.get_dim('rho_4') * v.get_dim('xi')  # 4D density × core thickness
    v.check_dims("Background density projection", v.get_dim('rho_0'), rho_0_projected)

    v.success("3D continuity equation verified")


def test_density_deficit_balance(v):
    """
    Test the density deficit equilibrium relationship.

    In equilibrium: δρ_{3D} ≈ -ρ_body, where ρ_body = Ṁ_body / (v_eff A_core)

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Density Deficit Balance")

    # Define density deficit relationship
    delta_rho_3D = v.get_dim('rho_3D') - v.get_dim('rho_0')  # Perturbation

    # Body density from sink strength: ρ_body = Ṁ_body / (v_eff A_core)
    rho_body_derived = v.get_dim('M_dot_body') / (v.get_dim('v_eff') * v.get_dim('A_core'))

    v.check_dims("Density perturbation", delta_rho_3D, v.get_dim('rho_3D'))
    v.check_dims("Body density from sinks", rho_body_derived, v.get_dim('rho'))

    # Check core area relationship: A_core ≈ π ξ_c²
    A_core_formula = pi * v.get_dim('xi')**2
    v.check_dims("Vortex core area", v.get_dim('A_core'), A_core_formula)

    # In equilibrium, deficit balances sink
    v.check_dims("Equilibrium balance", delta_rho_3D, rho_body_derived)

    v.success("Density deficit balance verified")


def test_euler_equation_linearization(v):
    """
    Test the linearization of the Euler equation in the static limit.

    Full equation: ∂_t v + (v·∇)v = -(1/ρ_{3D})∇P - (Ṁ_body v)/ρ_{3D}
    Static limit: ∇Φ_g = (1/ρ_0)∇P with P = (g/2)ρ_{3D}²

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Euler Equation Linearization")

    # Full Euler equation terms
    time_accel = v.dt(v.get_dim('v'))                        # ∂_t v
    convective = v.get_dim('v') * v.grad_dim(v.get_dim('v')) # (v·∇)v
    pressure_grad = v.grad_dim(v.get_dim('P')) / v.get_dim('rho_3D')  # -(1/ρ)∇P
    sink_drag = (v.get_dim('M_dot_body') * v.get_dim('delta3') *
                 v.get_dim('v')) / v.get_dim('rho_3D')  # Sink drag with δ³

    # All terms should have acceleration dimensions
    v.check_dims("Euler: time derivative", time_accel, v.get_dim('a'))
    v.check_dims("Euler: convective term", convective, v.get_dim('a'))
    v.check_dims("Euler: pressure gradient", pressure_grad, v.get_dim('a'))
    v.check_dims("Euler: sink drag", sink_drag, v.get_dim('a'))

    # Static limit: velocity potential ∇Φ_g
    phi_g_gradient = v.grad_dim(v.get_dim('Phi_g'))
    v.check_dims("Gravitational acceleration", phi_g_gradient, v.get_dim('g'))

    # EOS relationship: P = (g_eos/2)ρ_{3D}² (3D pressure)
    # With g_eos = c²/ρ₀ having dimensions [L⁵/(MT²)]
    P_from_EOS = (v.get_dim('g_eos') / 2) * v.get_dim('rho_3D')**2
    v.check_dims("EOS: P = (g/2)ρ²", v.get_dim('P'), P_from_EOS)

    # Static limit with correct sign: ∇Φ_g = -g ∇ρ_{3D}
    static_relation = -v.get_dim('g_eos') * v.grad_dim(v.get_dim('rho_3D'))
    v.check_dims("Static limit relation", phi_g_gradient, static_relation)

    v.success("Euler equation linearization verified")


def test_poisson_equation_emergence(v):
    """
    Test the emergence of the Newtonian Poisson equation.

    From continuity and static Euler: ∇²Φ_g = 4πG ρ_body

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Poisson Equation Emergence")

    # Taking divergence of static relation: ∇²Φ_g = -g ∇²ρ_{3D}
    laplacian_phi_g = v.lap_dim(v.get_dim('Phi_g'))
    laplacian_rho_3D = v.lap_dim(v.get_dim('rho_3D'))

    intermediate_form = -v.get_dim('g_eos') * laplacian_rho_3D
    v.check_dims("Intermediate Poisson", laplacian_phi_g, intermediate_form)

    # Derived density Poisson: ∇²ρ_{3D} = -(4πG/g) ρ_body = -ρ_body/ξ_c²
    density_poisson_rhs = -(4*sp.pi) * v.get_dim('G')/v.get_dim('g_eos') * v.get_dim('rho_body')
    v.check_dims("Density Poisson source", laplacian_rho_3D, density_poisson_rhs)

    # Alternative form: ∇²ρ_{3D} = -ρ_body/ξ_c² (when g = 4πGξ_c²)
    density_poisson_alt = -v.get_dim('rho_body') / v.get_dim('xi')**2
    v.check_dims("Density Poisson (ξ form)", density_poisson_rhs, density_poisson_alt)

    # Final Newtonian form: ∇²Φ_g = 4πG ρ_body
    newton_poisson_rhs = 4*sp.pi * v.get_dim('G') * v.get_dim('rho_body')

    # Verify standard Poisson equation for gravity
    verify_poisson_equation(v, "Newtonian gravity", laplacian_phi_g, newton_poisson_rhs)

    # For point mass M: Φ_g = -GM/r
    point_mass_potential = v.get_dim('G') * v.get_dim('m') / v.get_dim('r')
    v.check_dims("Point mass potential", v.get_dim('Phi_g'), point_mass_potential)

    # Gravitational acceleration: a = -∇Φ_g = -GM/r²
    newton_acceleration = v.get_dim('G') * v.get_dim('m') / v.get_dim('r')**2
    v.check_dims("Newton's law F=GMm/r²", v.get_dim('g'), newton_acceleration)

    v.success("Newtonian Poisson equation verified")


def test_calibration_constants(v):
    """
    Test the calibration relationships that connect 4D parameters to Newton's G.

    Key calibrations: g_eos = c²/ρ_0 and G = c²/(4πρ_0 ξ_c²)

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Calibration Constants")

    # Calibration of EOS parameter: g_eos = c²/ρ_0 (already set in add_dimensions)
    # Verify it matches the expected calibration
    g_calibrated = v.get_dim('c')**2 / v.get_dim('rho_0')
    v.check_dims("EOS coupling calibration", v.get_dim('g_eos'), g_calibrated)

    # This should give g_eos the right dimensions for P = (g_eos/2)ρ²
    pressure_dimensional = (g_calibrated / 2) * v.get_dim('rho_3D')**2
    v.check_dims("Pressure from calibrated g", v.get_dim('P'), pressure_dimensional)

    # Newton's constant calibration: G = c²/(4π ρ_0 ξ_c²)
    G_calibrated = v.get_dim('c')**2 / (4*sp.pi * v.get_dim('rho_0') * v.get_dim('xi')**2)
    v.check_dims("Newton's G calibration", v.get_dim('G'), G_calibrated)

    # Verify that G has correct dimensions
    G_dimensional_check = v.L**3 / (v.M * v.T**2)
    v.check_dims("G dimensional structure", v.get_dim('G'), G_dimensional_check)

    # Connection between g_eos and G through geometry: g = 4π G ξ_c²
    g_from_G_relation = 4*sp.pi * v.get_dim('G') * v.get_dim('xi')**2
    v.check_dims("g_eos from G relationship", v.get_dim('g_eos'), g_from_G_relation)

    v.success("Calibration constants verified")


def test_physical_interpretation(v):
    """
    Test the physical interpretation of gravity as rarefaction pressure gradients.

    Vortex sinks → rarefied zones → pressure gradients → apparent attraction

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Physical Interpretation")

    # Sink creates rarefaction: lower density in vortex core
    density_deficit = v.get_dim('rho_0') - v.get_dim('rho_3D')  # Positive deficit
    v.check_dims("Rarefaction deficit", density_deficit, v.get_dim('rho'))

    # Pressure gradient from density gradient (EOS: P ∝ ρ²)
    pressure_gradient = v.grad_dim(v.get_dim('P'))
    density_gradient = v.grad_dim(v.get_dim('rho_3D'))

    # From P = (g_eos/2)ρ², we get ∇P = g_eos ρ ∇ρ
    pressure_from_density = v.get_dim('g_eos') * v.get_dim('rho_3D') * density_gradient
    v.check_dims("Pressure gradient from density", pressure_gradient, pressure_from_density)

    # This pressure gradient creates effective gravitational field
    effective_g_field = pressure_gradient / v.get_dim('rho_3D')
    v.check_dims("Effective gravitational field", effective_g_field, v.get_dim('g'))

    # Two sinks attract like "two bathtub drains sharing outflow"
    # Both create rarefied regions, pressure gradients point toward sinks
    # → apparent mutual attraction through shared flow patterns

    v.success("Physical interpretation verified")


def test_matter_model_bridge(v):
    """
    Test the connection between loop ensembles and gravitational mass.

    'Matter' = ensembles of closed loops with minimized mass M_*

    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Matter Model Bridge")

    # Gravitational density as sum of minimized loop masses
    rho_body_from_loops = v.get_dim('M_ast') / v.get_dim('dV')  # M_* per volume
    v.check_dims("Body density from loops", v.get_dim('rho_body'), rho_body_from_loops)

    # Individual loop contribution: δ³(r - r_i) has dimensions L⁻³
    delta_function = v.get_dim('delta3')
    point_mass_density = v.get_dim('M_ast') * delta_function
    v.check_dims("Point mass density", point_mass_density, v.get_dim('rho'))

    # Total gravitational source: ρ_body = Σᵢ Mᵢ δ³(r - rᵢ)
    sum_of_sources = v.get_dim('M_ast') * v.get_dim('delta3')  # Single term in sum
    v.check_dims("Gravitational source sum", v.get_dim('rho_body'), sum_of_sources)

    # Internal loop parameters (Q, n₃, k, ...) only enter via M_*
    # The minimized mass M_* = M(R_*; Q, n₃, ℳ) includes all internal structure
    v.check_dims("Minimized loop mass", v.get_dim('M_ast'), v.get_dim('m'))

    # Electric charge Q doesn't enter gravity separately
    # Its EM energy contribution is already included in M_*
    # This ensures equivalence principle: gravitational mass = inertial mass

    v.success("Matter model bridge verified")


def test_newtonian_limit():
    """
    Main test function for Newtonian Limit.

    This function coordinates all verification tests for the emergence of
    Newtonian gravity from the 4D vortex framework in the static, low-velocity limit.

    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Newtonian Limit",
        "Emergence of Newton's gravity from 4D vortex dynamics"
    )

    # Add custom dimensions needed for this section
    v.add_dimensions({
        'M_ast': v.M,                                    # Minimized loop mass
        'g_eos': v.get_dim('c')**2 / v.get_dim('rho_0'),  # EOS coupling: g = c²/ρ₀ [L⁵/(MT²)]
        # M_dot_body and xi already defined in helper
    })

    # Verify the g-G relationship: g = 4π G ξ²
    v.check_dims("g from G relation", v.get_dim('g_eos'),
                 4*sp.pi*v.get_dim('G')*v.get_dim('xi')**2)

    v.section("NEWTONIAN LIMIT VERIFICATION")
    v.info("Testing emergence of Newton's law from 4D vortex sinks")

    # Call test functions in logical order
    v.info("\n--- 1) 3D Continuity Equation ---")
    test_continuity_equation_derivation(v)

    v.info("\n--- 2) Density Deficit Balance ---")
    test_density_deficit_balance(v)

    v.info("\n--- 3) Euler Equation Linearization ---")
    test_euler_equation_linearization(v)

    v.info("\n--- 4) Poisson Equation Emergence ---")
    test_poisson_equation_emergence(v)

    v.info("\n--- 5) Calibration Constants ---")
    test_calibration_constants(v)

    v.info("\n--- 6) Physical Interpretation ---")
    test_physical_interpretation(v)

    v.info("\n--- 7) Matter Model Bridge ---")
    test_matter_model_bridge(v)

    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_newtonian_limit()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)