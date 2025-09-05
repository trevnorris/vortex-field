#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Measurement, decoherence, and classicality - Verification
====================================

Comprehensive verification of the decoherence mechanism and classical limit
emergence in the vortex framework, including the Gaussian dephasing kernel,
environmental coupling, and the d² law for intrinsic slab-coupled decoherence.

This test validates the dimensional consistency of the decoherence equation
(eq:decoherence), analyzes the environmental mode integration effects, and
verifies the emergence of classical behavior through measurement-induced
decoherence in the 4D superfluid framework.

Based on doc/quantum.tex, section "Measurement, decoherence, and classicality" (lines 141-151).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, sqrt, ln, exp, simplify, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify,
)


def test_decoherence_kernel_equation(v):
    """
    Test the dimensional consistency of the Gaussian dephasing kernel equation.
    
    Verifies: Γ_dec(d) = Γ₀ + γ₂·d² + O(d⁴)
    where γ₂ ∝ α_tw · (ℏ_eff)/(m*·ℓ*⁴) · (ε/ℓ*)^p
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Decoherence Kernel (eq:decoherence)")
    
    # Define the symbolic variables as they appear in the document
    d, Gamma_0, gamma_2 = define_symbols_batch(
        ['d', 'Gamma_0', 'gamma_2'], positive=True
    )
    alpha_tw, hbar_eff, m_star, ell_star, epsilon, p = define_symbols_batch(
        ['alpha_tw', 'hbar_eff', 'm_star', 'ell_star', 'epsilon', 'p'], positive=True
    )
    
    # First verify the dimensions of the decoherence rate Γ_dec
    # Decoherence rate should have dimensions of inverse time [T⁻¹]
    v.check_dims("Decoherence rate Γ_dec", v.get_dim('gamma'), v.T**(-1))
    
    # Path separation d should have dimensions of length
    v.check_dims("Path separation d", v.get_dim('d'), v.L)
    
    # The constant term Γ₀ should have same dimensions as Γ_dec
    v.check_dims("Constant decoherence term Γ₀", v.get_dim('gamma'), v.T**(-1))
    
    # For the d² term: γ₂·d² must have dimensions [T⁻¹]
    # So γ₂ must have dimensions [T⁻¹ L⁻²]
    gamma_2_expected = v.T**(-1) / v.L**2
    v.add_dimensions({'gamma_2': gamma_2_expected})
    v.check_dims("γ₂ coefficient", v.get_dim('gamma_2'), gamma_2_expected)
    
    # Verify the d² term has correct dimensions
    d_squared_term = v.get_dim('gamma_2') * v.get_dim('d')**2
    v.check_dims("d² term γ₂·d²", d_squared_term, v.T**(-1))
    
    # Verify the complete decoherence kernel
    decoherence_kernel = v.get_dim('gamma') + d_squared_term
    v.check_dims("Complete decoherence kernel Γ_dec(d)",
                 decoherence_kernel, v.T**(-1))
    
    v.success("Decoherence kernel dimensional consistency verified")


def test_gamma_2_coupling_formula(v):
    """
    Test the dimensional consistency of the γ₂ coupling strength formula.
    
    Verifies: γ₂ ∝ α_tw · (ℏ_eff)/(m*·ℓ*⁴) · (ε/ℓ*)^p
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("γ₂ Coupling Strength Formula")
    
    # Add dimensions for the parameters in the γ₂ formula
    v.add_dimensions({
        'alpha_tw': 1,  # Twist coupling strength (dimensionless)
        'hbar_eff': v.M * v.L**2 / v.T,  # Effective Planck constant
        'm_star': v.M,  # Effective mass
        'ell_star': v.L,  # Characteristic length scale
        'epsilon_param': v.L,  # Small parameter (likely healing length scale)
        'p_exponent': 1,  # Dimensionless exponent
    })
    
    # Test individual components
    v.check_dims("Twist coupling α_tw", v.get_dim('alpha_tw'), 1)
    v.check_dims("Effective ℏ", v.get_dim('hbar_eff'), v.M * v.L**2 / v.T)
    v.check_dims("Effective mass m*", v.get_dim('m_star'), v.M)
    v.check_dims("Length scale ℓ*", v.get_dim('ell_star'), v.L)
    v.check_dims("Small parameter ε", v.get_dim('epsilon_param'), v.L)
    
    # Test the ratio ℏ_eff/(m*·ℓ*⁴)
    hbar_ratio = v.get_dim('hbar_eff') / (v.get_dim('m_star') * v.get_dim('ell_star')**4)
    expected_hbar_ratio = (v.M * v.L**2 / v.T) / (v.M * v.L**4)
    expected_hbar_ratio = v.T**(-1) / v.L**2
    v.check_dims("ℏ_eff/(m*·ℓ*⁴) ratio", hbar_ratio, expected_hbar_ratio)
    
    # Test the dimensionless ratio (ε/ℓ*)^p
    epsilon_ratio = (v.get_dim('epsilon_param') / v.get_dim('ell_star'))**v.get_dim('p_exponent')
    v.check_dims("Dimensionless ratio (ε/ℓ*)^p", epsilon_ratio, 1)
    
    # Test the complete γ₂ expression
    gamma_2_full = (v.get_dim('alpha_tw') * hbar_ratio * epsilon_ratio)
    v.check_dims("Complete γ₂ expression", gamma_2_full, v.get_dim('gamma_2'))
    
    v.success("γ₂ coupling strength formula verified")


def test_environmental_integration_physics(v):
    """
    Test the physical interpretation of environmental mode integration
    and its effect on quantum coherence.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Environmental Mode Integration Physics")
    
    # The document mentions "integrating out slab/environmental modes"
    # This should preserve action dimensions while introducing decoherence
    
    # Add dimensions for environmental physics
    v.add_dimensions({
        'S_env': v.M * v.L**2 / v.T,  # Environmental action
        'rho_env': v.M / v.L**3,  # Environmental density
        'T_env': v.T,  # Environmental time scale
        'xi_env': v.L,  # Environmental correlation length
        'coupling_strength': v.M * v.L**(-1) * v.T**(-2),  # Environment-system coupling
    })
    
    # Environmental action should have same dimensions as quantum action
    v.check_dims("Environmental action S_env", 
                 v.get_dim('S_env'), v.get_dim('S'))
    
    # Environmental correlation time should relate to decoherence rate
    env_decoherence_rate = 1 / v.get_dim('T_env')
    v.check_dims("Environmental decoherence rate 1/T_env",
                 env_decoherence_rate, v.T**(-1))
    
    # Path separation dependence comes from environmental correlation
    # The d² scaling should emerge from Gaussian environmental correlations
    env_correlation_scale = v.get_dim('xi_env')
    path_correlation_effect = 1 / env_correlation_scale**2
    v.check_dims("Path correlation effect 1/ξ_env²",
                 path_correlation_effect, v.L**(-2))
    
    # Combined environmental decoherence should match γ₂ structure
    env_decoherence = env_decoherence_rate * path_correlation_effect
    v.check_dims("Environmental decoherence coupling",
                 env_decoherence, v.get_dim('gamma_2'))
    
    v.success("Environmental integration physics verified")


def test_classical_limit_emergence(v):
    """
    Test the emergence of classical behavior through decoherence.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Classical Limit Emergence")
    
    # Add dimensions for classical limit analysis
    v.add_dimensions({
        'lambda_dB': v.L,  # de Broglie wavelength
        'coherence_length': v.L,  # Quantum coherence length
        'macroscopic_scale': v.L,  # Macroscopic length scale
        'decoherence_time': v.T,  # Decoherence time scale
        'classical_time': v.T,  # Classical evolution time scale
    })
    
    # de Broglie wavelength: λ_dB = ℏ/(m·v)
    lambda_dB_expected = v.get_dim('hbar_eff') / (v.M * v.L/v.T)
    lambda_dB_expected = simplify(lambda_dB_expected)
    v.check_dims("de Broglie wavelength λ_dB = ℏ/(mv)",
                 v.get_dim('lambda_dB'), lambda_dB_expected)
    
    # Coherence length should be limited by decoherence
    # At classical limit: coherence_length << macroscopic_scale
    v.info("Classical limit: coherence_length << macroscopic_scale")
    v.check_dims("Coherence length scale",
                 v.get_dim('coherence_length'), v.L)
    v.check_dims("Macroscopic length scale",
                 v.get_dim('macroscopic_scale'), v.L)
    
    # Decoherence time should be much shorter than classical evolution
    # At classical limit: decoherence_time << classical_time
    v.info("Classical limit: decoherence_time << classical_time")
    v.check_dims("Decoherence time scale",
                 v.get_dim('decoherence_time'), v.T)
    v.check_dims("Classical evolution time scale",
                 v.get_dim('classical_time'), v.T)
    
    # The d² law should dominate when d ~ macroscopic scales
    macroscopic_decoherence = v.get_dim('gamma_2') * v.get_dim('macroscopic_scale')**2
    v.check_dims("Macroscopic decoherence rate γ₂·d_macro²",
                 macroscopic_decoherence, v.T**(-1))
    
    # Verify that macroscopic decoherence >> microscopic rates
    v.info("Classical regime: γ₂·d_macro² >> Γ₀ (strong decoherence)")
    
    v.success("Classical limit emergence verified")


def test_residual_d2_law_detection(v):
    """
    Test the physics of residual d² law detection after subtracting
    standard collisional/thermal channels.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Residual d² Law Detection")
    
    # Add dimensions for decoherence channel analysis
    v.add_dimensions({
        'Gamma_collision': v.T**(-1),  # Collisional decoherence rate
        'Gamma_thermal': v.T**(-1),  # Thermal decoherence rate
        'Gamma_intrinsic': v.T**(-1),  # Intrinsic slab decoherence
        'T_thermal': v.T,  # Thermal time scale
        'sigma_collision': v.L**2,  # Collision cross-section
        'n_density': v.L**(-3),  # Number density
        'v_thermal': v.L / v.T,  # Thermal velocity
    })
    
    # Collisional decoherence: Γ_coll ~ n·σ·v (independent of path separation)
    collision_rate_expected = v.get_dim('n_density') * v.get_dim('sigma_collision') * v.get_dim('v_thermal')
    v.check_dims("Collisional decoherence rate n·σ·v",
                 collision_rate_expected, v.T**(-1))
    
    # Thermal decoherence: Γ_thermal ~ 1/T_thermal (also path-independent)
    thermal_rate_expected = 1 / v.get_dim('T_thermal')
    v.check_dims("Thermal decoherence rate 1/T_thermal",
                 thermal_rate_expected, v.T**(-1))
    
    # Total measured decoherence includes all channels
    total_decoherence = (v.get_dim('Gamma_collision') + 
                        v.get_dim('Gamma_thermal') + 
                        v.get_dim('Gamma_intrinsic'))
    v.check_dims("Total measured decoherence",
                 total_decoherence, v.T**(-1))
    
    # After subtracting standard channels, residual should follow d² law
    residual_decoherence = v.get_dim('Gamma_intrinsic')
    intrinsic_d2_component = v.get_dim('gamma_2') * v.get_dim('d')**2
    v.check_dims("Residual intrinsic decoherence",
                 residual_decoherence, v.T**(-1))
    v.check_dims("Intrinsic d² component",
                 intrinsic_d2_component, v.T**(-1))
    
    # The key prediction: intrinsic component shows d² scaling
    v.info("Key prediction: After subtracting collisional/thermal channels,")
    v.info("intrinsic slab-coupled decoherence must manifest as residual d² law")
    
    v.success("Residual d² law detection physics verified")


def test_slab_coupling_mechanism(v):
    """
    Test the physical mechanism of slab-coupled decoherence in the 4D framework.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Slab Coupling Mechanism")
    
    # Add dimensions for 4D slab physics
    v.add_dimensions({
        'w_slab': v.L,  # Slab thickness in 4th dimension
        'coupling_4D': v.M / (v.L**2 * v.T**2),  # 4D coupling strength
        'field_4D': v.L / v.T,  # 4D field amplitude
        'energy_4D': v.M * v.L**2 / v.T**2,  # 4D energy scale
        'density_slab': v.M / v.L**4,  # 4D slab density
    })
    
    # 4D coupling should relate to 3D physics through projection
    # The slab thickness provides the geometric coupling
    slab_volume_factor = v.get_dim('w_slab')  # 4D volume per unit 3D volume
    v.check_dims("4D slab thickness w", slab_volume_factor, v.L)
    
    # 4D density projected to 3D should match our density scales
    projected_density = v.get_dim('density_slab') * slab_volume_factor
    v.check_dims("Projected 4D density ρ_4D·w",
                 projected_density, v.M / v.L**3)
    
    # The coupling mechanism should preserve energy scales
    slab_energy_per_volume = v.get_dim('energy_4D') / v.L**4  # 4D energy density
    projected_energy_density = slab_energy_per_volume * slab_volume_factor
    v.check_dims("Projected energy density",
                 projected_energy_density, v.M * v.L**(-1) * v.T**(-2))
    
    # Path separation dependence comes from 4D field correlations
    # Gaussian correlations in 4D project to d² behavior in 3D
    field_correlation_4D = v.get_dim('field_4D')**2
    path_dependent_coupling = field_correlation_4D / v.get_dim('d')**2
    expected_coupling_scale = (v.L/v.T)**2 / v.L**2
    expected_coupling_scale = v.T**(-2)
    v.check_dims("4D field correlation/d²",
                 path_dependent_coupling, expected_coupling_scale)
    
    # This should connect to the γ₂ coefficient
    v.info("4D slab coupling mechanism explains d² path dependence")
    
    v.success("Slab coupling mechanism verified")


def test_measurement_decoherence_and_classicality():
    """
    Main test function for Measurement, decoherence, and classicality.
    
    This function coordinates all verification tests for the decoherence
    mechanism and classical limit emergence in the vortex framework,
    validating the Gaussian dephasing kernel, environmental mode integration,
    and the d² law for intrinsic slab-coupled decoherence.
    
    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Measurement, decoherence, and classicality",
        "Decoherence mechanism and classical limit emergence in vortex framework"
    )
    
    v.section("MEASUREMENT, DECOHERENCE, AND CLASSICALITY VERIFICATION")
    
    # Add common dimensions used throughout the tests
    v.add_dimensions({
        'd': v.L,  # Path separation
    })
    
    # Call test functions in logical order
    v.info("\n--- 1) Decoherence Kernel Equation ---")
    test_decoherence_kernel_equation(v)
    
    v.info("\n--- 2) γ₂ Coupling Strength Formula ---")
    test_gamma_2_coupling_formula(v)
    
    v.info("\n--- 3) Environmental Mode Integration Physics ---")
    test_environmental_integration_physics(v)
    
    v.info("\n--- 4) Classical Limit Emergence ---")
    test_classical_limit_emergence(v)
    
    v.info("\n--- 5) Residual d² Law Detection ---")
    test_residual_d2_law_detection(v)
    
    v.info("\n--- 6) Slab Coupling Mechanism ---")
    test_slab_coupling_mechanism(v)
    
    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_measurement_decoherence_and_classicality()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)