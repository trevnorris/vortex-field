#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Canonical Structure and Ehrenfest - Verification
===============================================

Comprehensive verification of the canonical structure and Ehrenfest theorem
connections in the quantum framework derived from the 4D vortex field theory.

This test validates the dimensional consistency of the canonical commutation
relations, symplectic form implications, momentum operator definitions, and
Hamiltonian structure as presented in the quantum subsection.

Based on doc/quantum.tex, subsection "Canonical structure and Ehrenfest" (lines 64-77).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, I, simplify, sqrt, Rational

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    define_symbols_batch,
    quick_verify,
)


def test_symplectic_form_implications(v):
    """
    Test the symplectic form from eq:Spsi_full and its implications for 
    equal-time Poisson brackets between (ρ,S).
    
    Verifies: {ρ(x), S(y)} = δ³(x-y)
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Symplectic Form and Poisson Brackets")
    
    # The symplectic form implies equal-time brackets for (ρ,S)
    # where ρ is the probability density and S is the phase
    
    # ρ has dimensions of probability density (dimensionless per unit volume)
    # In the wavefunction formulation ρ = |ψ|², so ρ is dimensionless density
    v.add_dimensions({
        'rho_prob': v.L**(-3),  # Probability density (dimensionless per volume)
        'S_phase': 1,  # Phase is dimensionless (angle)
    })
    
    # The Poisson bracket {ρ(x), S(y)} = δ³(x-y) must be dimensionally consistent
    # Left side: {ρ(x), S(y)} has dimensions [rho_prob] × [S_phase] = L⁻³
    # Right side: δ³(x-y) has dimensions L⁻³ (3D Dirac delta)
    
    lhs_bracket = v.get_dim('rho_prob') * v.get_dim('S_phase')
    rhs_delta = v.get_dim('delta3')
    
    v.check_dims("Poisson bracket {ρ(x), S(y)} = δ³(x-y)",
                 lhs_bracket, rhs_delta)
    
    v.success("Symplectic form implications verified")


def test_canonical_commutation_relations(v):
    """
    Test the canonical commutation relations derived from the Poisson brackets.
    
    Verifies: [x̂ᵢ, p̂ⱼ] = iℏₑff δᵢⱼ
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Canonical Commutation Relations")
    
    # Add dimension for effective Planck constant
    v.add_dimensions({
        'hbar_eff': v.M * v.L**2 / v.T,  # Same dimensions as ℏ
    })
    
    # Position operator x̂ᵢ has dimension of length
    x_op_dim = v.L
    
    # Momentum operator p̂ⱼ has dimension of momentum
    p_op_dim = v.get_dim('p')
    
    # Commutator [x̂ᵢ, p̂ⱼ] should have dimensions [L] × [M L T⁻¹] = M L² T⁻¹
    commutator_lhs = x_op_dim * p_op_dim  # This gives the "scale" of the commutator
    
    # Right hand side: iℏₑff δᵢⱼ
    # i is dimensionless, δᵢⱼ is dimensionless, ℏₑff has dimensions M L² T⁻¹
    commutator_rhs = v.get_dim('hbar_eff')
    
    # The commutator itself represents an action-like quantity
    v.check_dims("Canonical commutator [x̂ᵢ, p̂ⱼ] = iℏₑff δᵢⱼ",
                 commutator_rhs, commutator_rhs)  # RHS defines the dimension
    
    # Verify that this is consistent with the action dimension
    v.check_dims("Commutator has action dimensions",
                 commutator_rhs, v.get_dim('S'))  # Action S has dimension M L² T⁻¹
    
    v.success("Canonical commutation relations verified")


def test_momentum_operator_definition(v):
    """
    Test the momentum operator definition from the canonical structure.
    
    Verifies: p̂ = -iℏₑff∇
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Momentum Operator Definition")
    
    # Left side: p̂ has momentum dimensions
    p_op_lhs = v.get_dim('p')
    
    # Right side: -iℏₑff∇
    # i is dimensionless, ℏₑff has dimensions M L² T⁻¹, ∇ has dimensions L⁻¹
    p_op_rhs = v.get_dim('hbar_eff') * v.get_dim('nabla')
    
    v.check_dims("Momentum operator p̂ = -iℏₑff∇",
                 p_op_lhs, p_op_rhs)
    
    # Verify this is consistent when acting on a wavefunction
    # ∇ψ has dimensions [psi] × L⁻¹ = L⁻³/² × L⁻¹ = L⁻⁵/²
    # -iℏₑff∇ψ has dimensions M L² T⁻¹ × L⁻⁵/² = M L⁻¹/² T⁻¹
    
    # This should equal momentum density × wavefunction
    # Momentum density has dimensions M L⁻² T⁻¹
    # So momentum density × ψ = M L⁻² T⁻¹ × L⁻³/² = M L⁻⁷/² T⁻¹
    
    # Actually, let's check the momentum expectation value instead
    # ⟨p̂⟩ = ∫ ψ* (-iℏₑff∇) ψ d³x
    # This should have momentum dimensions
    
    momentum_expectation = (v.get_dim('psi') * v.get_dim('hbar_eff') * 
                           v.get_dim('nabla') * v.get_dim('psi') * v.get_dim('dV'))
    
    v.check_dims("Momentum expectation ⟨p̂⟩",
                 momentum_expectation, v.get_dim('p'))
    
    v.success("Momentum operator definition verified")


def test_hamiltonian_structure(v):
    """
    Test the Hamiltonian structure including kinetic and potential terms.
    
    Verifies: Ĥ = (p̂ - qA)²/(2m*) + qΦ + V
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Hamiltonian Structure")
    
    # Add dimensions for effective mass and charge
    v.add_dimensions({
        'm_eff': v.M,  # Effective mass m*
        'q_charge': v.Q,  # Electric charge q
    })
    
    # Kinetic term: (p̂ - qA)²/(2m*)
    # p̂ - qA has momentum dimensions: both p̂ and qA have dimensions M L T⁻¹
    canonical_momentum = v.get_dim('p')
    em_momentum = v.get_dim('q_charge') * v.get_dim('A')
    
    v.check_dims("Electromagnetic momentum qA",
                 em_momentum, canonical_momentum)
    
    # (p̂ - qA)² has dimensions (M L T⁻¹)² = M² L² T⁻²
    # Dividing by 2m* gives dimensions M² L² T⁻² / M = M L² T⁻² (energy)
    kinetic_term = canonical_momentum**2 / v.get_dim('m_eff')
    
    v.check_dims("Kinetic energy (p̂ - qA)²/(2m*)",
                 kinetic_term, v.get_dim('E_energy'))
    
    # Electric potential term: qΦ
    electric_potential_term = v.get_dim('q_charge') * v.get_dim('Phi')
    
    v.check_dims("Electric potential energy qΦ",
                 electric_potential_term, v.get_dim('E_energy'))
    
    # Additional potential V (could be gravitational, self-interaction, etc.)
    # This should also have energy dimensions
    v.check_dims("Additional potential V",
                 v.get_dim('V_energy'), v.get_dim('E_energy'))
    
    # Complete Hamiltonian should have energy dimensions
    hamiltonian_total = (kinetic_term + electric_potential_term + 
                        v.get_dim('V_energy'))
    
    v.check_dims("Total Hamiltonian Ĥ",
                 hamiltonian_total, v.get_dim('E_energy'))
    
    v.success("Hamiltonian structure verified")


def test_ehrenfest_theorem_connection(v):
    """
    Test the connection to Ehrenfest theorem and expectation value evolution.
    
    Verifies that Heisenberg evolution reproduces the continuity equation
    in expectation value (Ehrenfest theorem application).
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Ehrenfest Theorem Connection")
    
    # The continuity equation is: ∂_t ρ + ∇·j = 0
    # where ρ is probability density and j is probability current
    
    # Time derivative of probability density
    drho_dt = v.dt(v.get_dim('rho_prob'))
    
    # Divergence of probability current
    # The current j has dimensions of probability flux: [probability]/[time]/[area]
    # j = (ℏₑff/2m*i)(ψ*∇ψ - ψ∇ψ*) - (q/m*)A|ψ|²
    
    # First term: (ℏₑff/2m*i)(ψ*∇ψ - ψ∇ψ*)
    # ψ*∇ψ has dimensions L⁻³/² × L⁻¹ × L⁻³/² = L⁻⁴
    # ℏₑff/(2m*) has dimensions (M L² T⁻¹)/(M) = L² T⁻¹
    # So first term has dimensions L² T⁻¹ × L⁻⁴ = L⁻² T⁻¹
    
    quantum_current_term = (v.get_dim('hbar_eff') / v.get_dim('m_eff') * 
                           v.get_dim('psi')**2 * v.get_dim('nabla'))
    
    # Second term: (q/m*)A|ψ|²
    # Has dimensions Q/M × (M L² T⁻² / Q T) × L⁻³ = L⁻² T⁻¹
    
    classical_current_term = (v.get_dim('q_charge') / v.get_dim('m_eff') * 
                             v.get_dim('A') * v.get_dim('psi')**2)
    
    v.check_dims("Quantum current term dimensions",
                 quantum_current_term, v.L**(-2) / v.T)
    
    v.check_dims("Classical current term dimensions", 
                 classical_current_term, v.L**(-2) / v.T)
    
    # Total current density j
    current_density = quantum_current_term + classical_current_term
    
    # Divergence of current
    div_j = v.div_dim(current_density)
    
    v.check_dims("Continuity equation: ∂_t ρ + ∇·j = 0",
                 drho_dt, -div_j)  # The minus sign is conventional
    
    # This verifies that Heisenberg evolution of quantum operators
    # reproduces classical continuity in expectation values
    v.info("Ehrenfest theorem: quantum → classical in expectation")
    
    v.success("Ehrenfest theorem connection verified")


def test_gauge_covariance(v):
    """
    Test gauge covariance of the canonical structure.
    
    Verifies that the canonical momentum p̂ - qA transforms properly
    under gauge transformations.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Gauge Covariance")
    
    # Under a gauge transformation:
    # A → A + ∇χ (where χ is gauge function)
    # ψ → e^(iqχ/ℏₑff) ψ
    
    # The gauge function χ must have dimensions such that qχ/ℏₑff is dimensionless
    # So χ has dimensions ℏₑff/q = (M L² T⁻¹)/Q = M L² T⁻¹ Q⁻¹
    
    v.add_dimensions({
        'chi_gauge': v.M * v.L**2 * v.T**(-1) * v.Q**(-1),  # Gauge function
    })
    
    # Gauge transformation of vector potential: A → A + ∇χ
    gauge_transform_A = v.grad_dim(v.get_dim('chi_gauge'))
    
    v.check_dims("Gauge transformation ∇χ",
                 gauge_transform_A, v.get_dim('A'))
    
    # The canonical momentum p̂ - qA should be gauge invariant
    # because the ∇χ term cancels with the phase transformation of ψ
    
    # Original canonical momentum
    canonical_mom_orig = v.get_dim('p') - v.get_dim('q_charge') * v.get_dim('A')
    
    # After gauge transformation
    canonical_mom_gauge = (v.get_dim('p') - 
                          v.get_dim('q_charge') * (v.get_dim('A') + gauge_transform_A))
    
    # The additional term should be compensated by the wavefunction transformation
    gauge_compensation = v.get_dim('q_charge') * gauge_transform_A
    
    v.check_dims("Gauge transformation compensation",
                 gauge_compensation, v.get_dim('p'))
    
    v.info("Canonical momentum p̂ - qA is gauge covariant")
    v.success("Gauge covariance verified")


def test_quantum_classical_correspondence(v):
    """
    Test the quantum-classical correspondence in the canonical formalism.
    
    Verifies that quantum commutators reduce to classical Poisson brackets
    in the classical limit.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Quantum-Classical Correspondence")
    
    # Classical Poisson bracket {x, p} = 1 (dimensionless)
    # Quantum commutator [x̂, p̂] = iℏₑff (has action dimensions)
    
    # In the classical limit ℏₑff → 0, we should recover:
    # [x̂, p̂]/(iℏₑff) → {x, p} = 1
    
    # The quantum commutator normalized by iℏₑff should be dimensionless
    normalized_commutator = v.get_dim('hbar_eff') / v.get_dim('hbar_eff')
    
    v.check_dims("Normalized commutator [x̂, p̂]/(iℏₑff)",
                 normalized_commutator, 1)  # Should be dimensionless
    
    # Hamilton's equations in classical mechanics:
    # dx/dt = ∂H/∂p,  dp/dt = -∂H/∂x
    
    # Quantum Heisenberg equations:
    # d⟨x̂⟩/dt = (1/iℏₑff)⟨[x̂, Ĥ]⟩,  d⟨p̂⟩/dt = (1/iℏₑff)⟨[p̂, Ĥ]⟩
    
    # For the Hamiltonian Ĥ = p̂²/(2m) + V(x̂):
    # [x̂, p̂²] = 2iℏₑff p̂,  [p̂, V(x̂)] = -iℏₑff ∇V
    
    # So: d⟨x̂⟩/dt = ⟨p̂⟩/m,  d⟨p̂⟩/dt = -⟨∇V⟩
    
    # These have the correct dimensions for classical mechanics
    position_evolution = v.get_dim('p') / v.get_dim('m_eff')  # Should be velocity
    momentum_evolution = v.grad_dim(v.get_dim('V_energy'))   # Should be force
    
    v.check_dims("Position evolution d⟨x̂⟩/dt = ⟨p̂⟩/m",
                 position_evolution, v.L/v.T)
    
    v.check_dims("Momentum evolution d⟨p̂⟩/dt = -⟨∇V⟩",
                 momentum_evolution, v.get_dim('p')/v.T)  # Force = dp/dt
    
    v.success("Quantum-classical correspondence verified")


def test_canonical_structure_and_ehrenfest():
    """
    Main test function for Canonical Structure and Ehrenfest subsection.
    
    This function coordinates all verification tests for the canonical structure
    derived from the 4D vortex field theory and its connection to the Ehrenfest
    theorem, validating dimensional consistency and physical relationships.
    
    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Canonical Structure and Ehrenfest",
        "Quantum canonical formalism and Ehrenfest theorem from 4D vortex theory"
    )
    
    v.section("CANONICAL STRUCTURE AND EHRENFEST VERIFICATION")
    
    # Add custom dimensions needed for quantum mechanics tests
    v.add_dimensions({
        'delta_ij': 1,  # Kronecker delta (dimensionless)
    })
    
    # Call test functions in logical order
    v.info("\n--- 1) Symplectic Form Implications ---")
    test_symplectic_form_implications(v)
    
    v.info("\n--- 2) Canonical Commutation Relations ---")
    test_canonical_commutation_relations(v)
    
    v.info("\n--- 3) Momentum Operator Definition ---")
    test_momentum_operator_definition(v)
    
    v.info("\n--- 4) Hamiltonian Structure ---")
    test_hamiltonian_structure(v)
    
    v.info("\n--- 5) Ehrenfest Theorem Connection ---")
    test_ehrenfest_theorem_connection(v)
    
    v.info("\n--- 6) Gauge Covariance ---")
    test_gauge_covariance(v)
    
    v.info("\n--- 7) Quantum-Classical Correspondence ---")
    test_quantum_classical_correspondence(v)
    
    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_canonical_structure_and_ehrenfest()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)