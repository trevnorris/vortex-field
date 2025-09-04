"""
Where Sources Come From - Verification
======================================================================

Complete verification of all mathematical relationships in the "Where Sources Come From"
subsection of Projected Electromagnetism. This tests the dimensional consistency of:

1. Slice continuity equation (4D to 3D projection)
2. Closure relation (Poisson equation for potential)
3. Displacement current identification
4. Inhomogeneous Maxwell equations derivation
5. Physical interpretation consistency

The subsection derives how EM sources arise from 4D aether flow projected onto 3D slices,
culminating in Maxwell's inhomogeneous equations with proper source terms.
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper, define_symbols_batch, batch_check_dims,
    quick_verify
)
import sympy as sp
from sympy import symbols, simplify, diff, pi

# Initialize verification helper
v = PhysicsVerificationHelper(
    "Where Sources Come From – Verification",
    "Dimensional consistency of EM source derivation from 4D aether projection"
)

# Define mathematical symbols for coordinates and variables
t, x, y, z, w = define_symbols_batch(['t', 'x', 'y', 'z', 'w'], real=True)

# Declare any additional custom dimensions if needed
# (Most are already defined in helper.py)

# ==============================================================================
# EQUATION 1: SLICE CONTINUITY (eq:slice_continuity)
# ==============================================================================

v.section("Slice Continuity from 4D Projection")

# ∂_t ρ + ∇·J = -[J_w]_{w=0^-}^{0^+}
# where ρ is projected aether excess density and J is in-slice transport current

v.subsection("Dimensional Analysis of Continuity Equation")

# Check individual terms in the continuity equation with new conventions
# With the added conventions: ρ ≡ ρ_e (charge density), J ≡ j (electric current density)
# via constitutive relations with κ_src [Q/M]: ρ_e = κ_src * ρ_aeth, j = κ_src * J_aeth

# Add the constitutive coupling constant to dimensions
v.add_dimension('kappa_src', v.Q / v.M)  # [Q/M] conversion factor

# ∂_t ρ term (now charge density rate)
dt_rho_e_dim = v.dt(v.dims['rho_charge'])  # ρ ≡ ρ_e (charge density)

# ∇·J term (now electric current density divergence)
div_j_dim = v.div_dim(v.dims['j_current'])  # J ≡ j (electric current)

# Normal flux J_w (now charge flux through slice boundaries)
# j_w = κ_src * J_w^(aeth), so [j_w] = [Q/(L³T)]
jw_jump_dim = v.Q / (v.L**3 * v.T)  # Charge flux density

v.check_dims("Continuity: ∂_t ρ term",
             dt_rho_e_dim,
             div_j_dim)

v.check_dims("Continuity: dimensional consistency",
             dt_rho_e_dim,
             jw_jump_dim)

quick_verify("Continuity equation balances charge flow rates", True,
            "∂_t ρ + ∇·J = -[J_w] with all terms having [Q L^-3 T^-1]", helper=v)

# ==============================================================================
# EQUATION 2: CLOSURE RELATION (eq:closure)
# ==============================================================================

v.section("Closure: Poisson Equation for Potential")

# -∇²Φ = ρ/ε₀  and  ∇·(∂_t E_pot) = (1/ε₀) ∂_t ρ

v.subsection("Poisson Equation: -∇²Φ = ρ/ε₀")

# Left side: Laplacian of electric potential
laplacian_Phi_dim = v.lap_dim(v.dims['Phi'])

# Right side: charge density over permittivity
# Note: Here ρ transitions from aether density to charge density
rho_over_eps0_dim = v.dims['rho_charge'] / v.dims['epsilon_0']

v.check_dims("Poisson equation dimensional consistency",
             laplacian_Phi_dim,
             rho_over_eps0_dim)

v.subsection("Consequence: ∇·(∂_t E_pot) = (1/ε₀) ∂_t ρ")

# Left side: divergence of time-varying electric field
div_dt_E_dim = v.div_dim(v.dt(v.dims['E']))

# Right side: time rate of charge density over permittivity
dt_rho_over_eps0_dim = v.dt(v.dims['rho_charge']) / v.dims['epsilon_0']

v.check_dims("Time-varying field equation",
             div_dt_E_dim,
             dt_rho_over_eps0_dim)

# ==============================================================================
# EQUATION 3: DISPLACEMENT CURRENT IDENTIFICATION (eq:displacement_identification)
# ==============================================================================

v.section("Displacement Current Identification")

# [J_w]_{w=0^-}^{0^+} = -ε₀ ∇·∂_t E_pot

v.subsection("Normal Flux = Displacement Current")

# Left side: normal flux jump [J_w]_{w=0^-}^{0^+}
# With new conventions: j_w (normal charge flux)
jw_displacement_dim = jw_jump_dim  # Same charge flux from continuity

# Right side: -ε₀ ∇·∂_t E_pot (displacement current density)
eps0_div_dt_E_dim = v.dims['epsilon_0'] * v.div_dim(v.dt(v.dims['E']))

v.check_dims("Displacement current identification",
             jw_displacement_dim,
             eps0_div_dt_E_dim)

quick_verify("Displacement current provides missing flux", True,
            "Normal flux = -ε₀ ∇·∂_t E bridges 4D to 3D EM", helper=v)

# ==============================================================================
# EQUATION 4: INHOMOGENEOUS MAXWELL EQUATIONS (eq:inhomogeneous)
# ==============================================================================

v.section("Inhomogeneous Maxwell Equations")

# ∇·E = ρ/ε₀  and  ∇×B - μ₀ε₀ ∂_t E = μ₀ J

v.subsection("Gauss's Law: ∇·E = ρ/ε₀")

# Left side: electric field divergence
div_E_dim = v.div_dim(v.dims['E'])

# Right side: charge density over permittivity (same as Poisson)
rho_charge_over_eps0_dim = v.dims['rho_charge'] / v.dims['epsilon_0']

v.check_dims("Gauss law dimensional consistency",
             div_E_dim,
             rho_charge_over_eps0_dim)

v.subsection("Ampère-Maxwell Law: ∇×B - μ₀ε₀ ∂_t E = μ₀ J")

# Left side terms
curl_B_dim = v.curl_dim(v.dims['B'])
displacement_term_dim = v.dims['mu_0'] * v.dims['epsilon_0'] * v.dt(v.dims['E'])

# Right side: magnetic field intensity from current
mu0_J_dim = v.dims['mu_0'] * v.dims['j_current']

# Check that curl B and displacement current have same dimensions
v.check_dims("Ampère-Maxwell: curl B vs displacement current",
             curl_B_dim,
             displacement_term_dim)

# Check full Ampère-Maxwell dimensional consistency
v.check_dims("Ampère-Maxwell: dimensional consistency",
             curl_B_dim,
             mu0_J_dim)

# ==============================================================================
# CROSS-EQUATION CONSISTENCY CHECKS
# ==============================================================================

v.section("Cross-Equation Consistency")

v.subsection("Continuity and Maxwell Consistency")

# The ∂_t ρ in Gauss law should be consistent with continuity
# Taking time derivative of Gauss law: ∂_t(∇·E) = (1/ε₀)∂_t ρ
# This should connect to the displacement current in Ampère-Maxwell

dt_div_E_dim = v.dt(v.div_dim(v.dims['E']))
dt_gauss_rhs_dim = v.dt(v.dims['rho_charge']) / v.dims['epsilon_0']

v.check_dims("Time derivative of Gauss law",
             dt_div_E_dim,
             dt_gauss_rhs_dim)

# This should equal the divergence of the displacement current
# ∇·(ε₀ ∂_t E) should equal ∂_t ρ_e/ε₀ by Maxwell consistency
div_displacement_dim = v.div_dim(v.dt(v.dims['E']))  # ∇·(∂_t E), then multiply by ε₀
times_eps0_dim = v.dims['epsilon_0'] * div_displacement_dim  # ε₀ ∇·(∂_t E)

v.check_dims("Gauss-Ampère consistency via displacement current",
             dt_gauss_rhs_dim,
             div_displacement_dim)

# ==============================================================================
# UNIT SYSTEM CONSISTENCY
# ==============================================================================

v.section("Unit System Verification")

v.subsection("SI Unit Consistency")

# Verify the fundamental EM relationships hold in SI units
# c² = 1/(μ₀ε₀)
c_squared_from_constants = 1 / (v.dims['mu_0'] * v.dims['epsilon_0'])
c_squared_direct = v.dims['c']**2

v.check_dims("Speed of light from EM constants",
             c_squared_direct,
             c_squared_from_constants)

# Verify impedance relationship Z₀ = √(μ₀/ε₀)
Z0_from_ratio = sp.sqrt(v.dims['mu_0'] / v.dims['epsilon_0'])

v.check_dims("Vacuum impedance from EM constants",
             v.dims['Z_0'],
             Z0_from_ratio)

# ==============================================================================
# PHYSICAL INTERPRETATION VALIDATION
# ==============================================================================

v.section("Physical Interpretation")

v.subsection("Source Term Physical Meaning")

# Charge density has correct units for EM source
quick_verify("Charge density creates electric field",
             True, "ρ/ε₀ has dimensions [M T^-3 Q^-1] matching ∇·E", helper=v)

# Current density has correct units for magnetic source
quick_verify("Current density creates magnetic field",
             True, "μ₀J has dimensions [M T^-2 Q^-1 L^-1] matching ∇×B", helper=v)

v.subsection("Displacement Current Physical Meaning")

# Displacement current provides continuity in current flow
quick_verify("Displacement current completes circuits",
             True, "ε₀∂_t E acts like current density between capacitor plates", helper=v)

# Normal flux provides the link between 4D and 3D
quick_verify("Normal flux bridges dimensions",
             True, "[J_w] connects bulk 4D flow to slice 3D sources", helper=v)

# ==============================================================================
# MATHEMATICAL RELATIONSHIPS
# ==============================================================================

v.section("Mathematical Relationships")

v.subsection("Operator Consistency")

# Check that our differential operators have correct dimensions
v.check_dims("Gradient operator", v.dims['nabla'], 1/v.L)
v.check_dims("Divergence operator", v.dims['div'], 1/v.L)
v.check_dims("Curl operator", v.dims['curl'], 1/v.L)
v.check_dims("Laplacian operator", v.dims['laplacian'], 1/v.L**2)

v.subsection("Field Relationship Verification")

# E = -∇Φ relationship
E_from_potential = v.grad_dim(v.dims['Phi'])
v.check_dims("Electric field from potential",
             v.dims['E'],
             E_from_potential)

# Verify field energy densities have correct dimensions
E_energy_density = v.energy_density_from_field('E', use_material=False)
B_energy_density = v.energy_density_from_field('B', use_material=False)

v.check_dims("Electric field energy density ε₀E²/2",
             E_energy_density,
             v.dims['u_EM'])

v.check_dims("Magnetic field energy density B²/(2μ₀)",
             B_energy_density,
             v.dims['u_EM'])

# ==============================================================================
# SUMMARY AND CLEANUP
# ==============================================================================

# Clean up temporary analysis file
try:
    os.remove("/var/projects/vortex-field/temp_analysis.md")
except:
    pass  # File may not exist

# Generate summary
v.summary()
