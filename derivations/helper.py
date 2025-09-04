"""
Physics Verification Test Framework
==============================================================
Unified helper functions and classes for SymPy-based physics verifications.
Simplifies writing verification tests for complex physics relationships.

Supports multiple unit systems (SI, Gaussian, Heaviside-Lorentz) with
consistent dimensional analysis and smart symbol lookup with suggestions.
"""

import sympy as sp
from sympy import (symbols, Function, diff, simplify, solve, Eq, pi, sqrt,
                   limit, oo, exp, log, integrate, Matrix, sinh, cosh, tanh,
                   sech, atan, sin, cos, Rational, ln, Abs, I, re, im, sign,
                   factorial, Sum, Product)
from sympy.vector import CoordSys3D, gradient, divergence, curl
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
from enum import Enum
import difflib
import warnings
import sys
import argparse


class UnitSystem(Enum):
    """Supported unit systems for electromagnetic and gravitational equations."""
    SI = "SI"
    GAUSSIAN = "Gaussian"
    HEAVISIDE_LORENTZ = "Heaviside-Lorentz"
    # NATURAL = "Natural"  # Reserved for future: c = ħ = 1


class UnknownSymbolError(KeyError):
    """Custom exception for unknown symbols with suggestions."""
    pass


class DimensionalMismatchError(ValueError):
    """Exception raised when dimensions don't match in an expression."""
    pass


def _parse_global_quiet_flag():
    """
    Parse command line for --quiet flag without interfering with other arguments.
    Returns True if --quiet is present, False otherwise.
    """
    # Check if --quiet is in sys.argv without using full argument parsing
    # This allows scripts to have their own argument parsing if needed
    return '--quiet' in sys.argv or '-q' in sys.argv


class SymbolRegistry:
    """Registry for managing physics symbols with aliases and fuzzy matching."""

    # Canonical aliases mapping
    ALIASES = {
        # Field aliases
        "Efield": "E_field", "Bfield": "B_field",
        "Hfield": "H_field", "Dfield": "D_field",
        "Phi_el": "Phi_E", "Phi_em": "Phi_E",
        "Avec": "A", "A_vec": "A",

        # Density aliases
        "rho_e": "rho_charge", "rho_el": "rho_charge",
        "j_m": "j_mass", "Jm": "j_mass", "J_mass": "j_mass",
        "j": "j_current",  # Electric current density (canonical)
        "Jem": "j_current",  # Clear alias for EM current
        "rho4D": "rho_4", "rho3D": "rho",

        # Potential aliases
        "phig": "Phi_g", "Ag": "A_g", "Ag_vec": "A_g",

        # Common alternatives
        "vg": "v",
        "grav": "g",
        "omega": "omega_freq",
        "nu": "f", "frequency": "f",

        # Greek letter variants
        "lambda": "wavelength", "Lambda": "wavelength",
        "mu": "mu_material", "mu0": "mu_0",
        "epsilon": "epsilon_material", "eps0": "epsilon_0",
        "sigma_w": "sigma_kernel", "Sigma": "sigma_surface",

        # Additional common aliases
        "rho_m": "rho", "rho_q": "rho_charge",
        "J_ang": "J_angular", "Lagrangian": "mathcal_L",
        "sigma_e": "sigma_conductivity",
    }

    @classmethod
    def canonicalize(cls, name: str) -> str:
        """Convert alias to canonical name."""
        canonical = cls.ALIASES.get(name, name)
        return canonical

    @classmethod
    def get_symbol(cls, name: str, registry: dict, suggest: bool = True,
                   n_suggestions: int = 5) -> Any:
        """
        Look up symbol with fuzzy matching suggestions.

        Args:
            name: Symbol name to look up
            registry: Dictionary of available symbols
            suggest: Whether to provide suggestions for misses
            n_suggestions: Number of suggestions to provide

        Returns:
            Symbol dimension if found

        Raises:
            UnknownSymbolError: If symbol not found (with suggestions)
        """
        # Try canonical name
        canonical = cls.canonicalize(name)
        if canonical in registry:
            return registry[canonical]

        # Try exact match
        if name in registry:
            return registry[name]

        # Symbol not found
        if suggest:
            # Build pool of candidates
            pool = list(registry.keys()) + list(cls.ALIASES.keys())
            matches = difflib.get_close_matches(name, pool, n=n_suggestions, cutoff=0.6)

            if matches:
                suggestions = []
                for match in matches:
                    canonical_match = cls.canonicalize(match)
                    if canonical_match in registry:
                        suggestions.append(canonical_match)
                    elif match in registry:
                        suggestions.append(match)

                if suggestions:
                    unique_suggestions = list(set(suggestions))[:n_suggestions]
                    hint = f"\n  Did you mean: {', '.join(unique_suggestions)}?"
                else:
                    hint = ""
            else:
                hint = ""

            raise UnknownSymbolError(f"Unknown symbol '{name}'.{hint}")

        raise UnknownSymbolError(f"Unknown symbol '{name}'.")


class PhysicsVerificationHelper:
    """
    Main helper class for physics verification tests.
    Provides common verification patterns and dimensional analysis tools.
    Supports multiple unit systems and smart symbol lookup.
    """

    def __init__(self, section_name: str, description: str = "",
                 unit_system: UnitSystem = UnitSystem.SI, quiet: bool = None):
        """
        Initialize verification helper for a specific section.

        Args:
            section_name: Name/number of the section being tested
            description: Optional detailed description
            unit_system: Unit system to use (SI, Gaussian, Heaviside-Lorentz)
            quiet: If True, suppress informational output. If None, check command line for --quiet flag.
        """
        self.section_name = section_name
        self.description = description
        self.unit_system = unit_system

        # Auto-detect quiet mode from command line if not explicitly set
        if quiet is None:
            self.quiet = _parse_global_quiet_flag()
        else:
            self.quiet = quiet

        self.debug_mode = False  # Can be enabled for extra verbose output
        self.results = []
        self.failed_checks = []

        # Initialize fundamental dimension symbols
        self.L = symbols('L', positive=True)  # Length
        self.M = symbols('M', positive=True)  # Mass
        self.T = symbols('T', positive=True)  # Time
        self.Q = symbols('Q', positive=True)  # Charge
        self.Theta = symbols('Theta', positive=True)  # Temperature

        # Initialize standard dimensions dictionary
        self.dims = self._initialize_standard_dimensions()

        # Initialize equation prefactors for unit system
        self.prefactors = self._get_unit_system_prefactors()

        # Print header and run self-tests (unless in quiet mode)
        if not self.quiet:
            self.print_header()
            self._run_self_tests()

    def _get_unit_system_prefactors(self) -> Dict[str, Any]:
        """Get equation prefactors for the current unit system with proper dimensions."""
        if self.unit_system != UnitSystem.SI:
            warnings.warn("Dimensional checks are anchored to SI; "
                         "Gaussian/HL change equation prefactors only.")

        if self.unit_system == UnitSystem.SI:
            return {
                # Maxwell equation prefactors with dimensions
                'gauss_E_prefactor': 1 / self.dims['epsilon_0'],      # ∇·E = (1/ε₀)ρ
                'gauss_B_prefactor': 0,                               # ∇·B = 0
                'faraday_prefactor': -1,                              # ∇×E = -∂B/∂t
                'ampere_j_prefactor': self.dims['mu_0'],              # ∇×B = μ₀j + μ₀ε₀∂E/∂t
                'ampere_E_prefactor': self.dims['mu_0'] * self.dims['epsilon_0'],
                'c_factor': 1,                                        # No explicit c in SI
                'gem_time_factor': 1 / self.dims['c'],                # GEM: using Ag0 = Φg/c
                # Legacy flags
                'use_epsilon0': True,
                'use_mu0': True,
                'c_in_equations': False,
            }
        elif self.unit_system == UnitSystem.GAUSSIAN:
            return {
                # Maxwell equation prefactors with dimensions
                'gauss_E_prefactor': 4 * pi,                          # ∇·E = 4πρ
                'gauss_B_prefactor': 0,
                'faraday_prefactor': -1 / self.dims['c'],             # ∇×E = -(1/c)∂B/∂t
                'ampere_j_prefactor': 4 * pi / self.dims['c'],        # ∇×B = (4π/c)j + (1/c)∂E/∂t
                'ampere_E_prefactor': 1 / self.dims['c'],
                'c_factor': 1 / self.dims['c'],                       # Explicit c factor
                'gem_time_factor': 1 / self.dims['c'],                # GEM: using Ag0 = Φg/c
                # Legacy flags
                'use_epsilon0': False,
                'use_mu0': False,
                'c_in_equations': True,
            }
        elif self.unit_system == UnitSystem.HEAVISIDE_LORENTZ:
            return {
                # Maxwell equation prefactors with dimensions
                'gauss_E_prefactor': 1,                               # ∇·E = ρ
                'gauss_B_prefactor': 0,
                'faraday_prefactor': -1 / self.dims['c'],             # ∇×E = -(1/c)∂B/∂t
                'ampere_j_prefactor': 1,                              # ∇×B = j + (1/c)∂E/∂t
                'ampere_E_prefactor': 1 / self.dims['c'],
                'c_factor': 1 / self.dims['c'],                       # Explicit c factor
                'gem_time_factor': 1 / self.dims['c'],                # GEM: using Ag0 = Φg/c
                # Legacy flags
                'use_epsilon0': False,
                'use_mu0': False,
                'c_in_equations': True,
            }
        else:
            raise ValueError(f"Unknown unit system: {self.unit_system}")

    def _initialize_standard_dimensions(self) -> Dict[str, Any]:
        """Initialize the comprehensive dimensions dictionary."""
        dims = {}

        # Fundamental constants
        dims['c'] = self.L / self.T                                    # Speed of light
        dims['G'] = self.L**3 / (self.M * self.T**2)                  # Gravitational constant
        dims['hbar'] = self.M * self.L**2 / self.T                    # Reduced Planck constant
        dims['h'] = self.M * self.L**2 / self.T                       # Planck constant
        dims['e'] = self.Q                                            # Elementary charge
        dims['epsilon_0'] = self.M**(-1) * self.L**(-3) * self.T**2 * self.Q**2  # Vacuum permittivity
        dims['mu_0'] = self.M * self.L * self.Q**(-2)  # Vacuum permeability
        dims['Z_0'] = (dims['mu_0'] / dims['epsilon_0'])**(Rational(1,2))  # Vacuum impedance = sqrt(mu_0/epsilon_0)
        dims['k_B'] = self.M * self.L**2 * self.T**(-2) * self.Theta**(-1)  # Boltzmann constant
        dims['Phi0'] = self.M * self.L**2 * self.T**(-1) * self.Q**(-1)   # Magnetic flux quantum = h/(2e)

        # Define voltage early (needed for derived EM quantities)
        dims['V_volt'] = self.M * self.L**2 * self.T**(-2) * self.Q**(-1)  # Voltage

        # Mathematical operators
        dims['nabla'] = self.L**(-1)                                  # Gradient operator
        dims['div'] = self.L**(-1)                                    # Divergence
        dims['curl'] = self.L**(-1)                                   # Curl
        dims['laplacian'] = self.L**(-2)                              # Laplacian
        dims['delta3'] = self.L**(-3)                                 # 3D Dirac delta
        dims['delta2'] = self.L**(-2)                                 # 2D Dirac delta
        dims['delta4'] = self.L**(-4)                                 # 4D Dirac delta
        dims['delta_t'] = self.T**(-1)                                # Time Dirac delta

        # Geometric elements
        dims['dl'] = self.L                                           # Line element
        dims['dA'] = self.L**2                                        # Area element
        dims['dV'] = self.L**3                                        # Volume element
        dims['theta'] = 1                                             # Angle (dimensionless)
        dims['phi'] = 1                                               # Angle φ (dimensionless)

        # ==== 4-VECTORS AND THEIR COMPONENTS ====

        # EM 4-current J^μ = (cρ_charge, j_current) - FIXED TO CHARGE DIMENSIONS
        dims['J_mu'] = self.Q * self.L**(-2) * self.T**(-1)           # EM 4-current
        dims['J0'] = self.Q * self.L**(-2) * self.T**(-1)             # J^0 = cρ_charge
        dims['J'] = self.Q * self.L**(-2) * self.T**(-1)              # J^i (spatial current)
        dims['J1'] = dims['J']                                        # J^1 component
        dims['J2'] = dims['J']                                        # J^2 component
        dims['J3'] = dims['J']                                        # J^3 component

        # EM 4-potential A^μ = (Φ/c, A) - FIXED A0 TO BE Φ/c
        dims['A'] = dims['V_volt'] * self.T / self.L                 # Vector potential = V*s/m
        dims['A_mu'] = dims['A']                                     # 4-potential
        dims['A0'] = dims['A']                                       # A^0 = Φ/c (same dims as spatial A)
        dims['A1'] = dims['A']                                        # A^1 component
        dims['A2'] = dims['A']                                        # A^2 component
        dims['A3'] = dims['A']                                        # A^3 component

        # Gravitational 4-potential (GEM)
        dims['Ag_mu'] = self.L / self.T                               # Grav 4-potential
        dims['Ag0'] = self.L / self.T                                 # Ag^0 = Φ_g/c
        dims['Ag'] = self.L / self.T                                  # Ag^i (spatial)
        dims['Ag1'] = dims['Ag']
        dims['Ag2'] = dims['Ag']
        dims['Ag3'] = dims['Ag']

        # ==== 3D SLICE OBSERVABLES (keep simple names) ====

        # Electromagnetic fields and potentials
        dims['E'] = dims['V_volt'] / self.L                           # Electric field = V/m
        dims['B'] = dims['V_volt'] * self.T / self.L**2               # Magnetic field = V*s/m²
        dims['D'] = self.Q / self.L**2                                # Electric displacement
        dims['H'] = self.Q * self.T**(-1) * self.L**(-1)              # Magnetic field intensity
        dims['Phi'] = dims['V_volt']                                  # Electric potential = voltage
        dims['Phi_B'] = self.M * self.L**2 * self.T**(-1) * self.Q**(-1)  # Magnetic flux

        # Keep old aliases for backward compatibility
        dims['E_field'] = dims['E']
        dims['B_field'] = dims['B']
        dims['D_field'] = dims['D']
        dims['H_field'] = dims['H']
        dims['Phi_E'] = dims['Phi']

        # Charge and current densities (3D)
        dims['rho'] = self.M / self.L**3                              # Mass density (default)
        dims['rho_charge'] = self.Q * self.L**(-3)                    # Charge density
        dims['j_mass'] = self.M / (self.L**2 * self.T)               # Mass current density
        dims['j_current'] = self.Q * self.T**(-1) * self.L**(-2)      # Electric current density
        dims['j'] = dims['j_current']                                 # j refers to electric current by default

        # Gravitational/GEM fields (3D)
        dims['Phi_g'] = self.L**2 / self.T**2                         # Gravitoelectric potential
        dims['A_g'] = self.L / self.T                                 # Gravitomagnetic potential
        dims['E_g'] = self.L / self.T**2                              # Gravitoelectric field
        dims['B_g'] = self.T**(-1)                                    # Gravitomagnetic field
        dims['g'] = self.L / self.T**2                                # Gravitational acceleration

        # ==== VELOCITY AND FLOW FIELDS ====

        dims['v'] = self.L / self.T                                   # Velocity field
        dims['v_L'] = self.L / self.T                                 # Bulk/sheet velocity
        dims['v_eff'] = self.L / self.T                               # Effective velocity
        dims['v_theta'] = self.L / self.T                             # Azimuthal velocity
        dims['a'] = self.L / self.T**2                                # Acceleration

        # Velocity components
        dims['v_x'] = self.L / self.T
        dims['v_y'] = self.L / self.T
        dims['v_z'] = self.L / self.T
        dims['v_w'] = self.L / self.T
        dims['delta_v_x'] = self.L / self.T
        dims['delta_v_y'] = self.L / self.T
        dims['delta_v_z'] = self.L / self.T
        dims['delta_v_w'] = self.L / self.T

        # RENAMED FROM A_x/A_y/A_z to avoid confusion with EM vector potential
        dims['u_x'] = self.L / self.T                                 # Fluid velocity x
        dims['u_y'] = self.L / self.T                                 # Fluid velocity y
        dims['u_z'] = self.L / self.T                                 # Fluid velocity z

        # ==== SPATIAL COORDINATES AND SCALES ====

        dims['x'] = self.L
        dims['y'] = self.L
        dims['z'] = self.L
        dims['w'] = self.L                                            # Extra dimension (NOT frequency!)
        dims['r'] = self.L                                            # Radial coordinate
        dims['r_perp'] = self.L                                       # Perpendicular radius
        dims['rho_cyl'] = self.L                                      # Cylindrical radius
        dims['xi'] = self.L                                           # Core radius/healing length
        dims['lambda'] = self.L                                       # Wavelength
        dims['lambda_C'] = self.L                                     # Compton wavelength
        dims['lambda_abs'] = self.L                                   # Absorption length
        dims['L_univ'] = self.L                                       # Universal length scale
        dims['L_scale'] = self.L                                      # Generic length scale

        # ==== ENERGY AND THERMODYNAMICS ====

        dims['E_energy'] = self.M * self.L**2 / self.T**2             # Energy
        dims['E_bend'] = self.M * self.L**2 / self.T**2               # Bending energy
        dims['K_bend'] = self.M * self.L**3 / self.T**2               # Bending modulus
        dims['V_energy'] = self.M * self.L**2 / self.T**2             # Potential energy
        dims['U_pot'] = self.M * self.L**2 / self.T**2               # Potential energy (renamed from V)
        dims['V'] = dims['V_volt']                                    # V now refers to voltage by default
        dims['S'] = self.M * self.L**2 / self.T                       # Action
        dims['Delta_E'] = self.M * self.L**2 / self.T**2              # Energy barrier
        dims['P'] = self.M / (self.L * self.T**2)                     # Pressure (3D)
        dims['T_surface'] = self.M / self.T**2                        # Surface tension
        dims['sigma_surface'] = self.M / self.L**2                    # Surface mass density
        dims['S_poynting'] = self.M * self.T**(-3)                    # Poynting vector
        dims['EMF'] = dims['V_volt']  # Electromotive force = voltage

        # ==== DENSITY VARIANTS (with clearer naming) ====

        # 4D densities (use _4 suffix for 4-volume densities)
        dims['rho_4'] = self.M / self.L**4                            # 4D mass density
        dims['rho_4_bg'] = self.M / self.L**4                         # Background 4D density
        dims['delta_rho_4'] = self.M / self.L**4                      # 4D density perturbation


        # 3D densities (no suffix needed, these are default)
        dims['rho_0'] = self.M / self.L**3                            # Background 3D density
        dims['rho_body'] = self.M / self.L**3                         # Body/matter density
        dims['rho_3D'] = self.M / self.L**3                           # Explicit 3D density
        dims['rho_avg'] = self.M / self.L**3                          # Average density
        dims['rho_cosmo'] = self.M / self.L**3                        # Cosmic density
        dims['rho_m'] = self.M / self.L**3                            # Mass density
        dims['rho_e'] = self.Q * self.L**(-3)                         # Charge density

        # ==== QUANTUM AND WAVE QUANTITIES ====

        dims['psi'] = self.L**(-Rational(3, 2))                       # Wavefunction (standard 3D)
        dims['Psi_GP'] = self.L**(-Rational(3, 2))                    # GP order parameter (3D convention)
        dims['Psi_GP_4D'] = self.L**(-2)                             # 4D GP order parameter (|Psi|^2 ~ L^-4)
        dims['g_GP_4D'] = self.M * self.L**6 / self.T**2             # 4D GP coupling constant
        # Note: 4D symbols added per FIXES.md for dimensional consistency
        dims['omega'] = self.T**(-1)                                  # Angular frequency
        dims['omega_freq'] = self.T**(-1)                             # Angular frequency (alias)
        dims['f'] = self.T**(-1)                                      # Frequency
        dims['nu'] = self.T**(-1)                                     # Alternative frequency
        dims['k'] = self.L**(-1)                                      # Wavenumber
        dims['p'] = self.M * self.L / self.T                          # Momentum
        dims['wavelength'] = self.L                                   # Wavelength

        # ==== VORTEX AND CIRCULATION ====

        dims['Gamma'] = self.L**2 / self.T                            # Circulation quantum
        dims['kappa'] = self.L**2 / self.T                            # Quantum of circulation
        dims['M_dot'] = self.M / self.T                               # Sink strength/mass flow
        dims['M_dot_i'] = self.M / self.T                             # Individual sink strength
        dims['M_dot_body'] = self.M / self.T                          # Body sink coefficient (for use with δ³)
        dims['M_dot_density'] = self.M / (self.L**3 * self.T)        # Distributed sink density (without δ³)

        # ==== MATERIAL PROPERTIES ====

        dims['sigma_conductivity'] = self.M**(-1) * self.L**(-3) * self.T * self.Q**2  # T^1
        dims['sigma_kernel'] = self.L                                 # Gaussian kernel width
        dims['epsilon_material'] = self.M**(-1) * self.L**(-3) * self.T**2 * self.Q**2
        dims['epsilon'] = dims['epsilon_material']                    # Alias
        dims['mu_material'] = self.M * self.L * self.Q**(-2)
        dims['mu'] = dims['mu_material']                              # Alias

        # Dimensionless material parameters
        dims['epsilon_r'] = 1                                         # Relative permittivity
        dims['mu_r'] = 1                                              # Relative permeability
        dims['n_refr'] = 1                                            # Refractive index

        # ==== MISCELLANEOUS ====

        dims['m'] = self.M                                            # Mass
        dims['m_core'] = self.M / self.L**2                           # Core sheet density
        dims['J_angular'] = self.M * self.L**2 / self.T               # Angular momentum
        dims['g_GP'] = self.M * self.L**5 / self.T**2                 # GP interaction
        dims['mathcal_L'] = self.M * self.L**(-1) * self.T**(-2)      # Lagrangian density
        dims['gamma'] = self.T**(-1)                                  # Dissipation rate
        dims['t'] = self.T                                            # Time coordinate
        dims['A_core'] = self.L**2                                    # Core area
        dims['dA_w'] = self.L**2                                      # w-axis area element
        dims['kappa_geom'] = self.L**(-1)                            # Geometric curvature

        # Dimensionless quantities
        dims['epsilon_xi'] = 1                                        # Thinness parameter ξ/ρ
        dims['epsilon_kappa'] = 1                                     # Flatness parameter κρ
        dims['gamma_lorentz'] = 1                                     # Lorentz factor
        dims['beta'] = 1                                              # v/c
        dims['n_quantum'] = 1                                         # Quantum number
        dims['phi_golden'] = 1                                        # Golden ratio

        # Circuit quantities (optional, for EM tests)
        dims['I_current'] = self.Q / self.T                           # Current
        dims['I'] = dims['I_current']                                 # Alias for backward compatibility
        dims['R'] = dims['V_volt'] / dims['I']                        # Resistance = V/I
        dims['C_cap'] = self.Q**2 / dims['V_energy']                 # Capacitance = Q^2/Energy
        dims['L_ind'] = dims['V_volt'] * self.T / dims['I']           # Inductance = V*t/I

        # Chiral parameters
        dims['Omega_0'] = self.T**(-1)                               # Chiral coupling strength
        dims['tau_twist'] = self.L**(-1)                             # Twist density
        dims['theta_twist'] = 1                                      # Twist angle (dimensionless)

        # Membrane properties
        dims['n_bar'] = self.L**(-3)                                 # Membrane density
        dims['m_bar'] = self.M                                       # Average deficit mass

        # Additional legacy field theory potentials
        dims['Phi_4D'] = self.L**2 / self.T                          # 4D scalar velocity potential
        dims['B4_x'] = self.L**2 / self.T                            # 4D vector potential x
        dims['B4_y'] = self.L**2 / self.T                            # 4D vector potential y
        dims['B4_z'] = self.L**2 / self.T                            # 4D vector potential z
        dims['Phi_bar'] = self.L**3 / self.T                         # Integrated scalar potential
        dims['B4_x_bar'] = self.L**3 / self.T                        # Integrated vector x
        dims['B4_y_bar'] = self.L**3 / self.T                        # Integrated vector y
        dims['B4_z_bar'] = self.L**3 / self.T                        # Integrated vector z
        dims['Psi_scalar'] = self.L**2 / self.T**2                   # 3D scalar field potential
        dims['Psi_3D'] = self.L**2 / self.T**2                       # Alternative notation
        dims['Psi_field'] = self.L**2 / self.T**2                    # Field potential (avoid wavefunction collision)
        dims['Psi_global'] = self.L**2 / self.T**2                   # Global potential
        # Note: Use 'psi' (lowercase) for wavefunction to distinguish from field potentials

        # Mass current components (renamed to avoid confusion with EM J_mu)
        dims['Jm_x'] = dims['j_mass']                                # Mass current x component
        dims['Jm_y'] = dims['j_mass']                                # Mass current y component
        dims['Jm_z'] = dims['j_mass']                                # Mass current z component

        dims['A_em'] = dims['A']  # EM vector potential
        dims['P_4D'] = self.M / (self.L**2 * self.T**2)              # 4D pressure
        dims['delta_P'] = self.M / (self.L**2 * self.T**2)           # Pressure perturbation

        # 4D specific coordinates
        dims['r_4'] = self.L                                          # 4D radius
        dims['v_4D_x'] = self.L / self.T                             # 4D velocity components
        dims['v_4D_y'] = self.L / self.T
        dims['v_4D_z'] = self.L / self.T
        dims['v_4D_w'] = self.L / self.T

        # ==== NICE-TO-HAVE: EM ENERGY AND POWER DENSITIES ====

        dims['u_EM'] = self.M * self.L**(-1) * self.T**(-2)           # EM energy density: (ε₀E² + B²/μ₀)/2
        dims['u_4D'] = self.M * self.L**(-2) * self.T**(-2)           # 4D energy density (per 4D volume)
        dims['u_3D'] = self.M * self.L**(-1) * self.T**(-2)           # 3D energy density (per 3D volume)
        dims['qdot_Joule'] = self.M * self.L**(-1) * self.T**(-3)     # Joule heating: j·E
        dims['f_Lorentz'] = self.M * self.L**(-2) * self.T**(-2)      # Lorentz force density: ρE + j×B

        # ==== NICE-TO-HAVE: STRESS-ENERGY TENSOR COMPONENTS ====

        dims['T00'] = self.M * self.L**(-1) * self.T**(-2)            # Energy density T⁰⁰
        dims['T0i'] = self.M * self.L**(-2) * self.T**(-1)            # Momentum density T⁰ⁱ
        dims['Tij'] = self.M * self.L**(-1) * self.T**(-2)            # Stress tensor Tⁱʲ
        dims['kappa_Einstein'] = self.T**2 / (self.M * self.L)        # Einstein coupling: 8πG/c⁴

        # ==== NICE-TO-HAVE: ADDITIONAL EM/GRAVITY QUANTITIES ====

        dims['Maxwell_stress'] = self.M * self.L**(-1) * self.T**(-2) # Maxwell stress tensor
        dims['E_field_energy'] = self.M * self.L**(-1) * self.T**(-2) # Electric field energy density
        dims['B_field_energy'] = self.M * self.L**(-1) * self.T**(-2) # Magnetic field energy density

        # Keep E_energy for energy, E for electric field
        # REMOVED: dims['E'] = ... that was overwriting electric field

        return dims

    def _run_self_tests(self):
        """Run internal self-tests for fundamental relationships."""
        self.info("Running dimensional self-tests...")

        # Test fundamental speed of light relationship
        c_from_constants = 1 / sp.sqrt(self.dims['mu_0'] * self.dims['epsilon_0'])
        self.check_dims("c from μ₀ε₀", self.dims['c'], c_from_constants, verbose=False)

        # Test impedance relationship
        Z0_from_ratio = sp.sqrt(self.dims['mu_0'] / self.dims['epsilon_0'])
        self.check_dims("Z₀ = √(μ₀/ε₀)", self.dims['Z_0'], Z0_from_ratio, verbose=False)

        # Test voltage-current relationship for impedance
        Z0_from_VI = self.dims['V_volt'] / self.dims['I']
        self.check_dims("Z₀ as V/I", self.dims['Z_0'], Z0_from_VI, verbose=False)

        # Test D = εE relationship
        D_from_eps_E = self.dims['epsilon_0'] * self.dims['E']
        self.check_dims("D = ε₀E", self.dims['D'], D_from_eps_E, verbose=False)

        # Test B = μH relationship
        B_from_mu_H = self.dims['mu_0'] * self.dims['H']
        self.check_dims("B = μ₀H", self.dims['B'], B_from_mu_H, verbose=False)

        # Test capacitance as Q²/Energy
        C_dimensional = self.Q**2 / self.dims['E_energy']
        self.check_dims("C = Q²/U", self.dims['C_cap'], C_dimensional, verbose=False)

        # Test inductance as V·t/I
        L_dimensional = self.dims['V_volt'] * self.dims['t'] / self.dims['I']
        self.check_dims("L = V·t/I", self.dims['L_ind'], L_dimensional, verbose=False)

        # Test GEM relationships
        self.check_dims("GEM: B_g = curl A_g", self.dims['B_g'], self.curl_dim(self.dims['A_g']), verbose=False)
        self.check_dims("GEM: E_g from gradient", self.dims['E_g'], self.grad_dim(self.dims['Phi_g']), verbose=False)
        self.check_dims("GEM: E_g from time derivative", self.dims['E_g'], self.dt(self.dims['A_g']), verbose=False)

        self.info("Self-tests completed.\n")

    def get_dim(self, name: str, suggest: bool = True, strict: bool = False) -> Any:
        """
        Get dimension for a symbol with fuzzy matching suggestions.

        Args:
            name: Symbol name
            suggest: Whether to provide suggestions if not found
            strict: If True, disallow aliases and require canonical names

        Returns:
            Dimension expression

        Raises:
            UnknownSymbolError: If symbol not found
        """
        if strict:
            # In strict mode, only allow exact matches, no aliases
            if name in self.dims:
                return self.dims[name]
            if suggest:
                raise UnknownSymbolError(f"Unknown symbol '{name}' (strict mode - no aliases allowed)")
            raise UnknownSymbolError(f"Unknown symbol '{name}'")

        return SymbolRegistry.get_symbol(name, self.dims, suggest)

    # ==== PRINT UTILITY METHODS ====

    def info(self, message: str):
        """Print informational message (respects quiet mode)."""
        if not self.quiet:
            print(message)

    def success(self, message: str):
        """Print success message (respects quiet mode)."""
        if not self.quiet:
            print(f"✓ {message}")

    def warning(self, message: str):
        """Print warning message (always prints)."""
        print(f"⚠️ {message}")

    def error(self, message: str):
        """Print error message (always prints)."""
        print(f"✗ {message}")

    def debug(self, message: str):
        """Print debug message (only in verbose/debug mode)."""
        if not self.quiet and self.debug_mode:
            print(f"🔍 {message}")

    def section_header(self, title: str):
        """Print section header (respects quiet mode)."""
        if not self.quiet:
            print(f"\n--- {title} ---")

    def add_dimension(self, name: str, dimension: Any, allow_overwrite: bool = False):
        """
        Add a custom dimension to the dictionary with overwrite protection.

        Args:
            name: Symbol name
            dimension: Dimensional expression using self.L, self.M, self.T, self.Q, self.Theta
            allow_overwrite: If False (default), raise error if name already exists

        Raises:
            KeyError: If name already exists and allow_overwrite is False
        """
        if not allow_overwrite and name in self.dims:
            raise KeyError(f"Dimension '{name}' already defined. Use allow_overwrite=True to replace.")
        self.dims[name] = dimension

    def add_dimensions(self, new_dims: Dict[str, Any], allow_overwrite: bool = False):
        """
        Add multiple custom dimensions to the dictionary with overwrite protection.

        Args:
            new_dims: Dictionary of {name: dimension} pairs
            allow_overwrite: If False (default), raise error if any name already exists
        """
        for name, dimension in new_dims.items():
            self.add_dimension(name, dimension, allow_overwrite=allow_overwrite)

    def declare_dimensionless(self, *names: str):
        """
        Declare symbols as dimensionless to prevent dimensional analysis errors.

        Args:
            *names: Variable names to declare as dimensionless
        """
        for name in names:
            self.dims[name] = 1

    def dt(self, dim: Any) -> Any:
        """Get dimensions of time derivative."""
        return dim / self.T

    def dtt(self, dim: Any) -> Any:
        """Get dimensions of second time derivative."""
        return dim / self.T**2

    def dx(self, dim: Any) -> Any:
        """Get dimensions of spatial derivative."""
        return dim / self.L

    def dxx(self, dim: Any) -> Any:
        """Get dimensions of second spatial derivative."""
        return dim / self.L**2

    def assert_dimensionless(self, expr: Any, where: str = "") -> bool:
        """
        Assert that an expression is dimensionless (for use in transcendental functions).

        Args:
            expr: Expression to check
            where: Optional context string for error message

        Returns:
            True if dimensionless

        Raises:
            DimensionalMismatchError: If expression is not dimensionless
        """
        # Handle symbolic constants like pi, E, etc. - they are dimensionless by definition
        symbolic_constants = {pi, sp.E, sp.I}
        if expr in symbolic_constants:
            return True

        dim = self._extract_dimensions(expr, strict=True)  # Use strict=True to catch real mistakes
        is_dimensionless = simplify(dim) == 1

        if not is_dimensionless:
            context = f" in {where}" if where else ""
            raise DimensionalMismatchError(f"Expected dimensionless argument{context}, got [{dim}]")

        return True

    def exp_dimless(self, x: Any, where: str = "exp") -> Any:
        """
        Safe exponential function that ensures argument is dimensionless.

        Args:
            x: Argument to exponential
            where: Context for error messages

        Returns:
            exp(x) if x is dimensionless

        Raises:
            DimensionalMismatchError: If x is not dimensionless
        """
        self.assert_dimensionless(x, where)
        return sp.exp(x)

    def sin_dimless(self, x: Any, where: str = "sin") -> Any:
        """
        Safe sine function that ensures argument is dimensionless.

        Args:
            x: Argument to sine
            where: Context for error messages

        Returns:
            sin(x) if x is dimensionless

        Raises:
            DimensionalMismatchError: If x is not dimensionless
        """
        self.assert_dimensionless(x, where)
        return sp.sin(x)

    def cos_dimless(self, x: Any, where: str = "cos") -> Any:
        """
        Safe cosine function that ensures argument is dimensionless.

        Args:
            x: Argument to cosine
            where: Context for error messages

        Returns:
            cos(x) if x is dimensionless

        Raises:
            DimensionalMismatchError: If x is not dimensionless
        """
        self.assert_dimensionless(x, where)
        return sp.cos(x)

    def log_dimless(self, x: Any, where: str = "log") -> Any:
        """
        Safe logarithm function that ensures argument is dimensionless.

        Args:
            x: Argument to logarithm
            where: Context for error messages

        Returns:
            log(x) if x is dimensionless

        Raises:
            DimensionalMismatchError: If x is not dimensionless
        """
        self.assert_dimensionless(x, where)
        return sp.log(x)

    def validate_transcendentals(self, expr: Any, where: str = "expression") -> bool:
        """
        Walk an expression tree and validate that transcendental functions have dimensionless arguments.

        Args:
            expr: Expression to validate
            where: Context for error messages

        Returns:
            True if all transcendental functions have dimensionless arguments

        Raises:
            DimensionalMismatchError: If any transcendental has dimensional argument
        """
        transcendental_funcs = {sp.exp, sp.log, sp.sin, sp.cos, sp.tan,
                               sp.sinh, sp.cosh, sp.tanh, sp.asin, sp.acos, sp.atan,
                               sp.erf, sp.erfc, sp.Heaviside}

        def check_node(node):
            if hasattr(node, 'func') and node.func in transcendental_funcs:
                if len(node.args) > 0:
                    arg = node.args[0]
                    try:
                        self.assert_dimensionless(arg, f"{node.func.__name__} in {where}")
                    except DimensionalMismatchError:
                        raise

            # Recursively check all arguments
            if hasattr(node, 'args'):
                for arg in node.args:
                    check_node(arg)

        try:
            check_node(expr)
            return True
        except DimensionalMismatchError:
            raise

    def grad_dim(self, scalar_dim: Any) -> Any:
        """Get dimensions of gradient of a scalar field."""
        return self.dims['nabla'] * scalar_dim

    def div_dim(self, vector_dim: Any) -> Any:
        """Get dimensions of divergence of a vector field."""
        return self.dims['div'] * vector_dim

    def curl_dim(self, vector_dim: Any) -> Any:
        """Get dimensions of curl of a vector field."""
        return self.dims['curl'] * vector_dim

    def lap_dim(self, scalar_dim: Any) -> Any:
        """Get dimensions of Laplacian of a scalar field."""
        return self.dims['laplacian'] * scalar_dim

    def energy_density_from_field(self, field_name: str, use_material: bool = False) -> Any:
        """
        Get energy density dimensions for electromagnetic fields.

        Args:
            field_name: Field name ('E' for electric, 'B' for magnetic)
            use_material: If True, use material properties (ε, μ) instead of vacuum (ε₀, μ₀)

        Returns:
            Energy density dimensions

        Raises:
            KeyError: If field_name is not recognized
        """
        if field_name == 'E':
            eps = self.dims['epsilon'] if use_material else self.dims['epsilon_0']
            return eps * self.dims['E']**2
        elif field_name == 'B':
            mu = self.dims['mu'] if use_material else self.dims['mu_0']
            return self.dims['B']**2 / mu
        else:
            raise KeyError(f"Unknown field '{field_name}' for energy density calculation. Use 'E' or 'B'.")

    def poisson_dim(self, phi_dim: Any, source_dim: Any, name: str = "Poisson") -> bool:
        """
        Convenience method for Poisson equation dimensional check: ∇²φ = source.

        Args:
            phi_dim: Dimensions of potential φ
            source_dim: Dimensions of source term
            name: Name for the check

        Returns:
            True if dimensionally consistent
        """
        return self.check_dims(f"{name} equation", self.lap_dim(phi_dim), source_dim)

    def check_dims(self, name: str, expr1: Any, expr2: Any,
                   record: bool = True, verbose: bool = True,
                   si_only: bool = False) -> bool:
        """
        Check dimensional consistency between two expressions.

        Args:
            name: Description of the check
            expr1: First expression
            expr2: Second expression
            record: Whether to record in results
            verbose: Whether to print output
            si_only: If True, skip check when not in SI system

        Returns:
            True if dimensions match, False otherwise
        """
        # Guard against dimensional checks in non-SI systems
        if si_only and self.unit_system != UnitSystem.SI:
            if verbose:
                print(f"~ {name} - SKIPPED (SI-only check in {self.unit_system.value})")
            if record:
                self.results.append((name + " [SKIPPED]", True))
            return True
        # Extract dimensional structure ignoring numerical coefficients
        dim1 = self._extract_dimensions(expr1, strict=True)
        dim2 = self._extract_dimensions(expr2, strict=True)

        # Check for dimension mismatch sentinel
        if (isinstance(dim1, sp.Symbol) and str(dim1) == 'DIM_MISMATCH' or
            isinstance(dim2, sp.Symbol) and str(dim2) == 'DIM_MISMATCH'):
            check = False
            if verbose:
                print(f"✗ {name} - DIMENSIONAL MISMATCH IN EXPRESSION")
        else:
            # Check consistency
            check = self._check_dimensional_consistency(dim1, dim2)

        if record:
            self.results.append((name, check))
            if not check:
                self.failed_checks.append(name)

        if verbose:
            if check:
                self.success(name.replace("✓ ", ""))
            else:
                self.error(f"{name}")
                if str(dim1) != 'DIM_MISMATCH':
                    print(f"    Expected: [{dim1}]")
                    print(f"    Got:      [{dim2}]")

        return check

    def check_zero(self, name: str, expr: Any, record: bool = True,
                   verbose: bool = True) -> bool:
        """
        Check that an expression is identically zero with dimensional validation.

        Args:
            name: Description of the check
            expr: Expression that should be zero
            record: Whether to record in results
            verbose: Whether to print output

        Returns:
            True if expression is zero, False otherwise
        """
        dim = self._extract_dimensions(expr, strict=True)

        # Check if it's actually zero
        try:
            is_zero = simplify(expr) == 0
        except:
            is_zero = (expr == 0)

        # If it's truly zero, that's fine
        if is_zero:
            if record:
                self.results.append((name, True))
            if verbose:
                self.success(name)
            return True

        # If it's not zero but has ZERO dimension marker, something's wrong
        if str(dim) == 'ZERO':
            if record:
                self.results.append((name, False))
                self.failed_checks.append(name)
            if verbose:
                self.error(f"{name} - Expression marked as zero but evaluates to: {expr}")
            return False

        # Non-zero expression - warn about dimensions
        if record:
            self.results.append((name, False))
            self.failed_checks.append(name)
        if verbose:
            self.error(f"{name} - Non-zero expression with dimensions [{dim}]")
            print(f"    Expression: {expr}")

        return False

    def check_eq(self, name: str, lhs: Any, rhs: Any,
                 record: bool = True, verbose: bool = True) -> bool:
        """
        Check equation equality with simplification.

        Args:
            name: Description of the check
            lhs: Left hand side
            rhs: Right hand side
            record: Whether to record in results
            verbose: Whether to print output

        Returns:
            True if equations are equal, False otherwise
        """
        # Guard against misusing check_eq for pure dimension objects
        combined_symbols = set()
        if hasattr(lhs, 'free_symbols'):
            combined_symbols.update(lhs.free_symbols)
        if hasattr(rhs, 'free_symbols'):
            combined_symbols.update(rhs.free_symbols)

        base_dim_strs = {str(self.M), str(self.L), str(self.T), str(self.Q), str(self.Theta)}
        if combined_symbols and all(str(s) in base_dim_strs for s in combined_symbols):
            warnings.warn(f"{name}: use check_dims() for pure-dimension comparisons")
        try:
            check = simplify(lhs - rhs) == 0
        except:
            # Fallback for complex expressions
            check = (lhs == rhs)

        if record:
            self.results.append((name, check))
            if not check:
                self.failed_checks.append(name)

        if verbose:
            if check:
                self.success(name)
            else:
                self.error(name)
                print(f"    LHS: {lhs}")
                print(f"    RHS: {rhs}")
                print(f"    Diff: {simplify(lhs - rhs)}")

        return check

    def check_limit(self, name: str, expr: Any, var: Any, point: Any,
                    expected: Any, record: bool = True, verbose: bool = True) -> bool:
        """
        Check limit calculation.

        Args:
            name: Description of the check
            expr: Expression to take limit of
            var: Variable approaching limit
            point: Point to approach (can be oo, -oo, or number)
            expected: Expected limit value
            record: Whether to record in results
            verbose: Whether to print output

        Returns:
            True if limit matches expected, False otherwise
        """
        try:
            result = limit(expr, var, point)
            check = (result == expected)
        except Exception as e:
            result = f"Error: {e}"
            check = False

        if record:
            self.results.append((name, check))
            if not check:
                self.failed_checks.append(name)

        if verbose:
            if check:
                self.success(name)
            else:
                self.error(name)
                print(f"    Limit: {result}")
                print(f"    Expected: {expected}")

        return check

    def check_asymptotic(self, name: str, exact_expr: Any, approx_expr: Any,
                         var: Any, point: Any = oo, record: bool = True,
                         verbose: bool = True) -> bool:
        """
        Check asymptotic expansion validity.

        Args:
            name: Description of the check
            exact_expr: Exact expression
            approx_expr: Asymptotic approximation
            var: Variable for limit
            point: Point to check asymptotic behavior (default: oo)
            record: Whether to record in results
            verbose: Whether to print output

        Returns:
            True if asymptotic expansion is valid, False otherwise
        """
        try:
            difference = exact_expr - approx_expr
            limit_diff = limit(difference, var, point)
            check = (limit_diff == 0)
        except:
            check = False
            limit_diff = "Could not compute"

        if record:
            self.results.append((name, check))
            if not check:
                self.failed_checks.append(name)

        if verbose:
            if check:
                self.success(name)
            else:
                self.error(name)
                print(f"    Limit of difference: {limit_diff}")

        return check

    def verify_and_record(self, description: str, check_result: bool,
                          details: str = "") -> bool:
        """
        Generic verification with recording (compatibility with existing tests).

        Args:
            description: Description of the check
            check_result: Boolean result
            details: Optional details string

        Returns:
            The check_result
        """
        self.results.append((description, check_result))
        if not check_result:
            self.failed_checks.append((description, details))

        if check_result:
            self.success(description)
        else:
            self.error(description)
            if details:
                print(f"    Details: {details}")

        return check_result

    def section(self, title: str, width: int = 60):
        """Print section header."""
        if not self.quiet:
            print(f"\n{'='*width}")
            print(title)
            print('='*width)

    def subsection(self, title: str, width: int = 50):
        """Print subsection header."""
        if not self.quiet:
            print(f"\n{title}")
            print('-'*width)

    def print_header(self):
        """Print initial header for the test."""
        print("="*80)
        print(self.section_name)
        if self.description:
            print(self.description)
        if self.unit_system != UnitSystem.SI:
            print(f"Unit System: {self.unit_system.value}")
        print("="*80)

    def summary(self, show_failed: bool = True) -> float:
        """
        Print verification summary.

        Args:
            show_failed: Whether to list failed checks

        Returns:
            Success rate as percentage
        """
        passed = sum(1 for _, result in self.results if result)
        total = len(self.results)
        rate = (passed / total * 100) if total > 0 else 0

        self.section(f"VERIFICATION SUMMARY: {passed}/{total} passed ({rate:.1f}%)")

        if rate == 100:
            if not self.quiet:
                print("\n🎉 ALL VERIFICATIONS PASSED! 🎉")
        elif rate >= 95:
            if not self.quiet:
                print("\n✅ VERIFICATION SUBSTANTIALLY COMPLETE")
                print(f"   {passed}/{total} checks passed")
        elif rate >= 85:
            if not self.quiet:
                print("\n⚠️ VERIFICATION MOSTLY SUCCESSFUL")
                print(f"   {passed}/{total} checks passed")
        else:
            print(f"\n❌ VERIFICATION NEEDS ATTENTION")
            print(f"   Only {passed}/{total} checks passed")

        if show_failed and self.failed_checks:
            print(f"\nFailed checks ({len(self.failed_checks)}):")
            for item in self.failed_checks:
                if isinstance(item, tuple):
                    name, details = item
                    print(f"  ✗ {name}")
                    if details:
                        print(f"      {details}")
                else:
                    print(f"  ✗ {item}")

        return rate

    def _is_base_dim(self, symbol: Any) -> bool:
        """Check if a symbol is one of the base dimensions {M,L,T,Q,Theta}."""
        if not hasattr(symbol, 'name'):
            return False
        base_names = {str(self.M), str(self.L), str(self.T), str(self.Q), str(self.Theta)}
        return str(symbol) in base_names

    def _extract_dimensions(self, expr: Any, strict: bool = True) -> Any:
        """
        Extract dimensional structure from expression, ignoring numerical coefficients.
        IMPROVED: Checks all terms in sums for consistency.

        Args:
            expr: Expression to analyze
            strict: If True, raise on unknown symbols; if False, warn and treat as dimensionless

        Returns:
            Dimensional structure (or DIM_MISMATCH sentinel if inconsistent)
        """
        if expr is None:
            return None

        # Handle numeric zero as dimension-agnostic
        if expr == 0:
            return sp.Symbol('ZERO')

        # Handle pure numbers (including all SymPy numeric types and Python numbers)
        if getattr(expr, "is_Number", False) or isinstance(expr, (int, float, complex)):
            return 1

        # For symbolic expressions, try to separate numerical and dimensional parts
        try:
            # Check for unknown symbols that aren't base dimensions FIRST
            if isinstance(expr, sp.Symbol) and not self._is_base_dim(expr):
                # Check if this symbol exists in our dimensions registry
                if str(expr) not in self.dims:
                    msg = f"Unknown symbol '{expr}' found in dimensional analysis. " \
                          f"Consider adding to dimensions registry or use base dimensions M,L,T,Q,Theta."
                    if strict:
                        raise DimensionalMismatchError(msg)
                    warnings.warn(msg)
                    return 1  # Treat as dimensionless with warning

            # List of dimensionless constants to ignore
            dimensionless = [pi, sp.E, sp.sqrt(2), sp.sqrt(5), sp.log(2),
                            sp.Rational(1,2), sp.Rational(3,2)]

            if hasattr(expr, 'is_Mul'):
                if expr.is_Mul:
                    dimensional_part = 1
                    for factor in expr.args:
                        if not (factor.is_number or factor in dimensionless):
                            dimensional_part *= factor
                    return dimensional_part
                elif expr.is_Pow:
                    base_dim = self._extract_dimensions(expr.base, strict)
                    exp_val = expr.exp
                    if not (exp_val.is_Number or (hasattr(exp_val, 'is_rational') and exp_val.is_rational)):
                        raise DimensionalMismatchError(f"Non-numeric exponent '{exp_val}' in dimensional expression.")
                    return base_dim ** exp_val
                elif expr.is_Add:
                    # IMPROVED: Check all terms have same dimensions
                    if len(expr.args) > 0:
                        base_dim = self._extract_dimensions(expr.args[0], strict)
                        for term in expr.args[1:]:
                            term_dim = self._extract_dimensions(term, strict)
                            if not self._check_dimensional_consistency(base_dim, term_dim):
                                # Return sentinel indicating dimension mismatch
                                return sp.Symbol('DIM_MISMATCH')
                        return base_dim
                    return 1

            return expr

        except DimensionalMismatchError:
            # Re-raise dimensional mismatch errors (don't swallow them)
            raise
        except:
            return expr

    def _check_dimensional_consistency(self, dim1: Any, dim2: Any) -> bool:
        """
        Check if two dimensional expressions are consistent.
        IMPROVED: Uses ratio test primarily, which is more robust.

        Args:
            dim1: First dimensional expression
            dim2: Second dimensional expression

        Returns:
            True if consistent, False otherwise
        """
        if dim1 is None or dim2 is None:
            return False

        # Handle zero as dimension-agnostic
        if str(dim1) == 'ZERO' or str(dim2) == 'ZERO':
            return True

        try:
            # Try direct comparison
            if dim1 == dim2:
                return True

            # Try ratio (more robust than difference)
            ratio = simplify(dim1 / dim2)
            if ratio == 1:
                return True

            # Check if both are dimensionless (equal to 1)
            if simplify(dim1) == 1 and simplify(dim2) == 1:
                return True

            return False

        except:
            # If symbolic manipulation fails, do string comparison as last resort
            return str(dim1) == str(dim2)


# ==============================================================================
# CONVENIENCE FUNCTIONS
# ==============================================================================

def quick_verify(name: str, condition: bool, details: str = "", helper=None, expected_failure: bool = False) -> bool:
    """
    Quick one-line verification without recording.

    Args:
        name: Description of check
        condition: Boolean condition
        details: Optional details
        helper: Optional PhysicsVerificationHelper instance for quiet mode support
        expected_failure: If True, failures are expected (diagnostic tests) and won't print in quiet mode

    Returns:
        The condition value
    """
    message = name
    if details:
        message += f": {details}"

    if helper is not None:
        # Use helper's print methods to respect quiet mode
        if condition:
            helper.success(message)
        else:
            # For expected failures in quiet mode, only show if not quiet or if it's an unexpected success
            if expected_failure and helper.quiet:
                # Don't print expected failures in quiet mode
                pass  
            else:
                helper.error(message)
    else:
        # Fallback to direct print (backward compatibility)
        status = "✓" if condition else "✗"
        print(f"{status} {message}")

    return condition


def batch_check_dims(helper: PhysicsVerificationHelper,
                     checks: List[Tuple[str, Any, Any]]):
    """
    Run multiple dimensional consistency checks.

    Args:
        helper: PhysicsVerificationHelper instance
        checks: List of (name, expr1, expr2) tuples
    """
    for check in checks:
        if len(check) == 3:
            name, expr1, expr2 = check
            helper.check_dims(name, expr1, expr2)
        else:
            print(f"Warning: Invalid check format: {check}")


def batch_check_eqs(helper: PhysicsVerificationHelper,
                    equations: List[Tuple[str, Any, Any]]):
    """
    Run multiple equation checks.

    Args:
        helper: PhysicsVerificationHelper instance
        equations: List of (name, lhs, rhs) tuples
    """
    for eq in equations:
        if len(eq) == 3:
            name, lhs, rhs = eq
            helper.check_eq(name, lhs, rhs)
        else:
            print(f"Warning: Invalid equation format: {eq}")


def define_symbols_batch(names: List[str], **kwargs) -> tuple:
    """
    Define multiple symbols at once.

    Args:
        names: List of symbol names
        **kwargs: Keyword arguments for symbols (real=True, positive=True, etc.)

    Returns:
        Tuple of symbols
    """
    return symbols(' '.join(names), **kwargs)


def print_dimensional_analysis(expression: Any, helper: PhysicsVerificationHelper,
                               name: str = "Expression"):
    """
    Print detailed dimensional analysis of an expression.

    Args:
        expression: Expression to analyze
        helper: PhysicsVerificationHelper instance for dimension lookup
        name: Name to display
    """
    dim = helper._extract_dimensions(expression)
    print(f"\nDimensional analysis of {name}:")
    print(f"  Expression: {expression}")
    print(f"  Dimensions: [{dim}]")
    print(f"  Simplified: [{simplify(dim)}]")


def create_section_verifier(section_name: str, description: str = "") -> PhysicsVerificationHelper:
    """
    Factory function to create a verifier with standard setup.

    Args:
        section_name: Section name/number
        description: Optional description

    Returns:
        Configured PhysicsVerificationHelper instance
    """
    return PhysicsVerificationHelper(section_name, description)


# ==============================================================================
# SPECIALIZED VERIFICATION PATTERNS
# ==============================================================================

def verify_wave_equation(helper: PhysicsVerificationHelper, name: str,
                         time_term: Any, space_term: Any, source_term: Any = None) -> bool:
    """
    Verify a wave equation of the form: ∂_tt φ - c²∇²φ = source.

    Args:
        helper: PhysicsVerificationHelper instance
        name: Name of the equation
        time_term: ∂_tt φ term dimensions
        space_term: c²∇²φ term dimensions
        source_term: Optional source term dimensions

    Returns:
        True if all terms are dimensionally consistent
    """
    helper.subsection(f"{name} Wave Equation Verification")

    # Check time and space terms match
    time_space_check = helper.check_dims(
        f"{name}: time-space consistency",
        time_term, space_term
    )

    # If source term provided, check it matches too
    if source_term is not None:
        source_check = helper.check_dims(
            f"{name}: source term consistency",
            time_term, source_term
        )
        return time_space_check and source_check

    return time_space_check


def verify_conservation_law(helper: PhysicsVerificationHelper, name: str,
                           density_rate: Any, flux_div: Any, source: Any = None) -> bool:
    """
    Verify a conservation law: ∂_t ρ + ∇·J = source.

    Args:
        helper: PhysicsVerificationHelper instance
        name: Name of the conservation law
        density_rate: ∂_t ρ term dimensions
        flux_div: ∇·J term dimensions
        source: Optional source term dimensions

    Returns:
        True if conservation law is dimensionally consistent
    """
    helper.subsection(f"{name} Conservation Law")

    # Check density rate and flux divergence match
    continuity_check = helper.check_dims(
        f"{name}: continuity",
        density_rate, flux_div
    )

    # If source provided, check it too
    if source is not None:
        source_check = helper.check_dims(
            f"{name}: source term",
            density_rate, source
        )
        return continuity_check and source_check

    return continuity_check


def verify_poisson_em(helper: PhysicsVerificationHelper) -> bool:
    """
    Verify electromagnetic Poisson equations: ∇·E = ρ/ε₀ and ∇²Φ = -ρ/ε₀.

    Args:
        helper: PhysicsVerificationHelper instance

    Returns:
        True if both EM Poisson relationships are dimensionally consistent
    """
    helper.subsection("EM Poisson Equations")

    # Gauss law: ∇·E = ρ/ε₀
    gauss_check = helper.check_dims(
        "EM Gauss law: ∇·E = ρ/ε₀",
        helper.div_dim(helper.dims['E']),
        helper.dims['rho_charge'] / helper.dims['epsilon_0']
    )

    # Poisson equation: ∇²Φ = -ρ/ε₀ (sign ignored for dimensional analysis)
    poisson_check = helper.check_dims(
        "EM Poisson: ∇²Φ = -ρ/ε₀",
        helper.lap_dim(helper.dims['Phi']),
        helper.dims['rho_charge'] / helper.dims['epsilon_0']
    )

    return gauss_check and poisson_check


def verify_poisson_grav(helper: PhysicsVerificationHelper) -> bool:
    """
    Verify gravitational Poisson equation: ∇²Φ_g = 4πGρ.

    Args:
        helper: PhysicsVerificationHelper instance

    Returns:
        True if gravitational Poisson equation is dimensionally consistent
    """
    helper.subsection("Gravitational Poisson Equation")

    # Gravitational Poisson: ∇²Φ_g = 4πGρ (4π ignored for dimensional analysis)
    return helper.check_dims(
        "Grav Poisson: ∇²Φ_g = 4πGρ",
        helper.lap_dim(helper.dims['Phi_g']),
        helper.dims['G'] * helper.dims['rho']
    )


def verify_poisson_equation(helper: PhysicsVerificationHelper, name: str,
                           laplacian_term: Any, source_term: Any) -> bool:
    """
    Verify a Poisson equation: ∇²φ = source.

    Args:
        helper: PhysicsVerificationHelper instance
        name: Name of the equation
        laplacian_term: ∇²φ term dimensions
        source_term: Source term dimensions

    Returns:
        True if Poisson equation is dimensionally consistent
    """
    helper.subsection(f"{name} Poisson Equation")

    return helper.check_dims(
        f"{name}: Poisson consistency",
        laplacian_term, source_term
    )


# ==============================================================================
# TEST TEMPLATE GENERATOR
# ==============================================================================

def generate_test_template(section_number: str, section_title: str) -> str:
    """
    Generate a template for a new test file.

    Args:
        section_number: Section number (e.g., "2.3")
        section_title: Section title

    Returns:
        String containing template code
    """
    template = f'''"""
Section {section_number}: {section_title} - Verification
{"="*60}

Complete verification of all mathematical relationships in Section {section_number}.
"""

from helpers import (PhysicsVerificationHelper, quick_verify,
                     define_symbols_batch, batch_check_dims,
                     verify_wave_equation, verify_conservation_law,
                     verify_poisson_em, verify_poisson_grav)
import sympy as sp
from sympy import symbols, simplify, diff, integrate, limit, oo

# Initialize verification helper
v = PhysicsVerificationHelper(
    "Section {section_number}: {section_title}",
    "Complete mathematical verification"
)

# Define section-specific symbols
t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)

# Add any custom dimensions needed for this section
v.add_dimensions({{
    'custom_quantity': v.M * v.L / v.T**2,  # Example
}})

v.section("EQUATION 1: EXAMPLE VERIFICATION")

# Verification 1: Dimensional consistency
v.subsection("Dimensional Analysis")

lhs = v.get_dim('rho_4') / v.get_dim('t')
rhs = v.get_dim('rho_4') * v.get_dim('v') / v.get_dim('r')

v.check_dims("Equation 1 dimensional consistency", lhs, rhs)

# Verification 2: Mathematical relationship
v.subsection("Mathematical Verification")

expr1 = sp.sqrt(v.get_dim('G') * v.get_dim('rho_0'))  # This has dimension 1/T
expr2 = 1 / v.get_dim('t')  # Inverse time scale

v.check_eq("Dynamical time scale", expr1, expr2)

# Generate summary
success_rate = v.summary()

if success_rate == 100:
    print("\\n✅ Section {section_number} fully verified!")
'''

    return template


if __name__ == "__main__":
    # Example usage demonstration
    print("Physics Verification Helper Library")
    print("="*50)
    print("\nAvailable classes:")
    print("  - PhysicsVerificationHelper")
    print("  - UnitSystem (SI, GAUSSIAN, HEAVISIDE_LORENTZ)")
    print("  - SymbolRegistry")
    print("\nDimensions available:", len(PhysicsVerificationHelper("test").dims))
    print("\nReady for accurate physics verification!")
