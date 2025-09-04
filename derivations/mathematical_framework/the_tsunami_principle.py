#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verification tests for "The Tsunami Principle" subsection.

This module implements dimensional and mathematical verification for the core
bulk tsunami physics documented in "The Tsunami Principle" subsection of the
mathematical framework document.

Verifies:
- Bulk tsunami wave equation: ∂²_t δρ - v²_L ∇²₄ δρ = -M δ⁴(r₄) δ'(t)
- Retarded Green's function solution: δρ = M G^ret_(4)(r₄,t;v_L)
"""

import os
import sys

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem
)

def test_bulk_tsunami_wave_equation(v):
    """Test the core bulk tsunami wave equation from the paper."""
    # Paper equation: ∂²_t δρ(r₄,t) - v²_L ∇²₄ δρ(r₄,t) = -M δ⁴(r₄) δ'(t)

    # LHS terms should match
    time_term = v.dtt(v.get_dim('delta_rho_4'))
    space_term = (v.get_dim('v_L')**2) * v.lap_dim(v.get_dim('delta_rho_4'))
    v.check_dims("Bulk tsunami: time vs space terms", time_term, space_term)

    # RHS source term: [M δ⁴(r₄) δ'(t)] = M * L^-4 * T^-2 = M L^-4 T^-2
    source_term = v.get_dim('M_sink') * v.get_dim('delta4') * (v.T**-2)  # δ'(t) = T^-2
    v.check_dims("Bulk tsunami: LHS vs source", time_term, source_term)

    v.success("Bulk tsunami wave equation verified")


def test_retarded_greens_function_solution(v):
    """Test retarded Green's function solution dimensional consistency."""
    # Paper: δρ(r₄,t) = M G^ret_(4)(r₄,t;v_L)
    # Green's function must have dimensions to make RHS match LHS

    lhs_dims = v.get_dim('delta_rho_4')  # M L^-4
    mass_dims = v.get_dim('M_sink')      # M

    # Therefore G^ret must have dimensions L^-4
    greens_dims = lhs_dims / mass_dims
    expected_greens = v.L**-4
    v.check_dims("Retarded Green's function G^ret_(4)", expected_greens, greens_dims)

    v.success("Retarded Green's function solution verified")




def test_the_tsunami_principle():
    """
    Main test function verifying "The Tsunami Principle" subsection.

    Tests the core mathematical physics documented in the subsection:
    - Bulk wave equation: ∂²_t δρ - v²_L ∇²₄ δρ = -M δ⁴(r₄) δ'(t)
    - Retarded Green's function solution: δρ = M G^ret_(4)(r₄,t;v_L)

    Total: 2 focused tests of documented equations
    """
    v = PhysicsVerificationHelper(
        "The Tsunami Principle",
        "Bulk density re-equilibration vs observable wave propagation",
        unit_system=UnitSystem.SI
    )

    v.section_header("Testing The Tsunami Principle")

    # Add core dimensions for bulk tsunami physics
    v.add_dimensions({
        'delta_rho_4': v.M * (v.L**-4),          # 4D density perturbation
        'v_L': v.L / v.T,                        # bulk longitudinal speed
        'M_sink': v.M,                           # sink mass
        'r_4': v.L,                              # 4D radial coordinate
        'delta4': v.L**-4,                       # 4D Dirac delta function
    }, allow_overwrite=True)

    # Core documented physics from "The Tsunami Principle" subsection
    v.info("\n--- Core Bulk Tsunami Physics (from documented subsection) ---")
    v.section("Bulk tsunami wave equation and retarded solutions")
    test_bulk_tsunami_wave_equation(v)
    test_retarded_greens_function_solution(v)

    # Final summary
    v.summary()


if __name__ == "__main__":
    test_the_tsunami_principle()
