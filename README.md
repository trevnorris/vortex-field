# Aether-Vortex Field Equations

A unified fluid model for gravity in flat space, reimagining the luminiferous aether as a compressible superfluid in four-dimensional space.

## Overview

This repository contains the theoretical framework, mathematical derivations, and numerical calculations for an alternative theory of gravity based on superfluid dynamics. Rather than curved spacetime, gravity emerges from fluid mechanical phenomena - particles as vortex structures draining aether into an extra dimension, creating pressure gradients and flows that manifest as gravitational attraction.

## Repository Structure

```
aether-vortex/
├── doc/
│   └── vortex_field_equations.tex    # Main theoretical document (LaTeX)
└── scripts/
    ├── derivations/                   # Symbolic mathematics verification
    │   ├── 4d_superfluid_framework_and_projections.py
    │   ├── derivation_of_the_scalar_field_equation.py
    │   ├── derivation_of_the_vector_field_equation.py
    │   ├── physical_postulates.py
    │   └── unified_equations_and_force_law.py
    └── calculations/                  # Numerical predictions & validations
        ├── 1d_draining_flow.py       # Basic drainage flow analysis
        ├── 4d_drainage.py            # 4D vortex drainage calculations
        ├── frame_dragging.py         # Frame-dragging effects (Gravity Probe B)
        ├── mercury.py                # Mercury perihelion precession
        └── particle_mass.py          # Vortex-based particle mass calculations
```

## Theory Summary

The model postulates that:
- Space is filled with a 4D compressible superfluid (the "aether")
- Particles are quantized vortex structures extending into the 4th dimension
- Gravity arises from aether drainage creating local density variations
- The framework reproduces General Relativity's predictions while remaining in flat Euclidean space

Key features:
- **Dual wave modes**: Longitudinal compression waves (potentially faster than light in 4D bulk) and transverse waves (light) at c
- **Physical intuition**: Effects stem from tangible fluid mechanics rather than abstract curved manifolds
- **Testable predictions**: Including lab-scale frame-dragging and chromatic shifts in black hole photon spheres

## Components

### Documentation (`/doc`)
- **vortex_field_equations.tex**: Complete theoretical exposition including:
  - Physical postulates and mathematical framework
  - Derivation of field equations from fluid dynamics
  - Post-Newtonian expansions matching GR predictions
  - Validation against astronomical observations

### Derivations (`/scripts/derivations`)
Python scripts using SymPy for symbolic verification of theoretical results:
- Mathematical consistency checks
- Derivation of scalar and vector field equations
- Verification of physical postulates and relationships
- Unified force law derivations

### Calculations (`/scripts/calculations`)
Python scripts using NumPy/SciPy for numerical predictions:
- Mercury's perihelion advance (43 arcseconds/century)
- Frame-dragging effects matching Gravity Probe B observations
- Vortex drainage flow patterns
- Particle mass emergence from vortex core volumes

## Requirements

- LaTeX distribution (for compiling the theoretical document)
- Python 3.x with:
  - SymPy (for symbolic mathematics)
  - NumPy/SciPy (for numerical calculations)
  - Matplotlib (for visualizations)

## Author

Written by Trevor Norris

## Citation

If you use this work in your research, please cite:
```
Norris, T. (2025). The Aether-Vortex Field Equations: A Unified Fluid Model 
for Gravity in Flat Space. [Repository/Publication details]
```
