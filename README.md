# Aether-Vortex Field Equations

A geometric framework modeling particles as quantized vortices in four-dimensional space, yielding precise predictions for particle masses, gravitational phenomena, and quantum mechanics from minimal geometric inputs.

## Overview

This repository contains a mathematical framework that attempts to address three fundamental mysteries in physics: the origin of particle masses, the weakness of gravity, and quark confinement. By modeling particles as topological defects (quantized vortices) in a four-dimensional medium, the framework derives lepton masses accurate to fractions of a percent, reproduces general relativity, and naturally explains confinement.

**Key Achievement**: The model reduces the Standard Model's ~20 free parameters to just a few geometric inputs while making specific, testable predictions that could falsify the theory.

## Theory Summary

### Core Features
- **Particles as vortices**: Leptons are quantized vortex configurations in 4D space with n=1,2,3 windings
- **Golden-ratio scaling**: Mass hierarchy emerges from self-similarity: adding helical layers then rescaling by r → 1+1/r
- **Baryon confinement**: Baryons are single tri-phase vortex loops—"quarks" are inseparable wave phases, not particles
- **Emergent gravity**: Gravitational effects arise from vortex-induced flows in flat 4D Euclidean space

### Key Predictions
**Already confirmed:**
- Lepton masses: electron (exact), muon (-0.18% error), tau (+0.10% error)
- Mercury perihelion: 43.0"/century (observed: 42.98 ± 0.04)
- Binary pulsar decay matches observation to parts in 10¹²

**Testable predictions:**
- **4-lepton anomaly**: Excess production near √s = 33 GeV without narrow resonance
- **Threefold baryon structure**: F(q) ~ F₀(q) + F₃(q)cos(3φ) in nucleon form factors
- **No fourth stable lepton**: Would-be 16.48 GeV lepton exceeds stability threshold
- **Intrinsic decoherence**: Γ(d) = Γ₀ + γ₂d² scaling with path separation

## Components

### Documentation (`/doc`)
- **mathematical_framework.tex**: Complete theoretical exposition and mathematical derivations
- **introduction.tex**: Author's motivation, problem overview, and key predictions

### Derivations (`/derivations`)
SymPy-based symbolic mathematics verification:
- **mathematical_framework/**: Core theoretical derivations, field equations, projection mechanisms
- **emergent_particle_masses/**: Particle physics applications, mass calculations, photon theory

### Calculations (`/calculations`) 
NumPy/SciPy-based numerical predictions and validations:
- **mathematical_framework/**: Vortex dynamics, stability analysis, parameter studies
- **emergent_particle_masses/**: Particle mass predictions and experimental comparisons

## What Makes This Different

This framework inverts the typical approach to unification:
- **Parameter reduction**: Dozens of Standard Model parameters → few geometric inputs
- **Testable predictions**: Specific energies and signatures that can falsify the theory  
- **Mathematical simplicity**: Undergraduate fluid mechanics rather than advanced geometry
- **Conceptual clarity**: Mysterious phenomena like confinement become topological necessities

**Philosophical stance**: This is primarily a tool for discovering mathematical patterns. Whether nature actually employs 4D vortices matters less than the precise correspondences the framework reveals and the constraints it places on any successful theory.

## Requirements

- LaTeX distribution (for compiling documentation)
- Python 3.x with SymPy, NumPy/SciPy, Matplotlib

## Testing

- Run all tests: `python derivations/run_tests.py`
- Run all tests (quiet output): `python derivations/run_tests.py --quiet`

## Author's Note

*"I am not a physicist. I am a computer programmer who set out to test AI capabilities with a weekend experiment. This paper was never supposed to exist... After weeks of testing increasingly complex phenomena, the model continued delivering precise results. I present this work not as a claim to have found 'the answer,' but as a discovery of remarkable mathematical patterns that demand explanation."*

—From the introduction

## Status

**Current focus**: Completing derivations and calculations. The framework successfully explains lepton masses, gravitational phenomena, and baryon confinement, with several testable predictions awaiting experimental verification.
