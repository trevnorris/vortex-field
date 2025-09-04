# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains a mathematical framework modeling particles as quantized vortices in four-dimensional space. The framework attempts to derive particle masses, gravitational phenomena, and quantum mechanics from minimal geometric inputs, reducing the Standard Model's ~20 free parameters to just a few geometric constants.

## Development Commands

### LaTeX Documentation
- Compile main document: `cd doc && pdflatex main.tex`
- Multiple passes may be required for references: `cd doc && pdflatex main.tex && pdflatex main.tex`

### Python Environment
- Python 3.x required with: SymPy, NumPy, SciPy, Matplotlib
- No package manager configuration files present - dependencies managed manually
- Run individual Python scripts directly: `python script_name.py`

### Testing
- Individual derivation verification: run specific Python files in `derivations/` directories
- When running tests that use `derivations/helper.py` you should end the test with `--quiet` so that it only displays the necessary output.
- The tests in `derivations/` are to make sure the math is correct, so when writing new tests it's more important to have the test match what's in the paper rather than pass. So that we can find where there are actual issues with the derivations/dimensions/units.
- When updating tests never use terms like "updated" or "corrected". Tests should be timeless, and there's no need to mention that we've made changes to the test, in the test.

### Finding Subsections
- Within doc are multiple `.tex` files.
- To find the list of subsection, in order to find the range to read from when working with or reading from, a specific subsection run the command `grep -n '^\\subsection' -- doc/<filename>` where `<filename>` is the specific file you are looking at.

## Code Architecture

### Core Framework Structure

**`/derivations/`** - SymPy-based symbolic mathematics verification
- `helper.py` - Unified physics verification framework with dimensional analysis and unit system support
- `helper.test.py` - Comprehensive test suite for verification framework
- `mathematical_framework/` - Core theoretical derivations (projection mechanisms, conservation laws, foundational postulates)
- `emergent_particle_masses/` - Particle physics applications and mass calculations

**`/calculations/`** - NumPy/SciPy-based numerical predictions and validations
- `mathematical_framework/` - Vortex dynamics, stability analysis, parameter studies
- `emergent_particle_masses/` - Numerical mass predictions and experimental comparisons
- Individual photon deflection calculations at root level

**`/doc/`** - LaTeX documentation
- `main.tex` - Master document with formatting setup
- Individual `.tex` files for each major section (mathematical framework, particle masses, gravity, etc.)

### Key Components

**Physics Verification Framework (`derivations/helper.py`)**
- Multi-unit system support (SI, Gaussian, Heaviside-Lorentz)
- Dimensional analysis and smart symbol lookup with suggestions
- Verification patterns for wave equations, conservation laws, Poisson equations
- Section-based organization with batch checking capabilities
