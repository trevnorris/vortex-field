# Demonstration Scripts Specification for 4D Vortex Framework

## Overview

This document specifies demonstration scripts to complement the mathematical verification already completed in `mathematical_framework.py`. While the SymPy verification proves all ~100 mathematical relationships in the framework, these scripts will provide numerical demonstrations, visualizations, and pedagogical tools to make the abstract mathematics concrete and accessible.

### Context
The framework models particles as topological defects (vortices) in a 4D compressible medium, projecting to 3D dynamics. Key features include:
- 4-fold circulation enhancement from geometric projection
- Emergent gravitational-like dynamics
- Golden ratio stability in vortex configurations
- Resolution of preferred frame paradox through Machian principles

---

## Script Specifications

### 1. **hemisphere_integral_numerical.py**

**Purpose**: Numerically demonstrate the critical integral ∫₀^∞ dw'/(ρ²+w'²)^(3/2) = 1/ρ² that underlies the 4-fold enhancement.

**Mathematical Background**:
- This integral represents the contribution from hemispherical vortex sheet projections
- It's key to understanding why each hemisphere contributes exactly Γ to the total circulation

**Required Functionality**:
```python
def hemisphere_integral(rho, w_max=None, n_points=1000):
    """
    Numerically compute ∫₀^w_max dw'/(ρ²+w'²)^(3/2)

    Parameters:
    - rho: radial distance in the 3D plane
    - w_max: upper integration limit (None for adaptive)
    - n_points: number of integration points

    Returns:
    - integral value
    - estimated error
    - convergence data for plotting
    """

def demonstrate_convergence():
    """
    Show convergence to 1/ρ² for multiple ρ values
    Generate plots showing:
    - Integral value vs upper limit
    - Error vs upper limit (log scale)
    - Comparison across different ρ values
    """
```

**Expected Outputs**:
- Convergence plots showing approach to analytical value
- Table of numerical vs analytical values for ρ ∈ [0.1, 1, 10]
- Error analysis showing exponential convergence

---

### 2. **greens_function_visualization.py**

**Purpose**: Visualize the 4D Green's function and its 3D projection, demonstrating causality preservation.

**Mathematical Background**:
- 4D wave equation: ∂²φ/∂t² - v_L²∇₄²φ = S(r₄,t)
- Green's function: G₄(t,r₄) with support potentially outside light cone for v_L > c
- Projection: G_proj(t,r) = ∫dw G₄(t,√(r²+w²))

**Required Functionality**:
```python
def compute_greens_4d(t, r, w, v_L, xi):
    """
    Compute 4D Green's function at given spacetime point
    Include both prompt (delta function) and tail contributions
    """

def project_to_3d(G_4d_func, t, r, w_range):
    """
    Integrate 4D Green's function over w to get 3D projection
    Show how bulk modes at v_L > c project to observable modes at c
    """

def visualize_lightcone():
    """
    Create 2D plots (t vs r) showing:
    - 4D Green's function support (extends beyond r = ct if v_L > c)
    - 3D projected Green's function (confined to r ≤ ct)
    - Smearing effects from finite ξ
    """
```

**Expected Outputs**:
- Heatmap of G₄(t,r,w=0) showing bulk propagation
- Heatmap of G_proj(t,r) showing observable confinement to c
- Animation showing how projection enforces causality
- Quantitative demonstration of smearing Δt ~ ξ²/(2rv_L)

---

### 3. **vortex_flow_field.py**

**Purpose**: Visualize the complete 4D flow field around a vortex sheet and its 3D projections.

**Mathematical Background**:
- Vortex sheet in 4D: 2D surface with circulation Γ
- Four contributions: direct, upper/lower hemispheres, induced w-flow
- Each contributes Γ for total 4Γ

**Required Functionality**:
```python
def velocity_field_4d(x, y, z, w, vortex_params):
    """
    Compute 4D velocity field components (v_x, v_y, v_z, v_w)
    Include all contributions:
    - Azimuthal flow from sheet
    - Drainage velocity v_w
    - Induced components
    """

def project_velocity_field(v_4d, w_slice=0):
    """
    Project 4D velocity to 3D slice at given w
    Show how 4 separate circulation contributions sum
    """

def create_flow_visualization():
    """
    Generate:
    - Streamline plots in (x,y) plane at w=0
    - Vector field plots showing all 4 contributions
    - 3D visualization of drainage into w-direction
    - Animation rotating around vortex
    """
```

**Expected Outputs**:
- 2D streamline plot showing 4Γ total circulation
- Separated plots for each contribution
- 3D visualization with color-coded flow magnitude
- Quantitative verification: line integrals = Γ for each component

---

### 4. **golden_ratio_stability.py**

**Purpose**: Demonstrate why the golden ratio φ emerges as the unique stable configuration for hierarchical vortex structures.

**Mathematical Background**:
- Energy functional: E ∝ (x-1)²/2 - ln(x) where x = R_{n+1}/R_n
- Minimization yields x² = x + 1, giving φ = (1+√5)/2
- Rational ratios lead to resonant instabilities

**Required Functionality**:
```python
def vortex_interaction_energy(radius_ratio, n_levels=3):
    """
    Compute total energy for hierarchical vortex configuration
    Include:
    - Elastic energy from deformation
    - Topological constraints
    - Resonance penalties for rational ratios
    """

def simulate_dynamics(initial_ratios, time_steps):
    """
    Evolve vortex configuration under energy minimization
    Show:
    - Rational ratios develop instabilities
    - System converges to golden ratio
    - Resonance catastrophes for p/q ratios
    """

def visualize_energy_landscape():
    """
    Create plots showing:
    - Energy vs radius ratio
    - Stability analysis (second derivative)
    - Phase space trajectories
    - Fibonacci spiral emergence
    """
```

**Expected Outputs**:
- Energy landscape plot with minimum at φ
- Time evolution showing convergence to golden ratio
- Comparison of stable (φ) vs unstable (rational) configurations
- Animation of resonance catastrophe for rational ratios

---

### 5. **wave_propagation_demo.py**

**Purpose**: Numerically solve the scalar wave equation with position-dependent v_eff(r) to show gravitational-like effects.

**Mathematical Background**:
- Wave equation: (1/v_eff²)∂²Ψ/∂t² - ∇²Ψ = 4πGρ
- Near mass: v_eff ≈ c(1 - GM/(2c²r))
- Predicts wave slowing, deflection, redshift

**Required Functionality**:
```python
def setup_velocity_profile(mass_distribution):
    """
    Create v_eff(r) profile based on mass distribution
    Include:
    - Point masses
    - Extended bodies
    - Background cosmological density
    """

def solve_wave_equation(initial_conditions, v_eff_profile, t_max):
    """
    Numerical PDE solver for wave equation
    Use finite difference or spectral methods
    Handle variable coefficients from v_eff(r)
    """

def analyze_propagation():
    """
    Extract and visualize:
    - Wave packet trajectories
    - Time delays from v_eff < c
    - Deflection angles (gravitational lensing analog)
    - Frequency shifts (redshift analog)
    """
```

**Expected Outputs**:
- Animation of wave packet propagating past mass
- Quantitative comparison with GR predictions
- Shapiro delay demonstration
- Lensing angle vs impact parameter plot

---

### 6. **reconnection_dynamics.py**

**Purpose**: Simulate vortex reconnection events and the associated energy barriers and drainage bursts.

**Mathematical Background**:
- Energy barrier: ΔE ≈ ρ₄D⁰Γ²ξ²ln(L/ξ)/(4π)
- Reconnection releases trapped flux into bulk
- Acts as "valve" for drainage regulation

**Required Functionality**:
```python
def compute_reconnection_energy(vortex_config, separation):
    """
    Calculate energy barrier for vortex reconnection
    Include:
    - Core overlap energy
    - Topological constraint relaxation
    - Logarithmic factors from long-range interactions
    """

def simulate_reconnection_event(initial_state, parameters):
    """
    Time-evolve through reconnection
    Track:
    - Core positions and shapes
    - Energy vs time
    - Flux release into w-direction
    - Sound wave emission
    """

def visualize_reconnection():
    """
    Create visualizations showing:
    - Before/during/after vortex configurations
    - Energy barrier crossing
    - Drainage burst in w-direction
    - Radiated wave patterns
    """
```

**Expected Outputs**:
- Energy vs separation showing barrier
- Time series of reconnection event
- Quantification of released flux
- Parameter study of barrier height vs Γ, ξ, L

---

### 7. **machian_frame_demo.py**

**Purpose**: Illustrate how distributed vortex sinks eliminate any global rest frame, resolving the preferred frame problem.

**Mathematical Background**:
- Background potential: Ψ_bg ∝ r² (outward acceleration)
- Cosmic inflows: Ψ_cosmic ∝ -r² (inward from distant matter)
- Balance yields local inertial frames, no global rest

**Required Functionality**:
```python
def compute_flow_field(sink_distribution):
    """
    Calculate net flow at each point from all sinks
    Include:
    - Local vortices (discrete sinks)
    - Cosmological background (continuous)
    - Boundary conditions at cosmic scale
    """

def find_balance_surfaces():
    """
    Identify surfaces where net flow vanishes
    These define local inertial frames
    Show they form complex, position-dependent pattern
    """

def demonstrate_no_global_frame():
    """
    Visualize:
    - Flow fields from different sink distributions
    - How balance surfaces shift with distribution
    - Absence of any universal rest configuration
    - Connection to Mach's principle
    """
```

**Expected Outputs**:
- Vector field plots showing complex flow patterns
- Highlighted balance surfaces (local rest frames)
- Animation showing how frames shift with observer motion
- Quantitative measure of frame-dragging effects

---

### 8. **parameter_robustness.py**

**Purpose**: Comprehensive parameter study showing which results are geometric/topological (parameter-independent) vs which scale with physical parameters.

**Mathematical Background**:
- Geometric: 4-fold factor, golden ratio
- Scaling: Energy barriers, timescales, propagation speeds
- Mixed: Gravitational effects (geometric structure, physical scales)

**Required Functionality**:
```python
def parameter_sweep(calculation_func, param_ranges, fixed_params):
    """
    Systematic variation of parameters
    Test ranges spanning multiple orders of magnitude
    Identify which outputs change and how
    """

def analyze_sensitivities():
    """
    For each output quantity:
    - Determine functional dependence on parameters
    - Classify as geometric invariant or scaling quantity
    - Compare with theoretical predictions
    """

def generate_robustness_report():
    """
    Create comprehensive report showing:
    - Which results are truly parameter-independent
    - Scaling laws for parameter-dependent quantities
    - Confidence bounds on key predictions
    """
```

**Expected Outputs**:
- Matrix plot of sensitivities
- Log-log plots showing scaling relationships
- Table of invariant vs scaling quantities
- Executive summary of robustness findings

---

## Implementation Guidelines

### General Requirements
1. **Documentation**: Each script should include comprehensive docstrings explaining the physics and mathematics
2. **Visualization**: Use matplotlib/plotly for publication-quality figures with proper labels and units
3. **Performance**: Optimize numerical computations for reasonable run times (< 5 minutes per script)
4. **Validation**: Include self-consistency checks and comparison with analytical results where available
5. **Accessibility**: Provide command-line interfaces with sensible defaults for easy exploration

### Output Standards
- All plots should be saved in high resolution (300+ DPI)
- Numerical results should include error estimates
- Animations should be exportable as MP4/GIF
- Each script should generate a summary report of key findings

### Dependencies
- NumPy for numerical computation
- SciPy for integration and PDE solving
- Matplotlib/Plotly for visualization
- SymPy for symbolic verification (where needed)
- Optional: Numba/Cython for performance-critical sections

### Testing
Each script should include:
- Unit tests for core functions
- Integration tests for full workflows
- Regression tests against known results
- Parameter validation and error handling

---

## Priority Order

1. **4_fold_enhancement.py** - Already complete, serves as template
2. **hemisphere_integral_numerical.py** - Core mathematical demonstration
3. **vortex_flow_field.py** - Visual proof of 4-fold mechanism
4. **wave_propagation_demo.py** - Shows gravitational correspondence
5. **golden_ratio_stability.py** - Demonstrates emergent geometry
6. **greens_function_visualization.py** - Resolves causality concerns
7. **machian_frame_demo.py** - Addresses conceptual issues
8. **reconnection_dynamics.py** - Shows dynamical processes
9. **parameter_robustness.py** - Validates entire framework

This prioritization moves from core mathematical demonstrations to more complex dynamical simulations, building understanding progressively.
