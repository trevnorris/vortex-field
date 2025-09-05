# Physics Verification Helper API Documentation

## Overview

The `helper.py` module provides a comprehensive framework for dimensional analysis and mathematical verification of physics equations. It supports multiple unit systems (SI, Gaussian, Heaviside-Lorentz) and provides smart symbol lookup with fuzzy matching.

## Core Classes

### PhysicsVerificationHelper

Main verification class for physics tests.

#### Initialization
```python
PhysicsVerificationHelper(section_name, description="", unit_system=UnitSystem.SI, quiet=None)
```
- `section_name`: Name/number of the section being tested
- `description`: Optional detailed description
- `unit_system`: UnitSystem.SI, UnitSystem.GAUSSIAN, or UnitSystem.HEAVISIDE_LORENTZ
- `quiet`: If True, suppress output. If None, auto-detect from --quiet flag

#### Core Attributes
- `dims`: Dictionary of all available physical dimensions
- `prefactors`: Unit system-specific equation prefactors
- `results`: List of verification results
- `failed_checks`: List of failed checks
- `L, M, T, Q, Theta`: Base dimensional symbols

### UnitSystem (Enum)
- `UnitSystem.SI`
- `UnitSystem.GAUSSIAN`
- `UnitSystem.HEAVISIDE_LORENTZ`

### Custom Exceptions
- `UnknownSymbolError`: Raised for unknown symbols (includes suggestions)
- `DimensionalMismatchError`: Raised for dimensional inconsistencies

## Symbol and Dimension Management

### Getting Dimensions
```python
get_dim(name, suggest=True, strict=False)
```
Gets dimension for a symbol with fuzzy matching suggestions.

### Adding Custom Dimensions
```python
add_dimension(name, dimension, allow_overwrite=False)
add_dimensions(new_dims_dict, allow_overwrite=False)
declare_dimensionless(*names)
```

### Available Physical Quantities

#### Fundamental Constants
- `c`: Speed of light [L/T]
- `G`: Gravitational constant [L³/(M·T²)]
- `hbar`, `h`: Planck constants [M·L²/T]
- `e`: Elementary charge [Q]
- `epsilon_0`: Vacuum permittivity [M⁻¹·L⁻³·T²·Q²]
- `mu_0`: Vacuum permeability [M·L·Q⁻²]
- `k_B`: Boltzmann constant [M·L²·T⁻²·Θ⁻¹]

#### Electromagnetic Fields and Potentials
- `E`: Electric field [V/m]
- `B`: Magnetic field [V·s/m²]
- `D`: Electric displacement [Q/m²]
- `H`: Magnetic field intensity [Q·s⁻¹·m⁻¹]
- `A`: Vector potential [V·s/m]
- `Phi`: Electric potential [V]

#### 4-Vectors
- `J_mu`, `J0`, `J1`, `J2`, `J3`: EM 4-current components
- `A_mu`, `A0`, `A1`, `A2`, `A3`: EM 4-potential components

#### Gravitational/GEM Fields
- `Phi_g`: Gravitoelectric potential [L²/T²]
- `A_g`: Gravitomagnetic potential [L/T]
- `E_g`: Gravitoelectric field [L/T²]
- `B_g`: Gravitomagnetic field [T⁻¹]
- `g`: Gravitational acceleration [L/T²]

#### Densities and Currents
- `rho`: Mass density [M/L³]
- `rho_charge`: Charge density [Q/L³]
- `j_current`: Electric current density [Q·T⁻¹·L⁻²]
- `j_mass`: Mass current density [M/(L²·T)]

#### Quantum and Wave Quantities
- `psi`: Wavefunction [L⁻³/²]
- `omega`: Angular frequency [T⁻¹]
- `k`: Wavenumber [L⁻¹]
- `p`: Momentum [M·L/T]
- `lambda`: Wavelength [L]

#### Energy and Thermodynamics
- `E_energy`: Energy [M·L²/T²]
- `P`: Pressure [M/(L·T²)]
- `u_EM`: EM energy density [M·L⁻¹·T⁻²]

#### Coordinates and Geometry
- `x`, `y`, `z`, `r`: Spatial coordinates [L]
- `t`: Time [T]
- `theta`, `phi`: Angles [dimensionless]

## Verification Methods

### Dimensional Consistency
```python
check_dims(name, expr1, expr2, record=True, verbose=True, si_only=False)
```
Check if two expressions have the same dimensions.

### Equation Verification
```python
check_eq(name, lhs, rhs, record=True, verbose=True)
```
Check if LHS equals RHS after simplification.

### Zero Checks
```python
check_zero(name, expr, record=True, verbose=True)
```
Verify an expression is identically zero.

### Limit Calculations
```python
check_limit(name, expr, var, point, expected, record=True, verbose=True)
```
Check limit calculations.

### Asymptotic Analysis
```python
check_asymptotic(name, exact_expr, approx_expr, var, point=oo, record=True, verbose=True)
```
Verify asymptotic approximations.

## Dimensional Analysis Utilities

### Derivative Operations
```python
dt(dim)    # Time derivative: dim/T
dtt(dim)   # Second time derivative: dim/T²
dx(dim)    # Spatial derivative: dim/L
dxx(dim)   # Second spatial derivative: dim/L²
```

### Field Operations
```python
grad_dim(scalar_dim)   # Gradient: nabla * scalar_dim
div_dim(vector_dim)    # Divergence: div * vector_dim
curl_dim(vector_dim)   # Curl: curl * vector_dim
lap_dim(scalar_dim)    # Laplacian: laplacian * scalar_dim
```

### Safe Transcendental Functions
```python
exp_dimless(x, where="exp")
sin_dimless(x, where="sin")
cos_dimless(x, where="cos")
log_dimless(x, where="log")
assert_dimensionless(expr, where="")
validate_transcendentals(expr, where="expression")
```

## Output and Formatting

### Print Methods
```python
info(message)          # Respects quiet mode
success(message)       # Respects quiet mode
warning(message)       # Always prints
error(message)         # Always prints
debug(message)         # Only in debug mode
```

### Section Headers
```python
section(title, width=60)
subsection(title, width=50)
section_header(title)
```

### Results Summary
```python
summary(show_failed=True)  # Returns success rate percentage
```

## Specialized Verification Patterns

### Wave Equations
```python
verify_wave_equation(helper, name, time_term, space_term, source_term=None)
```
Verify ∂_tt φ - c²∇²φ = source

### Conservation Laws
```python
verify_conservation_law(helper, name, density_rate, flux_div, source=None)
```
Verify ∂_t ρ + ∇·J = source

### Poisson Equations
```python
verify_poisson_em(helper)           # EM Poisson equations
verify_poisson_grav(helper)         # Gravitational Poisson
verify_poisson_equation(helper, name, laplacian_term, source_term)
```

## Convenience Functions

### Quick Verification
```python
quick_verify(name, condition, details="", helper=None, expected_failure=False)
```
One-line verification without recording.

### Batch Operations
```python
batch_check_dims(helper, checks)      # List of (name, expr1, expr2) tuples
batch_check_eqs(helper, equations)    # List of (name, lhs, rhs) tuples
define_symbols_batch(names, **kwargs) # Define multiple symbols
```

### Analysis Tools
```python
print_dimensional_analysis(expression, helper, name="Expression")
create_section_verifier(section_name, description="")
```

## Usage Patterns

### Basic Test Structure
```python
from helper import PhysicsVerificationHelper

# Initialize
v = PhysicsVerificationHelper("Section 2.3", "Test Description")

# Add custom dimensions if needed
v.add_dimensions({
    'custom_field': v.M * v.L / v.T**2
})

# Run checks
v.check_dims("Field consistency", v.get_dim('E'), v.get_dim('V') / v.get_dim('x'))
v.check_eq("Maxwell equation", lhs_expr, rhs_expr)

# Generate summary
success_rate = v.summary()
```

### Quiet Mode Support
Tests automatically detect `--quiet` flag from command line. Use quiet mode for automated testing:
```bash
python test_file.py --quiet
```

### Unit System Support
```python
# Test in different unit systems
v_si = PhysicsVerificationHelper("Test", unit_system=UnitSystem.SI)
v_gauss = PhysicsVerificationHelper("Test", unit_system=UnitSystem.GAUSSIAN)
v_hl = PhysicsVerificationHelper("Test", unit_system=UnitSystem.HEAVISIDE_LORENTZ)
```

## Error Handling

The framework provides detailed error messages with suggestions:
- Unknown symbols trigger fuzzy matching suggestions
- Dimensional mismatches show expected vs actual dimensions
- Transcendental function arguments are validated for dimensionlessness

## Example Test Template

```python
from helper import PhysicsVerificationHelper, define_symbols_batch
import sympy as sp

v = PhysicsVerificationHelper("Section X.Y", "Description")

# Define symbols
t, x, y, z = define_symbols_batch(['t', 'x', 'y', 'z'], real=True)
E, B, rho = sp.symbols('E B rho')

# Verification
v.section("MAIN VERIFICATION")
v.check_dims("Gauss law",
             v.div_dim(v.get_dim('E')),
             v.get_dim('rho_charge') / v.get_dim('epsilon_0'))

success_rate = v.summary()
```