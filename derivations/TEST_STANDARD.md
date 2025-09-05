# Physics Verification Test Standard

This document defines the standardized structure for all physics verification tests in the derivations framework. Following this standard ensures compatibility with the automated test runner and consistent code organization.

## Standard Test Structure

### File Template

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
[Section Name] - Verification
====================================

Brief description of what this test module verifies, including the specific
mathematical relationships, equations, or physical principles being validated.

Based on doc/[document].tex, section [reference] (lines X-Y if specific).
"""

import os
import sys
import sympy as sp
from sympy import symbols, pi, simplify  # Import specific symbols as needed

# Add parent directory to path to import helper
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper import (
    PhysicsVerificationHelper,
    UnitSystem,  # If needed for non-SI tests
    verify_wave_equation,  # Import specific verification patterns as needed
    quick_verify,
    # ... other imports as needed
)


# Helper functions (if needed for complex tests)
def test_specific_aspect_1(v):
    """
    Test a specific aspect of the physics being verified.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Specific Aspect 1")
    
    # Test implementation here
    # All tests use the shared 'v' instance
    v.check_dims("description", expr1, expr2)
    # ... more tests ...
    
    v.success("Specific aspect 1 verified")


def test_specific_aspect_2(v):
    """
    Test another specific aspect.
    
    Args:
        v: PhysicsVerificationHelper instance
    """
    v.subsection("Specific Aspect 2")
    
    # Test implementation here
    # ... tests ...
    
    v.success("Specific aspect 2 verified")


def test_[module_name]():
    """
    Main test function for [Section Name].
    
    This function coordinates all verification tests for the section,
    calling helper functions as needed and providing a single entry point.
    
    Returns:
        float: Success rate (0-100) from verification summary
    """
    # Initialize verification helper
    v = PhysicsVerificationHelper(
        "Section Name",
        "Brief description of what's being verified"
    )
    
    v.section("MAIN SECTION TITLE")
    
    # Define any custom symbols if needed
    # (Use define_symbols_batch for multiple symbols)
    # x, y, z = define_symbols_batch(['x', 'y', 'z'], real=True)
    
    # Add custom dimensions if needed
    # v.add_dimensions({
    #     'custom_quantity': v.M * v.L / v.T**2,
    # })
    
    # Call test functions in logical order
    v.info("\n--- 1) First Test Group ---")
    test_specific_aspect_1(v)
    
    v.info("\n--- 2) Second Test Group ---") 
    test_specific_aspect_2(v)
    
    # Add more test calls as needed...
    
    # Return success rate for test runner integration
    return v.summary()


if __name__ == "__main__":
    success_rate = test_[module_name]()
    # Exit with non-zero code if tests failed (for CI/automation)
    if success_rate < 100.0:
        sys.exit(1)
```

## Key Requirements

### 1. Function Naming Convention
- **Main function**: `test_[module_name]()` where `[module_name]` exactly matches the Python filename (without `.py`)
- **Helper functions**: `test_[descriptive_name](v)` taking the helper instance as argument

### 2. Return Value
- The main function **must** return `v.summary()` 
- This returns the success rate (0-100) for integration with the test runner

### 3. Helper Function Pattern
For complex tests with multiple aspects:

```python
def test_wave_equations(v):
    """Test wave equation dimensional consistency."""
    v.subsection("Wave Equations")
    # ... test implementation using 'v' ...
    v.success("Wave equations verified")

def test_constants(v): 
    """Test electromagnetic constant relationships."""
    v.subsection("EM Constants")
    # ... test implementation using 'v' ...
    v.success("Constants verified")

def test_fixing_constants_and_waves():
    """Main test coordinating all sub-tests."""
    v = PhysicsVerificationHelper(
        "Fixing the Constants and Waves",
        "EM wave equations and constant relationships"
    )
    
    v.section("FIXING THE CONSTANTS AND WAVES VERIFICATION")
    
    # Call helper functions in sequence
    v.info("\n--- 1) Wave Equations ---")
    test_wave_equations(v)
    
    v.info("\n--- 2) EM Constants ---") 
    test_constants(v)
    
    return v.summary()
```

### 4. Import Structure
- Always use the standard path manipulation for helper import
- Import specific symbols from SymPy rather than `import sympy as sp` when possible
- Group imports logically (standard library, sympy, helper)

### 5. Documentation
- Module docstring with clear description and document reference
- Function docstrings for main function and complex helper functions
- Inline comments for non-obvious physics or mathematics

### 6. Error Handling
- Main block should exit with code 1 if success rate < 100%
- This enables CI/automation to detect test failures

## Examples by Complexity

### Simple Test (Direct Implementation)
```python
def test_simple_relationship():
    """Test a straightforward dimensional relationship."""
    v = PhysicsVerificationHelper("Simple Test", "Basic dimensional check")
    
    v.section("SIMPLE RELATIONSHIP TEST")
    
    # Direct test implementation
    lhs = v.get_dim('E') * v.get_dim('B')
    rhs = v.get_dim('S_poynting') * v.get_dim('mu_0')
    v.check_dims("Poynting vector relationship", lhs, rhs)
    
    return v.summary()
```

### Complex Test (Multiple Helper Functions)
```python
def test_conservation_laws_detailed():
    """Test comprehensive conservation law framework."""
    v = PhysicsVerificationHelper("Conservation Laws", "Mass, momentum, energy")
    
    v.section("CONSERVATION LAWS VERIFICATION")
    
    # Multiple organized test groups
    test_mass_continuity(v)
    test_momentum_balance(v) 
    test_energy_conservation(v)
    test_noether_currents(v)
    
    return v.summary()
```

## Migration Guide

### Converting Existing Tests

1. **Pattern 1 (Direct execution)**: Wrap existing code in main function
2. **Pattern 2 (Function with return)**: Usually already compliant, just rename function
3. **Pattern 3 (Function without return)**: Add `return` to `v.summary()` call

### Checklist for Conversion
- [ ] Main function named `test_[module_name]()`
- [ ] Function returns `v.summary()`
- [ ] Helper functions take `v` parameter
- [ ] Main block handles exit codes
- [ ] Imports follow standard pattern
- [ ] Documentation includes section reference

## Test Runner Integration

The standardized structure enables a simple test runner:

```python
def run_test_module(module_path):
    module = import_module(module_path)
    test_function = getattr(module, f"test_{module_name}")
    success_rate = test_function()
    return success_rate
```

This eliminates the need for complex pattern detection and fallback mechanisms.

## Benefits

1. **Consistency**: All tests follow the same structure
2. **Maintainability**: Easy to understand and modify tests  
3. **Automation**: Simple integration with test runners and CI
4. **Organization**: Clear separation between test logic and coordination
5. **Reliability**: No guesswork about how to run tests

Following this standard ensures that all physics verification tests integrate seamlessly with the automated testing infrastructure while maintaining clear, readable code that accurately validates the mathematical framework.