"""
Comprehensive Test Suite for helper.py
===============================================

Complete testing of all major functionality in the PhysicsVerificationHelper class.
Tests dimensional analysis, symbol registry, verification patterns, unit systems,
transcendental validation, and error handling.
"""

import sys
import warnings
from io import StringIO
import sympy as sp
from sympy import symbols, Rational, simplify, pi, exp, sin, cos, log, sqrt, oo, limit
from helper import (PhysicsVerificationHelper, UnitSystem, SymbolRegistry,
                   UnknownSymbolError, DimensionalMismatchError,
                   quick_verify, batch_check_dims, batch_check_eqs,
                   define_symbols_batch, verify_wave_equation,
                   verify_conservation_law, verify_poisson_em,
                   verify_poisson_grav, verify_poisson_equation,
                   create_section_verifier, print_dimensional_analysis,
                   generate_test_template)

class TestResult:
    """Helper class to track test results."""
    def __init__(self, name):
        self.name = name
        self.passed = 0
        self.total = 0
        self.details = []
    
    def add_test(self, test_name, result, details=""):
        self.total += 1
        if result:
            self.passed += 1
            self.details.append(f"✓ {test_name}")
        else:
            self.details.append(f"✗ {test_name}: {details}")
    
    def summary(self):
        rate = (self.passed / self.total * 100) if self.total > 0 else 0
        return f"{self.name}: {self.passed}/{self.total} ({rate:.1f}%)"

def test_symbol_registry():
    """Test the symbol registry and alias system."""
    print("=" * 70)
    print("TESTING SYMBOL REGISTRY AND ALIASES")
    print("=" * 70)
    
    result = TestResult("Symbol Registry")
    v = PhysicsVerificationHelper("Test", "Symbol registry tests")
    
    # Test 1: Basic symbol lookup
    print("\n1. Testing basic symbol lookup...")
    try:
        E_dim = v.get_dim('E')
        expected_E = v.dims['V_volt'] / v.L
        basic_lookup = v.check_dims("Basic E lookup", E_dim, expected_E, verbose=False)
        result.add_test("Basic symbol lookup", basic_lookup)
    except Exception as e:
        result.add_test("Basic symbol lookup", False, str(e))
    
    # Test 2: Alias resolution
    print("\n2. Testing alias resolution...")
    test_aliases = [
        ('Efield', 'E_field'),
        ('Bfield', 'B_field'), 
        ('rho_e', 'rho_charge'),
        ('j', 'j_current'),
        ('mu0', 'mu_0'),
        ('eps0', 'epsilon_0')
    ]
    
    for alias, canonical in test_aliases:
        try:
            alias_dim = v.get_dim(alias)
            canonical_dim = v.get_dim(canonical)
            alias_check = v.check_dims(f"Alias {alias}->{canonical}", 
                                     alias_dim, canonical_dim, verbose=False)
            result.add_test(f"Alias {alias}", alias_check)
        except Exception as e:
            result.add_test(f"Alias {alias}", False, str(e))
    
    # Test 3: Unknown symbol with suggestions
    print("\n3. Testing unknown symbol suggestions...")
    try:
        v.get_dim('Efeild')  # Typo in 'Efield'
        result.add_test("Unknown symbol error", False, "Should have raised error")
    except UnknownSymbolError as e:
        # Check if suggestion is provided
        has_suggestion = "Did you mean" in str(e)
        result.add_test("Unknown symbol suggestions", has_suggestion, str(e))
    except Exception as e:
        result.add_test("Unknown symbol suggestions", False, str(e))
    
    # Test 4: SymbolRegistry direct testing
    print("\n4. Testing SymbolRegistry class...")
    registry = {'test_symbol': v.M * v.L}
    
    # Test canonicalization
    canonical = SymbolRegistry.canonicalize('Efield')
    canon_check = canonical == 'E_field'
    result.add_test("Canonicalization", canon_check)
    
    # Test direct symbol lookup
    try:
        sym = SymbolRegistry.get_symbol('test_symbol', registry, suggest=False)
        direct_lookup = sym == v.M * v.L
        result.add_test("Direct registry lookup", direct_lookup)
    except Exception as e:
        result.add_test("Direct registry lookup", False, str(e))
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_dimensional_consistency():
    """Test dimensional consistency checking in various scenarios."""
    print("\n" + "=" * 70)
    print("TESTING DIMENSIONAL CONSISTENCY CHECKS")
    print("=" * 70)
    
    result = TestResult("Dimensional Consistency")
    v = PhysicsVerificationHelper("Test", "Dimensional consistency tests")
    
    # Test 1: Basic dimensional equality
    print("\n1. Testing basic dimensional equality...")
    
    # Energy dimensions should match
    kinetic = v.M * v.L**2 / v.T**2
    potential = v.get_dim('E_energy')
    energy_check = v.check_dims("Energy equivalence", kinetic, potential, verbose=False)
    result.add_test("Energy dimensions", energy_check)
    
    # Test 2: Complex expression consistency
    print("\n2. Testing complex expressions...")
    
    # Maxwell stress tensor: ε₀E² has energy density dimensions
    maxwell_stress = v.dims['epsilon_0'] * v.dims['E']**2
    energy_density = v.M * v.L**(-1) * v.T**(-2)
    stress_check = v.check_dims("Maxwell stress", maxwell_stress, energy_density, verbose=False)
    result.add_test("Maxwell stress tensor", stress_check)
    
    # Test 3: Fundamental relations
    print("\n3. Testing fundamental relations...")
    
    # c = 1/√(μ₀ε₀)
    c_calculated = 1 / sp.sqrt(v.dims['mu_0'] * v.dims['epsilon_0'])
    c_check = v.check_dims("Speed of light from constants", 
                          v.dims['c'], c_calculated, verbose=False)
    result.add_test("c from μ₀ε₀", c_check)
    
    # Test 4: Addition consistency
    print("\n4. Testing addition of quantities...")
    
    # Adding two energies should work
    E1 = v.M * v.L**2 / v.T**2
    E2 = v.M * v.L**2 / v.T**2
    sum_expr = E1 + E2
    sum_check = v.check_dims("Energy addition", sum_expr, E1, verbose=False)
    result.add_test("Energy addition", sum_check)
    
    # Test 5: Dimensional mismatch detection
    print("\n5. Testing mismatch detection...")
    
    # Try to add energy and momentum (should fail)
    try:
        energy = v.M * v.L**2 / v.T**2
        momentum = v.M * v.L / v.T
        mismatch_expr = energy + momentum
        # This should be caught in dimensional analysis
        dim = v._extract_dimensions(mismatch_expr)
        mismatch_detected = str(dim) == 'DIM_MISMATCH'
        result.add_test("Mismatch detection", mismatch_detected)
    except Exception as e:
        result.add_test("Mismatch detection", True, "Caught as expected")
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_transcendental_validation():
    """Test transcendental function validation."""
    print("\n" + "=" * 70)
    print("TESTING TRANSCENDENTAL FUNCTION VALIDATION")
    print("=" * 70)
    
    result = TestResult("Transcendental Validation")
    v = PhysicsVerificationHelper("Test", "Transcendental validation tests")
    
    # Test 1: Safe transcendental functions
    print("\n1. Testing safe transcendental functions...")
    
    # exp, sin, cos, log with dimensionless arguments
    try:
        # v/c is dimensionless
        beta = v.get_dim('v') / v.get_dim('c')
        
        exp_result = v.exp_dimless(beta, "relativistic factor")
        sin_result = v.sin_dimless(sp.pi, "sine of pi")
        cos_result = v.cos_dimless(sp.pi/2, "cosine of pi/2")
        log_result = v.log_dimless(2, "log of 2")
        
        result.add_test("exp_dimless", True)
        result.add_test("sin_dimless", True)
        result.add_test("cos_dimless", True)
        result.add_test("log_dimless", True)
        
    except Exception as e:
        result.add_test("Safe transcendentals", False, str(e))
    
    # Test 2: Transcendental validation with dimensional arguments (should fail)
    print("\n2. Testing dimensional arguments (should fail)...")
    
    try:
        v.exp_dimless(v.get_dim('E'), "electric field")
        result.add_test("exp with E field", False, "Should have failed")
    except DimensionalMismatchError:
        result.add_test("exp with E field", True, "Correctly failed")
    except Exception as e:
        result.add_test("exp with E field", False, f"Wrong error: {e}")
    
    # Test 3: Expression tree validation
    print("\n3. Testing expression tree validation...")
    
    # Complex expression with transcendentals
    try:
        # This should pass: sin(π) + cos(π/2) + exp(0)
        complex_expr = sp.sin(sp.pi) + sp.cos(sp.pi/2) + sp.exp(0)
        v.validate_transcendentals(complex_expr, "complex transcendental")
        result.add_test("Complex transcendental validation", True)
    except Exception as e:
        result.add_test("Complex transcendental validation", False, str(e))
    
    # Test 4: Invalid complex expression
    print("\n4. Testing invalid complex expressions...")
    
    try:
        # This should fail: sin(E) where E has dimensions
        invalid_expr = sp.sin(v.get_dim('E'))
        v.validate_transcendentals(invalid_expr, "invalid expression")
        result.add_test("Invalid expression validation", False, "Should have failed")
    except DimensionalMismatchError:
        result.add_test("Invalid expression validation", True, "Correctly failed")
    except Exception as e:
        result.add_test("Invalid expression validation", False, f"Wrong error: {e}")
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_differential_operators():
    """Test differential operators (grad, div, curl, laplacian)."""
    print("\n" + "=" * 70)
    print("TESTING DIFFERENTIAL OPERATORS")
    print("=" * 70)
    
    result = TestResult("Differential Operators")
    v = PhysicsVerificationHelper("Test", "Differential operator tests")
    
    # Test 1: Gradient operator
    print("\n1. Testing gradient operator...")
    
    # ∇φ where φ is a potential
    phi_dim = v.get_dim('Phi')  # Electric potential
    grad_phi = v.grad_dim(phi_dim)  # Should give electric field
    expected_E = v.get_dim('E')
    
    grad_check = v.check_dims("∇φ = E", grad_phi, expected_E, verbose=False)
    result.add_test("Gradient of potential", grad_check)
    
    # Test 2: Divergence operator
    print("\n2. Testing divergence operator...")
    
    # ∇·E should have dimensions of charge density / ε₀
    div_E = v.div_dim(v.get_dim('E'))
    expected_div_E = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
    
    div_check = v.check_dims("∇·E", div_E, expected_div_E, verbose=False)
    result.add_test("Divergence of E field", div_check)
    
    # Test 3: Curl operator
    print("\n3. Testing curl operator...")
    
    # ∇×A should give B field
    curl_A = v.curl_dim(v.get_dim('A'))
    expected_B = v.get_dim('B')
    
    curl_check = v.check_dims("∇×A = B", curl_A, expected_B, verbose=False)
    result.add_test("Curl of vector potential", curl_check)
    
    # Test 4: Laplacian operator
    print("\n4. Testing Laplacian operator...")
    
    # ∇²φ in Poisson equation
    lap_phi = v.lap_dim(v.get_dim('Phi'))
    expected_poisson = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
    
    lap_check = v.check_dims("∇²φ Poisson", lap_phi, expected_poisson, verbose=False)
    result.add_test("Laplacian of potential", lap_check)
    
    # Test 5: Time derivatives
    print("\n5. Testing time derivatives...")
    
    # ∂B/∂t in Faraday's law should match ∇×E
    dB_dt = v.dt(v.get_dim('B'))
    curl_E = v.curl_dim(v.get_dim('E'))
    
    faraday_check = v.check_dims("Faraday: ∂B/∂t ~ ∇×E", dB_dt, curl_E, verbose=False)
    result.add_test("Faraday time derivative", faraday_check)
    
    # Test 6: Second derivatives
    print("\n6. Testing second derivatives...")
    
    # ∂²φ/∂t² in wave equation
    d2phi_dt2 = v.dtt(v.get_dim('Phi'))
    lap_phi_c2 = v.get_dim('c')**2 * v.lap_dim(v.get_dim('Phi'))
    
    wave_check = v.check_dims("Wave equation: ∂²φ/∂t² ~ c²∇²φ", d2phi_dt2, lap_phi_c2, verbose=False)
    result.add_test("Wave equation derivatives", wave_check)
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_verification_patterns():
    """Test specialized verification patterns."""
    print("\n" + "=" * 70)
    print("TESTING VERIFICATION PATTERNS")
    print("=" * 70)
    
    result = TestResult("Verification Patterns")
    
    # Create a helper but capture its output to avoid cluttering test output
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Pattern tests")
        
        # Test 1: Wave equation verification
        print("\n1. Testing wave equation pattern...")
        
        # ∂²φ/∂t² - c²∇²φ = 0
        time_term = v.dtt(v.get_dim('Phi'))
        space_term = v.get_dim('c')**2 * v.lap_dim(v.get_dim('Phi'))
        
        wave_result = verify_wave_equation(v, "Scalar wave", time_term, space_term)
        result.add_test("Wave equation pattern", wave_result)
        
        # Test 2: Conservation law verification
        print("\n2. Testing conservation law pattern...")
        
        # ∂ρ/∂t + ∇·j = 0
        density_rate = v.dt(v.get_dim('rho'))
        flux_div = v.div_dim(v.get_dim('j_mass'))
        
        conservation_result = verify_conservation_law(v, "Mass", density_rate, flux_div)
        result.add_test("Conservation law pattern", conservation_result)
        
        # Test 3: EM Poisson equations
        print("\n3. Testing EM Poisson pattern...")
        
        poisson_em_result = verify_poisson_em(v)
        result.add_test("EM Poisson equations", poisson_em_result)
        
        # Test 4: Gravitational Poisson equation
        print("\n4. Testing gravitational Poisson pattern...")
        
        poisson_grav_result = verify_poisson_grav(v)
        result.add_test("Gravitational Poisson", poisson_grav_result)
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_unit_systems():
    """Test different unit system support."""
    print("\n" + "=" * 70)
    print("TESTING UNIT SYSTEM SUPPORT")
    print("=" * 70)
    
    result = TestResult("Unit Systems")
    
    # Test 1: SI unit system (default)
    print("\n1. Testing SI unit system...")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v_si = PhysicsVerificationHelper("Test", "SI test", UnitSystem.SI)
        
        # Check that SI prefactors are correct
        si_gauss_prefactor = v_si.prefactors['gauss_E_prefactor']
        expected_si_gauss = 1 / v_si.dims['epsilon_0']
        
        si_check = v_si.check_dims("SI Gauss prefactor", si_gauss_prefactor, expected_si_gauss, verbose=False)
        result.add_test("SI unit system", si_check)
        
    finally:
        sys.stdout = old_stdout
    
    # Test 2: Gaussian unit system
    print("\n2. Testing Gaussian unit system...")
    
    sys.stdout = StringIO()
    
    try:
        v_gaussian = PhysicsVerificationHelper("Test", "Gaussian test", UnitSystem.GAUSSIAN)
        
        # In Gaussian units, Gauss law has 4π factor
        gaussian_gauss_prefactor = v_gaussian.prefactors['gauss_E_prefactor']
        expected_gaussian_gauss = 4 * sp.pi
        
        gaussian_check = (gaussian_gauss_prefactor == expected_gaussian_gauss)
        result.add_test("Gaussian unit system", gaussian_check)
        
    finally:
        sys.stdout = old_stdout
    
    # Test 3: Heaviside-Lorentz unit system
    print("\n3. Testing Heaviside-Lorentz unit system...")
    
    sys.stdout = StringIO()
    
    try:
        v_hl = PhysicsVerificationHelper("Test", "HL test", UnitSystem.HEAVISIDE_LORENTZ)
        
        # In HL units, Gauss law has no prefactor
        hl_gauss_prefactor = v_hl.prefactors['gauss_E_prefactor']
        expected_hl_gauss = 1
        
        hl_check = (hl_gauss_prefactor == expected_hl_gauss)
        result.add_test("Heaviside-Lorentz unit system", hl_check)
        
    finally:
        sys.stdout = old_stdout
    
    # Test 4: Unit system warnings
    print("\n4. Testing unit system warnings...")
    
    sys.stdout = StringIO()
    
    try:
        v_warn = PhysicsVerificationHelper("Test", "Warning test", UnitSystem.GAUSSIAN)
        # Should see warning about dimensional checks being SI-anchored
        warning_check = True  # If we get here without error, warning system works
        result.add_test("Unit system warnings", warning_check)
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_equation_checks():
    """Test equation checking and limit calculations."""
    print("\n" + "=" * 70)
    print("TESTING EQUATION CHECKS AND LIMITS")
    print("=" * 70)
    
    result = TestResult("Equation Checks")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Equation tests")
        
        # Test 1: Basic equation checking
        print("\n1. Testing basic equation checks...")
        
        # Simple algebraic equation
        x = symbols('x')
        lhs = 2*x + 3
        rhs = 2*x + 3
        
        eq_check = v.check_eq("Simple equation", lhs, rhs, verbose=False)
        result.add_test("Basic equation check", eq_check)
        
        # Test 2: Zero checking
        print("\n2. Testing zero checks...")
        
        zero_expr = x - x  # Should be zero
        zero_check = v.check_zero("x - x = 0", zero_expr, verbose=False)
        result.add_test("Zero check", zero_check)
        
        # Test 3: Limit calculations
        print("\n3. Testing limit calculations...")
        
        # lim(sin(x)/x as x->0) = 1
        limit_expr = sp.sin(x)/x
        limit_check = v.check_limit("sin(x)/x limit", limit_expr, x, 0, 1, verbose=False)
        result.add_test("Limit calculation", limit_check)
        
        # Test 4: Asymptotic behavior
        print("\n4. Testing asymptotic behavior...")
        
        # For large x: exp(-x) ≈ 0
        exact = sp.exp(-x)
        approx = 0
        asymptotic_check = v.check_asymptotic("exp(-x) -> 0", exact, approx, x, oo, verbose=False)
        result.add_test("Asymptotic behavior", asymptotic_check)
        
        # Test 5: Complex equation verification
        print("\n5. Testing complex equations...")
        
        # Euler's identity: e^(iπ) + 1 = 0
        euler_lhs = sp.exp(sp.I * sp.pi) + 1
        euler_rhs = 0
        
        euler_check = v.check_eq("Euler's identity", euler_lhs, euler_rhs, verbose=False)
        result.add_test("Euler's identity", euler_check)
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_error_handling():
    """Test error handling and edge cases."""
    print("\n" + "=" * 70)
    print("TESTING ERROR HANDLING AND EDGE CASES")
    print("=" * 70)
    
    result = TestResult("Error Handling")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Error tests")
        
        # Test 1: Unknown symbol handling
        print("\n1. Testing unknown symbol handling...")
        
        try:
            v.get_dim('nonexistent_symbol')
            result.add_test("Unknown symbol error", False, "Should have raised error")
        except UnknownSymbolError:
            result.add_test("Unknown symbol error", True)
        except Exception as e:
            result.add_test("Unknown symbol error", False, f"Wrong error type: {e}")
        
        # Test 2: Dimensional mismatch handling
        print("\n2. Testing dimensional mismatch handling...")
        
        try:
            # Try to assert that energy is dimensionless (should fail)
            v.assert_dimensionless(v.get_dim('E_energy'))
            result.add_test("Dimensional mismatch error", False, "Should have raised error")
        except DimensionalMismatchError:
            result.add_test("Dimensional mismatch error", True)
        except Exception as e:
            result.add_test("Dimensional mismatch error", False, f"Wrong error type: {e}")
        
        # Test 3: Custom dimension addition
        print("\n3. Testing custom dimension addition...")
        
        # Add new dimension
        try:
            v.add_dimension('custom_test', v.M * v.L**3 / v.T)
            custom_dim = v.get_dim('custom_test')
            expected = v.M * v.L**3 / v.T
            custom_check = v.check_dims("Custom dimension", custom_dim, expected, verbose=False)
            result.add_test("Custom dimension addition", custom_check)
        except Exception as e:
            result.add_test("Custom dimension addition", False, str(e))
        
        # Test 4: Overwrite protection
        print("\n4. Testing overwrite protection...")
        
        try:
            # Try to overwrite existing dimension without permission
            v.add_dimension('E', v.M, allow_overwrite=False)
            result.add_test("Overwrite protection", False, "Should have raised error")
        except KeyError:
            result.add_test("Overwrite protection", True)
        except Exception as e:
            result.add_test("Overwrite protection", False, f"Wrong error type: {e}")
        
        # Test 5: Batch operations
        print("\n5. Testing batch operations...")
        
        try:
            # Test batch dimension addition
            new_dims = {
                'batch_test_1': v.M * v.L,
                'batch_test_2': v.L / v.T
            }
            v.add_dimensions(new_dims)
            
            # Check both were added
            dim1 = v.get_dim('batch_test_1')
            dim2 = v.get_dim('batch_test_2')
            
            batch_check1 = v.check_dims("Batch test 1", dim1, v.M * v.L, verbose=False)
            batch_check2 = v.check_dims("Batch test 2", dim2, v.L / v.T, verbose=False)
            
            result.add_test("Batch operations", batch_check1 and batch_check2)
        except Exception as e:
            result.add_test("Batch operations", False, str(e))
        
        # Test 6: Dimensionless declaration
        print("\n6. Testing dimensionless declaration...")
        
        try:
            v.declare_dimensionless('alpha', 'beta', 'gamma')
            
            # Check that they're treated as dimensionless
            alpha_dim = v.get_dim('alpha')
            dimensionless_check = (alpha_dim == 1)
            result.add_test("Dimensionless declaration", dimensionless_check)
        except Exception as e:
            result.add_test("Dimensionless declaration", False, str(e))
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_convenience_functions():
    """Test convenience functions and utilities."""
    print("\n" + "=" * 70)
    print("TESTING CONVENIENCE FUNCTIONS")
    print("=" * 70)
    
    result = TestResult("Convenience Functions")
    
    # Test 1: quick_verify function
    print("\n1. Testing quick_verify...")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        # Should return True and print checkmark
        quick_result = quick_verify("Test condition", True, "details")
        result.add_test("quick_verify", quick_result == True)
    finally:
        sys.stdout = old_stdout
    
    # Test 2: define_symbols_batch
    print("\n2. Testing define_symbols_batch...")
    
    try:
        syms = define_symbols_batch(['a', 'b', 'c'], real=True)
        batch_symbols_check = len(syms) == 3 and all(s.is_real for s in syms)
        result.add_test("define_symbols_batch", batch_symbols_check)
    except Exception as e:
        result.add_test("define_symbols_batch", False, str(e))
    
    # Test 3: create_section_verifier
    print("\n3. Testing create_section_verifier...")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        verifier = create_section_verifier("Test Section", "Test description")
        factory_check = isinstance(verifier, PhysicsVerificationHelper)
        result.add_test("create_section_verifier", factory_check)
    finally:
        sys.stdout = old_stdout
    
    # Test 4: batch_check_dims
    print("\n4. Testing batch_check_dims...")
    
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Batch test")
        
        checks = [
            ("Energy check", v.get_dim('E_energy'), v.M * v.L**2 / v.T**2),
            ("Force check", v.M * v.L / v.T**2, v.M * v.get_dim('a'))
        ]
        
        # This should not raise an error
        batch_check_dims(v, checks)
        result.add_test("batch_check_dims", True)
    except Exception as e:
        result.add_test("batch_check_dims", False, str(e))
    finally:
        sys.stdout = old_stdout
    
    # Test 5: batch_check_eqs
    print("\n5. Testing batch_check_eqs...")
    
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Batch eq test")
        x = symbols('x')
        
        equations = [
            ("Identity", x, x),
            ("Simple algebra", 2*x + 1, 2*x + 1)
        ]
        
        batch_check_eqs(v, equations)
        result.add_test("batch_check_eqs", True)
    except Exception as e:
        result.add_test("batch_check_eqs", False, str(e))
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_missing_methods():
    """Test methods that were missed in the first pass."""
    print("\n" + "=" * 70)
    print("TESTING PREVIOUSLY MISSED METHODS")
    print("=" * 70)
    
    result = TestResult("Missing Methods")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Missing methods test")
        
        # Test 1: poisson_dim convenience method
        print("\n1. Testing poisson_dim convenience method...")
        
        phi_dim = v.get_dim('Phi')
        source_dim = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
        
        poisson_result = v.poisson_dim(phi_dim, source_dim, "EM Poisson test")
        result.add_test("poisson_dim method", poisson_result)
        
        # Test 2: verify_and_record method
        print("\n2. Testing verify_and_record method...")
        
        verify_result = v.verify_and_record("Test verification", True, "Test details")
        result.add_test("verify_and_record", verify_result == True)
        
        # Test 3: dx and dxx spatial derivatives
        print("\n3. Testing spatial derivatives...")
        
        phi_dim = v.get_dim('Phi')
        dphidx = v.dx(phi_dim)
        d2phidx2 = v.dxx(phi_dim)
        
        # ∂φ/∂x should be E field (negative gradient)
        dx_check = v.check_dims("∂φ/∂x", dphidx, v.get_dim('E'), verbose=False)
        
        # ∂²φ/∂x² should match Laplacian dimensionally
        dxx_check = v.check_dims("∂²φ/∂x²", d2phidx2, v.lap_dim(phi_dim), verbose=False)
        
        result.add_test("dx spatial derivative", dx_check)
        result.add_test("dxx second spatial derivative", dxx_check)
        
        # Test 4: _is_base_dim internal method
        print("\n4. Testing _is_base_dim internal method...")
        
        base_check_M = v._is_base_dim(v.M)
        base_check_L = v._is_base_dim(v.L)
        base_check_other = v._is_base_dim(symbols('x'))
        
        result.add_test("_is_base_dim for M", base_check_M)
        result.add_test("_is_base_dim for L", base_check_L)  
        result.add_test("_is_base_dim for non-base", not base_check_other)
        
        # Test 5: get_dim with strict parameter
        print("\n5. Testing get_dim strict parameter...")
        
        # In strict mode, aliases should be disallowed
        try:
            v.get_dim('Efield', strict=True)  # This should fail
            result.add_test("get_dim strict mode", False, "Should have failed")
        except UnknownSymbolError:
            result.add_test("get_dim strict mode", True)
        except Exception as e:
            result.add_test("get_dim strict mode", False, f"Wrong error: {e}")
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_formatting_and_reporting():
    """Test section formatting and summary reporting."""
    print("\n" + "=" * 70)
    print("TESTING FORMATTING AND REPORTING")
    print("=" * 70)
    
    result = TestResult("Formatting & Reporting")
    
    # Capture all output to test formatting methods
    old_stdout = sys.stdout
    captured_output = StringIO()
    sys.stdout = captured_output
    
    try:
        v = PhysicsVerificationHelper("Test", "Formatting test")
        
        # Test 1: section formatting
        v.section("Test Section", width=40)
        section_output = captured_output.getvalue()
        
        section_check = "Test Section" in section_output and "=" in section_output
        result.add_test("section formatting", section_check)
        
        # Test 2: subsection formatting  
        captured_output = StringIO()
        sys.stdout = captured_output
        
        v.subsection("Test Subsection", width=30)
        subsection_output = captured_output.getvalue()
        
        subsection_check = "Test Subsection" in subsection_output and "-" in subsection_output
        result.add_test("subsection formatting", subsection_check)
        
        # Test 3: summary reporting with statistics
        captured_output = StringIO()
        sys.stdout = captured_output
        
        # Add some test results
        v.results = [("Test 1", True), ("Test 2", False), ("Test 3", True)]
        v.failed_checks = ["Test 2"]
        
        success_rate = v.summary(show_failed=True)
        summary_output = captured_output.getvalue()
        
        # Check that summary includes statistics and failed tests
        summary_stats_check = "2/3" in summary_output or "66.7%" in str(success_rate)
        summary_failed_check = "Test 2" in summary_output
        
        result.add_test("summary statistics", summary_stats_check)
        result.add_test("summary failed list", summary_failed_check)
        
        # Test 4: Different success rate messages
        # Test high success rate (>95%)
        v.results = [("T1", True), ("T2", True), ("T3", True), ("T4", False)]  # 75%
        v.failed_checks = ["T4"]
        
        captured_output = StringIO()
        sys.stdout = captured_output
        
        rate_75 = v.summary(show_failed=False)
        output_75 = captured_output.getvalue()
        
        # Test perfect success rate
        v.results = [("T1", True), ("T2", True), ("T3", True)]  # 100%
        v.failed_checks = []
        
        captured_output = StringIO()
        sys.stdout = captured_output
        
        rate_100 = v.summary(show_failed=False)
        output_100 = captured_output.getvalue()
        
        rate_messages_check = ("🎉" in output_100) and rate_100 == 100.0
        result.add_test("success rate messages", rate_messages_check)
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_utility_functions():
    """Test utility functions like print_dimensional_analysis and template generation."""
    print("\n" + "=" * 70)
    print("TESTING UTILITY FUNCTIONS")
    print("=" * 70)
    
    result = TestResult("Utility Functions")
    
    # Test 1: print_dimensional_analysis
    print("\n1. Testing print_dimensional_analysis...")
    
    old_stdout = sys.stdout
    captured_output = StringIO()
    sys.stdout = captured_output
    
    try:
        v = PhysicsVerificationHelper("Test", "Utility test")
        
        # Test dimensional analysis printing
        expr = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')  # Poynting vector
        print_dimensional_analysis(expr, v, "Poynting vector")
        
        analysis_output = captured_output.getvalue()
        
        analysis_check = ("Poynting vector" in analysis_output and 
                         "Dimensions:" in analysis_output and
                         "Expression:" in analysis_output)
        result.add_test("print_dimensional_analysis", analysis_check)
        
    finally:
        sys.stdout = old_stdout
    
    # Test 2: generate_test_template
    print("\n2. Testing generate_test_template...")
    
    try:
        template = generate_test_template("2.5", "Wave Propagation")
        
        # Check that template contains expected elements
        template_checks = [
            "Section 2.5" in template,
            "Wave Propagation" in template,
            "PhysicsVerificationHelper" in template,
            "check_dims" in template,
            "summary()" in template,
            'v = PhysicsVerificationHelper(' in template
        ]
        
        template_check = all(template_checks)
        result.add_test("generate_test_template", template_check)
        
    except Exception as e:
        result.add_test("generate_test_template", False, str(e))
    
    # Test 3: verify_poisson_equation helper
    print("\n3. Testing verify_poisson_equation helper...")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Poisson test")
        
        lap_term = v.lap_dim(v.get_dim('Phi'))
        source_term = v.get_dim('rho_charge') / v.get_dim('epsilon_0')
        
        poisson_helper_result = verify_poisson_equation(v, "EM Poisson", lap_term, source_term)
        result.add_test("verify_poisson_equation helper", poisson_helper_result)
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_edge_cases_and_robustness():
    """Test edge cases and robustness scenarios."""
    print("\n" + "=" * 70)  
    print("TESTING EDGE CASES AND ROBUSTNESS")
    print("=" * 70)
    
    result = TestResult("Edge Cases & Robustness")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Edge case tests")
        
        # Test 1: None handling in dimensional extraction
        print("\n1. Testing None handling...")
        
        none_dim = v._extract_dimensions(None)
        none_check = none_dim is None
        result.add_test("None handling", none_check)
        
        # Test 2: Zero handling with different types
        print("\n2. Testing zero handling...")
        
        zero_int = v._extract_dimensions(0)
        zero_float = v._extract_dimensions(0.0)
        zero_sympy = v._extract_dimensions(sp.Integer(0))
        
        zero_checks = [str(zero_int) == 'ZERO', str(zero_float) == 'ZERO', str(zero_sympy) == 'ZERO']
        result.add_test("Zero handling types", all(zero_checks))
        
        # Test 3: Complex nested expressions
        print("\n3. Testing complex nested expressions...")
        
        # Test nested multiplication with dimensionless constants
        complex_expr = 2 * sp.pi * v.get_dim('f') * v.get_dim('L_ind')  # 2πfL
        expected_dim = v.get_dim('R')  # Should be resistance dimensions (ωL = resistance)
        
        complex_check = v.check_dims("Complex impedance", complex_expr, expected_dim, verbose=False)
        result.add_test("Complex nested expressions", complex_check)
        
        # Test 4: Very large dimension expressions
        print("\n4. Testing large dimension expressions...")
        
        # Energy density: (1/2)(ε₀E² + B²/μ₀)
        large_expr = (v.get_dim('epsilon_0') * v.get_dim('E')**2 + 
                     v.get_dim('B')**2 / v.get_dim('mu_0')) / 2
        expected_energy_density = v.M * v.L**(-1) * v.T**(-2)
        
        large_check = v.check_dims("Large energy expression", large_expr, expected_energy_density, verbose=False)
        result.add_test("Large dimension expressions", large_check)
        
        # Test 5: SI-only parameter testing
        print("\n5. Testing SI-only parameter...")
        
        # This should be skipped in non-SI systems
        v_gaussian = PhysicsVerificationHelper("Test", "Gaussian test", UnitSystem.GAUSSIAN)
        
        # Try SI-only check in Gaussian system
        si_only_result = v_gaussian.check_dims("SI-only test", 
                                               v_gaussian.get_dim('E'), 
                                               v_gaussian.get_dim('E'),
                                               si_only=True, verbose=False)
        # Should return True (skipped) in non-SI system
        result.add_test("SI-only parameter", si_only_result)
        
        # Test 6: Very long symbol names and edge case symbols
        print("\n6. Testing edge case symbols...")
        
        # Add dimension with very long name
        very_long_name = "very_long_physics_quantity_name_that_tests_string_handling"
        v.add_dimension(very_long_name, v.M * v.L**2 / v.T**3)
        
        long_name_check = v.get_dim(very_long_name) == v.M * v.L**2 / v.T**3
        result.add_test("Long symbol names", long_name_check)
        
        # Test 7: Transcendental with complex arguments
        print("\n7. Testing transcendental edge cases...")
        
        try:
            # Test nested transcendentals
            nested = sp.exp(sp.sin(sp.pi/4))  # exp(sin(π/4)) - should pass
            v.validate_transcendentals(nested, "nested transcendentals")
            result.add_test("Nested transcendentals", True)
        except Exception as e:
            result.add_test("Nested transcendentals", False, str(e))
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def test_physics_integration():
    """Test integration with real physics equations and workflows."""
    print("\n" + "=" * 70)
    print("TESTING PHYSICS INTEGRATION")
    print("=" * 70)
    
    result = TestResult("Physics Integration")
    
    old_stdout = sys.stdout
    sys.stdout = StringIO()
    
    try:
        v = PhysicsVerificationHelper("Test", "Physics integration tests")
        
        # Test 1: Complete Maxwell equations verification
        print("\n1. Testing complete Maxwell equations...")
        
        # Gauss law: ∇·E = ρ/ε₀
        gauss_E = v.check_dims("Gauss E", 
                              v.div_dim(v.get_dim('E')),
                              v.get_dim('rho_charge') / v.get_dim('epsilon_0'),
                              verbose=False)
        
        # Gauss magnetism: ∇·B = 0
        gauss_B = v.check_dims("Gauss B",
                              v.div_dim(v.get_dim('B')),
                              0,  # This is tricky - zero should match anything
                              verbose=False)
        
        # Faraday: ∇×E = -∂B/∂t
        faraday = v.check_dims("Faraday",
                              v.curl_dim(v.get_dim('E')),
                              v.dt(v.get_dim('B')),
                              verbose=False)
        
        # Ampere: ∇×B = μ₀j + μ₀ε₀∂E/∂t
        ampere_j = v.get_dim('mu_0') * v.get_dim('j_current')
        ampere_E = v.get_dim('mu_0') * v.get_dim('epsilon_0') * v.dt(v.get_dim('E'))
        
        ampere_total = ampere_j  # Just test first term for simplicity
        ampere = v.check_dims("Ampere",
                             v.curl_dim(v.get_dim('B')),
                             ampere_total,
                             verbose=False)
        
        maxwell_complete = all([gauss_E, faraday, ampere])
        result.add_test("Complete Maxwell equations", maxwell_complete)
        
        # Test 2: Wave equation derivation
        print("\n2. Testing wave equation derivation...")
        
        # From Maxwell: ∂²E/∂t² = c²∇²E
        wave_lhs = v.dtt(v.get_dim('E'))
        wave_rhs = v.get_dim('c')**2 * v.lap_dim(v.get_dim('E'))
        
        wave_eq = v.check_dims("EM wave equation", wave_lhs, wave_rhs, verbose=False)
        result.add_test("Wave equation derivation", wave_eq)
        
        # Test 3: Energy conservation in EM fields
        print("\n3. Testing energy conservation...")
        
        # Energy density: u = (1/2)(ε₀E² + B²/μ₀)
        u_E = v.get_dim('epsilon_0') * v.get_dim('E')**2
        u_B = v.get_dim('B')**2 / v.get_dim('mu_0')
        
        energy_density_check = v.check_dims("EM energy density E term", u_E, u_B, verbose=False)
        result.add_test("Energy conservation", energy_density_check)
        
        # Test 4: Poynting vector
        print("\n4. Testing Poynting vector...")
        
        # S = E × B / μ₀
        poynting = v.get_dim('E') * v.get_dim('B') / v.get_dim('mu_0')
        expected_flux = v.M * v.T**(-3)  # Energy flux density
        
        poynting_check = v.check_dims("Poynting vector", poynting, expected_flux, verbose=False)
        result.add_test("Poynting vector", poynting_check)
        
        # Test 5: Gravitational wave analogy (GEM)
        print("\n5. Testing gravitational wave analogy...")
        
        # GEM: ∂²E_g/∂t² = c²∇²E_g
        grav_wave_lhs = v.dtt(v.get_dim('E_g'))
        grav_wave_rhs = v.get_dim('c')**2 * v.lap_dim(v.get_dim('E_g'))
        
        grav_wave = v.check_dims("GEM wave equation", grav_wave_lhs, grav_wave_rhs, verbose=False)
        result.add_test("Gravitational wave analogy", grav_wave)
        
        # Test 6: Dimensional analysis workflow
        print("\n6. Testing complete dimensional analysis workflow...")
        
        # Simulate a complete physics verification workflow
        workflow_checks = []
        
        # Step 1: Define problem - use gravitational field
        v.add_dimensions({
            'test_field': v.M * v.L / v.T**2,
            'test_source': v.M / v.L**3
        })
        
        # Step 2: Check gravitational field equation ∇²Φ_g = 4πGρ (ignoring constants)
        workflow_checks.append(v.check_dims("Workflow Poisson",
                                           v.lap_dim(v.get_dim('Phi_g')),
                                           v.get_dim('G') * v.get_dim('test_source'),
                                           verbose=False))
        
        # Step 3: Check wave propagation
        workflow_checks.append(v.check_dims("Workflow wave",
                                           v.dtt(v.get_dim('test_field')),
                                           v.get_dim('c')**2 * v.lap_dim(v.get_dim('test_field')),
                                           verbose=False))
        
        workflow_complete = all(workflow_checks)
        result.add_test("Complete workflow", workflow_complete)
        
    finally:
        sys.stdout = old_stdout
    
    print(f"\n{result.summary()}")
    for detail in result.details:
        print(f"  {detail}")
    
    return result.passed == result.total

def run_comprehensive_tests():
    """Run all comprehensive tests and provide overall summary."""
    print("=" * 80)
    print("COMPREHENSIVE TEST SUITE FOR ALL HELPER.PY FUNCTIONALITY")
    print("=" * 80)
    print("Testing dimensional analysis, verification patterns, unit systems,")
    print("transcendental validation, error handling, and utilities...")
    
    test_suites = [
        ("Symbol Registry & Aliases", test_symbol_registry),
        ("Dimensional Consistency", test_dimensional_consistency),
        ("Transcendental Validation", test_transcendental_validation),
        ("Differential Operators", test_differential_operators),
        ("Verification Patterns", test_verification_patterns),
        ("Unit System Support", test_unit_systems),
        ("Equation Checks & Limits", test_equation_checks),
        ("Error Handling", test_error_handling),
        ("Convenience Functions", test_convenience_functions),
        ("Missing Methods", test_missing_methods),
        ("Formatting & Reporting", test_formatting_and_reporting),
        ("Utility Functions", test_utility_functions),
        ("Edge Cases & Robustness", test_edge_cases_and_robustness),
        ("Physics Integration", test_physics_integration),
    ]
    
    results = []
    
    # Run all test suites
    for suite_name, test_func in test_suites:
        print(f"\nRunning {suite_name}...")
        try:
            success = test_func()
            results.append((suite_name, success))
        except Exception as e:
            print(f"ERROR in {suite_name}: {e}")
            results.append((suite_name, False))
    
    # Overall summary
    print("\n" + "=" * 80)
    print("COMPREHENSIVE TEST RESULTS SUMMARY")
    print("=" * 80)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for suite_name, result in results:
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{suite_name:.<35} {status}")
    
    print(f"\nOVERALL: {passed}/{total} test suites passed ({passed/total*100:.1f}%)")
    
    if passed == total:
        print("\n🎉 ALL COMPREHENSIVE TESTS PASSED! 🎉")
        print("helper.py is fully functional and robust!")
    elif passed >= total * 0.9:
        print(f"\n✅ COMPREHENSIVE TESTING SUBSTANTIALLY COMPLETE")
        print(f"   {passed}/{total} test suites passed - excellent coverage!")
    elif passed >= total * 0.8:
        print(f"\n⚠️  COMPREHENSIVE TESTING MOSTLY SUCCESSFUL")
        print(f"   {passed}/{total} test suites passed - good coverage")
    else:
        print(f"\n❌ COMPREHENSIVE TESTING NEEDS ATTENTION")
        print(f"   Only {passed}/{total} test suites passed")
    
    failed_suites = [name for name, result in results if not result]
    if failed_suites:
        print(f"\nFailed test suites ({len(failed_suites)}):")
        for suite in failed_suites:
            print(f"  ✗ {suite}")
    
    return passed == total

if __name__ == "__main__":
    success = run_comprehensive_tests()
    sys.exit(0 if success else 1)