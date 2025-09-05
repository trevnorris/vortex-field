# Quantum Test Remediation Plan

## Problem Analysis

The current quantum tests have a **critical flaw**: they verify dimensional consistency of individual terms but **fail to test the actual mathematical relationships and equations** from the paper. This results in meaningless 100% pass rates that don't validate the physics.

### Specific Issues Found

1. **Missing Equation Testing**: Tests check `v.check_dims("term", expr, expected_dim)` but not `v.check_eq("equation", lhs, rhs)`
2. **No Mathematical Relationships**: Tests don't verify that terms actually combine correctly in equations
3. **False Confidence**: 100% pass rates don't reflect mathematical correctness
4. **Dimensional-Only Focus**: Tests miss coefficient errors, sign errors, and structural mistakes

### Example from Madelung Test

**Current approach (wrong)**:
```python
v.check_dims("∂_t ρ", drho_dt_dim, v.M * v.L**(-3) * v.T**(-1))
v.check_dims("∇·(ρ velocity field)", div_current_dim, v.M * v.L**(-3) * v.T**(-1))
```

**Correct approach (missing)**:
```python
# Test actual continuity equation: ∂_t ρ + ∇·(ρv) = 0
lhs = v.dt(rho_density)
rhs = -v.div_dim(current_density)
v.check_eq("Continuity equation", lhs, rhs)
```

## Remediation Strategy

### Phase 1: Test Analysis and Categorization

For each of the 11 quantum test files, systematically identify:

1. **Equations from paper** that need verification
2. **Current dimensional checks** that are insufficient
3. **Missing mathematical relationships** that should be tested
4. **Derivation steps** that need verification

### Phase 2: Enhanced Test Implementation

Replace or supplement dimensional checks with:

1. **Equation Verification**: Use `v.check_eq()` to test actual mathematical relationships
2. **Derivation Testing**: Step-by-step verification of mathematical derivations
3. **Coefficient Validation**: Ensure numerical factors and signs are correct
4. **Consistency Checking**: Verify equations are self-consistent and properly related
5. **Limit Testing**: Check classical limits, scaling behavior, and asymptotic forms

### Phase 3: Agent-Based Remediation

Deploy specialized agents to fix each test file with:

1. **Detailed equation analysis** from the source subsection
2. **Mathematical relationship identification**
3. **Enhanced test implementation** with proper equation checking
4. **Validation against the paper** to ensure correctness

## Implementation Plan

### Step 1: Create Enhanced Test Templates

Develop patterns for proper equation testing:

```python
# Pattern 1: Test complete equations
def test_continuity_equation(v):
    """Test ∂_t ρ + ∇·(ρ(∇S - qA)/m_*) = 0"""

    # Define terms
    drho_dt = v.dt(v.get_dim('rho'))
    velocity = (v.grad_dim(v.get_dim('S_phase')) - v.get_dim('q') * v.get_dim('A')) / v.get_dim('m_star')
    div_current = v.div_dim(v.get_dim('rho') * velocity)

    # Test the actual equation
    v.check_eq("Continuity equation", drho_dt, -div_current)

# Pattern 2: Test derivation steps
def test_polar_decomposition_derivation(v):
    """Test that ψ = √ρ e^(iS/ℏ) leads to correct real/imag separation"""

    # Define polar form
    psi_polar = sqrt(rho) * exp(I * S / hbar_eff)

    # Test Schrödinger equation with polar form
    schrodinger_lhs = I * hbar_eff * diff(psi_polar, t)
    schrodinger_rhs = -hbar_eff**2/(2*m_star) * laplacian(psi_polar) + V * psi_polar

    # Verify equation holds symbolically
    v.check_eq("Schrödinger with polar form", schrodinger_lhs, schrodinger_rhs)

    # Test separation into real/imaginary parts yields expected equations
    eq_real = re(schrodinger_lhs - schrodinger_rhs)
    eq_imag = im(schrodinger_lhs - schrodinger_rhs)

    v.check_eq("Real part yields Hamilton-Jacobi", eq_real, expected_HJ_equation)
    v.check_eq("Imaginary part yields continuity", eq_imag, expected_continuity_equation)

# Pattern 3: Test consistency relationships
def test_quantum_potential_consistency(v):
    """Test that Q[ρ] has correct form and satisfies required relationships"""

    # Define quantum potential
    Q_rho = -hbar_eff**2/(2*m_star) * laplacian(sqrt(rho))/sqrt(rho)

    # Test dimensional consistency
    v.check_dims("Quantum potential", Q_rho, v.get_dim('E_energy')/v.get_dim('rho'))

    # Test classical limit: Q[ρ] → 0 as ℏ_eff → 0
    Q_classical_limit = limit(Q_rho, hbar_eff, 0)
    v.check_eq("Classical limit", Q_classical_limit, 0)

    # Test WKB limit: Q[ρ] ≪ kinetic energy for slowly varying ρ
    # (This would require more sophisticated analysis)
```

### Step 2: Agent Deployment Strategy

Deploy agents sequentially with enhanced instructions:

```
AGENT INSTRUCTION TEMPLATE:
=========================

CRITICAL MISSION: Fix the fundamental testing flaw in [subsection_name].py

PROBLEM: Current test only checks dimensional consistency. It does NOT verify the actual mathematical equations and relationships from the paper.

YOUR TASK:
1. READ the subsection [subsection_name] from doc/quantum.tex carefully
2. IDENTIFY every equation, derivation, and mathematical relationship
3. REWRITE the test to verify ACTUAL EQUATIONS using v.check_eq(), not just dimensions
4. TEST mathematical derivation steps, not just final dimensional consistency
5. ENSURE the test would catch coefficient errors, sign errors, and structural mistakes

EXAMPLES OF WHAT TO CHANGE:

❌ WRONG (current approach):
v.check_dims("∂_t ρ term", dt_rho, v.M * v.L**(-3) * v.T**(-1))
v.check_dims("∇·J term", div_J, v.M * v.L**(-3) * v.T**(-1))

✅ CORRECT (what you should do):
lhs = v.dt(rho)
rhs = -v.div_dim(current_density)
v.check_eq("Continuity equation: ∂_t ρ + ∇·J = 0", lhs + rhs, 0)

REQUIREMENT: The new test should be able to detect if there are actual mathematical errors in the paper's equations, not just dimensional inconsistencies.
```

### Step 3: Validation and Quality Control

For each remediated test:

1. **Run the enhanced test** and verify it tests actual equations
2. **Introduce deliberate errors** to confirm the test catches them
3. **Compare against paper equations** to ensure accuracy
4. **Document any mathematical issues** discovered in the paper

### Step 4: Test Integration and Documentation

1. **Update TEST_STANDARD.md** with enhanced testing requirements
2. **Create example templates** for future equation-based testing
3. **Document any discrepancies** found between tests and paper
4. **Establish quality metrics** beyond just pass/fail rates

## Success Metrics

### Quantitative Measures
- **Number of actual equations tested** (vs. just dimensional checks)
- **Mathematical relationships verified** (vs. individual term dimensions)
- **Derivation steps validated** (vs. final result dimensions only)
- **Error detection capability** (test with deliberately introduced errors)

### Qualitative Measures
- **Confidence in mathematical framework** increases
- **Understanding of theoretical gaps** improves
- **Foundation for numerical work** strengthens
- **Experimental comparison readiness** enhances
