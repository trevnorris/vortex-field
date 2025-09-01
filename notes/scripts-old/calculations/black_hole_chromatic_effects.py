import sympy as sp
import numpy as np

# Part 1: Symbolic verification of the derivation
r, L, alpha = sp.symbols('r L alpha')
M = 1  # Normalized units

f = 1 - 2 * M / r
g = L**2 / r**2 + alpha / r
V = f * g
V_prime = sp.diff(V, r)

# Equations for circular orbit
eq1 = V - 1  # V(r_ph) = 1
eq2 = V_prime  # V'(r_ph) = 0

# Solve for L^2 from eq1
L2_from_eq1 = sp.solve(eq1, L**2)[0]  # L**2 = r**2 * (1/f - alpha / r)

# Substitute into eq2 and simplify to get the equation for r
eq2_sub = eq2.subs(L**2, L2_from_eq1)
eq_for_r = sp.simplify(eq2_sub * r**4)  # Multiply by r^4 to clear denominators, for readability

# The derived equation should be alpha * (r - 2)**2 - 2 * r**2 * (r - 3) = 0
derived_eq = alpha * (r - 2)**2 - 2 * r**2 * (r - 3)
print("Derived equation for r_ph (should be 0):", derived_eq)
print("SymPy-verified equation from substitution (expanded):", sp.expand(eq_for_r))  # Verify match (up to sign or factor)

# Symbolic solution for the cubic equation
cubic = 2 * r**3 - (6 + alpha) * r**2 + 4 * alpha * r - 4 * alpha
roots = sp.solve(cubic, r)

# Part 2: Numerical computation for sample alpha values
def compute_for_alpha(alpha_val):
    # Solve cubic numerically for r > 3
    cubic_num = lambda x: 2 * x**3 - (6 + alpha_val) * x**2 + 4 * alpha_val * x - 4 * alpha_val
    # Use numerical root finding (bisection, assuming one root >3)
    from scipy.optimize import bisect
    r_ph = bisect(cubic_num, 3.0001, 10)  # Search in (3+,10)

    # Compute b = sqrt(- r_ph^3 (r_ph -4) / (r_ph -2)^2 )
    b = np.sqrt( - (r_ph**3 * (r_ph - 4)) / (r_ph - 2)**2 )

    # Vacuum values
    r_vac = 3.0
    b_vac = 3 * np.sqrt(3)  # ≈5.19615

    # Shifts
    shift_r = (r_ph - r_vac) / r_vac * 100
    shift_b = (b - b_vac) / b_vac * 100

    return r_ph, b, shift_r, shift_b

# Sample calculations (replace with your alpha values)
alphas = [0, 0.05, 0.1, 0.2]
results = []
for a in alphas:
    if a == 0:
        results.append((3.0, 3*np.sqrt(3), 0.0, 0.0))
    else:
        results.append(compute_for_alpha(a))

# Print table
print("\n| α | r_ph | b | Shift from vacuum (%) |")
print("|---|------|---|-----------------------|")
for i, a in enumerate(alphas):
    r_ph, b, shift_r, shift_b = results[i]
    print(f"| {a} | {r_ph:.4f} | {b:.3f} | {shift_b:.2f} |")

# Optional: Verify with direct numerical solve of system for a sample alpha
# For alpha=0.1, solve V=1, V'=0 numerically for r and L
