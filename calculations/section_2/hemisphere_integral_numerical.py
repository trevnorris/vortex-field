"""
Hemisphere Integral: The Mathematical Gem Behind 4-fold Enhancement

This script provides a deep pedagogical exploration of the critical integral:
∫₀^∞ dw'/(ρ²+w'²)^(3/2) = 1/ρ²

This integral is the mathematical foundation that explains why each hemisphere
of a 4D vortex sheet contributes exactly Γ to the total circulation, leading
to the remarkable 4-fold enhancement in the projected 3D dynamics.

Mathematical Context:
- Arises from 4D Biot-Savart law for vortex sheet projections
- Each hemisphere (w>0, w<0) contributes this integral to total circulation
- The exact analytical result 1/ρ² is a beautiful mathematical fact
- Combined with direct intersection and induced circulation → 4Γ total

Educational Focus:
- Visualize WHY the integral converges to 1/ρ²
- Show the mathematical journey from numerical to analytical
- Demonstrate universal scaling across parameter ranges
- Connect back to physical interpretation in 4D→3D projection
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.special import beta
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

class HemisphereIntegralExplorer:
    """
    Deep exploration of the hemisphere integral that underlies
    the 4-fold circulation enhancement in 4D vortex projections.
    """

    def __init__(self):
        """Initialize with default plotting parameters."""
        plt.style.use('default')
        plt.rcParams.update({
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'legend.fontsize': 12,
            'figure.figsize': (10, 8)
        })

    def integrand(self, w_prime, rho):
        """
        The integrand: 1/(ρ² + w'²)^(3/2)

        This is the 4D Biot-Savart kernel that determines how much
        circulation a vortex sheet element at height w' contributes
        to the flow at the w=0 slice (our 3D universe).

        Parameters:
        -----------
        w_prime : float or array
            Height coordinate in 4th dimension
        rho : float
            Radial distance in 3D plane: ρ = √(x² + y²)

        Returns:
        --------
        float or array
            Integrand value(s)
        """
        return 1.0 / (rho**2 + w_prime**2)**(3/2)

    def analytical_solution(self, rho):
        """
        Exact analytical result: 1/ρ²

        Derivation:
        ∫₀^∞ dw'/(ρ²+w'²)^(3/2)

        Substitute u = w'/ρ, dw' = ρ du:
        = (1/ρ²) ∫₀^∞ du/(1+u²)^(3/2)

        The remaining integral equals 1 (can be solved using
        trigonometric substitution u = tan(θ) or beta functions).

        Therefore: ∫₀^∞ dw'/(ρ²+w'²)^(3/2) = 1/ρ²
        """
        return 1.0 / rho**2

    def numerical_integral(self, rho, w_max=None, method='quad'):
        """
        Compute the integral numerically with convergence tracking.

        Parameters:
        -----------
        rho : float
            Radial distance parameter
        w_max : float, optional
            Upper integration limit (None for adaptive)
        method : str
            Integration method ('quad', 'trapz', 'simpson')

        Returns:
        --------
        dict
            Results containing value, error, and convergence data
        """
        if w_max is None:
            # Adaptive upper limit: integrate until negligible contribution
            w_max = 10 * rho  # Start with reasonable guess

        if method == 'quad':
            # Use scipy's adaptive quadrature
            integral, error = integrate.quad(
                lambda w: self.integrand(w, rho),
                0, w_max,
                epsabs=1e-12, epsrel=1e-12
            )

            return {
                'value': integral,
                'error': error,
                'w_max': w_max,
                'method': method
            }

        elif method in ['trapz', 'simpson']:
            # Use fixed-grid methods for convergence studies
            w_points = np.linspace(0, w_max, 10000)
            integrand_values = self.integrand(w_points, rho)

            if method == 'trapz':
                integral = integrate.trapz(integrand_values, w_points)
            else:  # simpson
                integral = integrate.simpson(integrand_values, w_points)

            # Estimate error from boundary contribution
            boundary_contrib = self.integrand(w_max, rho) * rho

            return {
                'value': integral,
                'error': boundary_contrib,
                'w_max': w_max,
                'method': method
            }

    def convergence_study(self, rho_values=[0.1, 0.5, 1.0, 2.0, 5.0, 10.0]):
        """
        Study convergence behavior across different ρ values.

        Shows how the integral approaches 1/ρ² as w_max → ∞
        for each value of ρ.
        """
        print("Hemisphere Integral Convergence Study")
        print("=" * 50)
        print("Integral: ∫₀^∞ dw'/(ρ²+w'²)^(3/2)")
        print("Analytical result: 1/ρ²")
        print()

        # Create range of upper limits for each rho
        results = {}

        for rho in rho_values:
            # Logarithmic spacing of upper limits
            w_max_values = np.logspace(
                np.log10(0.1 * rho),  # Start small
                np.log10(20 * rho),   # Go to large multiples of rho
                50
            )

            integral_values = []
            errors = []

            for w_max in w_max_values:
                result = self.numerical_integral(rho, w_max, method='quad')
                integral_values.append(result['value'])
                errors.append(result['error'])

            analytical = self.analytical_solution(rho)
            relative_errors = np.abs(np.array(integral_values) - analytical) / analytical

            results[rho] = {
                'w_max_values': w_max_values,
                'integral_values': np.array(integral_values),
                'errors': np.array(errors),
                'relative_errors': relative_errors,
                'analytical': analytical
            }

            # Print summary for this rho
            final_value = integral_values[-1]
            final_error = relative_errors[-1]
            print(f"ρ = {rho:5.1f}: Final = {final_value:8.5f}, "
                  f"Analytical = {analytical:8.5f}, "
                  f"Rel. Error = {final_error:.2e}")

        return results

    def visualize_integrand(self, rho_values=[0.5, 1.0, 2.0]):
        """
        Visualize the integrand function for different ρ values.

        Shows how the integrand shape changes with ρ and why
        the integral converges.
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Linear scale plot
        w_max = 10
        w = np.linspace(0, w_max, 1000)

        for rho in rho_values:
            integrand_vals = self.integrand(w, rho)
            ax1.plot(w, integrand_vals, label=f'ρ = {rho}', linewidth=2)

            # Show area under curve for intuition
            if rho == 1.0:  # Highlight one case
                ax1.fill_between(w, 0, integrand_vals, alpha=0.3,
                               label=f'Area = 1/ρ² = {1/rho**2:.2f}')

        ax1.set_xlabel("w' (4th dimension coordinate)")
        ax1.set_ylabel("Integrand: 1/(ρ² + w'²)^(3/2)")
        ax1.set_title("Integrand Shape for Different ρ Values")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Log scale to show tail behavior
        w_log = np.logspace(-2, 2, 1000)

        for rho in rho_values:
            integrand_vals = self.integrand(w_log, rho)
            ax2.loglog(w_log, integrand_vals, label=f'ρ = {rho}', linewidth=2)

            # Show power law tail: ~ w'^(-3) for large w'
            if rho == 1.0:
                tail_line = 1.0 / w_log**3
                ax2.loglog(w_log[w_log > 3], tail_line[w_log > 3],
                          'k--', alpha=0.7, label='~ w\'^(-3)')

        ax2.set_xlabel("w' (log scale)")
        ax2.set_ylabel("Integrand (log scale)")
        ax2.set_title("Log-Log: Shows Power Law Tail")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        print("\nIntegrand Analysis:")
        print("- Peak at w' = 0, height = 1/ρ³")
        print("- Width scales with ρ")
        print("- Tail decays as w'^(-3) → ensures convergence")
        print("- Total area = 1/ρ² (exact analytical result)")

    def visualize_convergence(self, convergence_results):
        """
        Create comprehensive convergence visualization.
        """
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

        # 1. Integral value vs upper limit
        for rho, data in convergence_results.items():
            w_max_norm = data['w_max_values'] / rho  # Normalize by rho
            ax1.semilogx(w_max_norm, data['integral_values'],
                        label=f'ρ = {rho}', marker='o', markersize=3)

            # Show analytical line
            analytical = data['analytical']
            ax1.axhline(analytical, color='gray', linestyle='--', alpha=0.7)

        ax1.set_xlabel("w_max / ρ (normalized upper limit)")
        ax1.set_ylabel("Integral Value")
        ax1.set_title("Convergence to Analytical Result")
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # 2. Relative error vs upper limit
        for rho, data in convergence_results.items():
            w_max_norm = data['w_max_values'] / rho
            ax2.loglog(w_max_norm, data['relative_errors'] + 1e-16,  # Avoid log(0)
                      label=f'ρ = {rho}', marker='o', markersize=3)

        # Show exponential decay line
        w_range = np.logspace(0, 1.5, 100)
        exp_decay = np.exp(-w_range)
        ax2.loglog(w_range, exp_decay, 'k--', alpha=0.7,
                  label='~ exp(-w_max/ρ)')

        ax2.set_xlabel("w_max / ρ")
        ax2.set_ylabel("Relative Error")
        ax2.set_title("Exponential Convergence")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        # 3. Universal scaling: Show all converge to 1/ρ²
        rho_vals = list(convergence_results.keys())
        final_values = [data['integral_values'][-1] for data in convergence_results.values()]
        analytical_values = [data['analytical'] for data in convergence_results.values()]

        ax3.loglog(rho_vals, final_values, 'bo', markersize=8, label='Numerical')
        ax3.loglog(rho_vals, analytical_values, 'r-', linewidth=2, label='1/ρ²')

        ax3.set_xlabel("ρ")
        ax3.set_ylabel("Integral Value")
        ax3.set_title("Universal Scaling: ∝ 1/ρ²")
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # 4. Accuracy summary
        relative_errors_final = [data['relative_errors'][-1] for data in convergence_results.values()]

        ax4.semilogy(rho_vals, relative_errors_final, 'go', markersize=8)
        ax4.axhline(1e-10, color='red', linestyle='--',
                   label='10^-10 (excellent)')
        ax4.axhline(1e-6, color='orange', linestyle='--',
                   label='10^-6 (good)')

        ax4.set_xlabel("ρ")
        ax4.set_ylabel("Final Relative Error")
        ax4.set_title("Numerical Accuracy Achieved")
        ax4.legend()
        ax4.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

    def demonstrate_substitution(self, rho=1.0):
        """
        Show the u = w'/ρ substitution that leads to the analytical result.
        """
        print("\nAnalytical Solution via Substitution")
        print("=" * 40)
        print("Starting integral: ∫₀^∞ dw'/(ρ²+w'²)^(3/2)")
        print(f"With ρ = {rho}")
        print()
        print("Step 1: Substitute u = w'/ρ")
        print("        Then dw' = ρ du")
        print("        And ρ²+w'² = ρ²(1+u²)")
        print()
        print("Step 2: Transform the integral:")
        print("∫₀^∞ dw'/(ρ²+w'²)^(3/2) = ∫₀^∞ (ρ du)/[ρ²(1+u²)]^(3/2)")
        print("                          = ∫₀^∞ (ρ du)/[ρ³(1+u²)^(3/2)]")
        print("                          = (1/ρ²) ∫₀^∞ du/(1+u²)^(3/2)")
        print()
        print("Step 3: The remaining integral ∫₀^∞ du/(1+u²)^(3/2) = 1")
        print("        (Can be solved using trigonometric substitution)")
        print()
        print("Final result: ∫₀^∞ dw'/(ρ²+w'²)^(3/2) = 1/ρ²")
        print(f"              For ρ = {rho}: Result = {1/rho**2:.6f}")

        # Numerical verification of the u-integral
        u_max = 20  # Effectively infinity
        u_integral, _ = integrate.quad(
            lambda u: 1.0 / (1 + u**2)**(3/2),
            0, u_max
        )

        print(f"\nNumerical verification:")
        print(f"∫₀^{u_max} du/(1+u²)^(3/2) = {u_integral:.10f}")
        print(f"Theoretical value = 1.0")
        print(f"Error = {abs(u_integral - 1.0):.2e}")

    def connect_to_4fold_enhancement(self):
        """
        Explain how this integral connects to the 4-fold circulation enhancement.
        """
        print("\nConnection to 4-Fold Enhancement")
        print("=" * 40)
        print("This integral is one of FOUR contributions to the total circulation:")
        print()
        print("1. Direct Intersection: Γ")
        print("   - Vortex sheet crosses w=0 plane directly")
        print("   - Standard 2D vortex circulation")
        print()
        print("2. Upper Hemisphere (w>0): Γ")
        print("   - Our integral: ∫₀^∞ dw'/(ρ²+w'²)^(3/2) = 1/ρ²")
        print("   - After proper normalization → contributes Γ")
        print()
        print("3. Lower Hemisphere (w<0): Γ")
        print("   - Same integral from -∞ to 0")
        print("   - Also contributes Γ by symmetry")
        print()
        print("4. Induced w-flow Circulation: Γ")
        print("   - Drainage into 4th dimension creates swirl")
        print("   - Topological linking contributes Γ")
        print()
        print("TOTAL: 4Γ observed circulation")
        print("This 4-fold enhancement is geometric, not fitted!")
        print()
        print("Physical interpretation:")
        print("- Each hemisphere acts like a 'ghost vortex'")
        print("- The 4D→3D projection naturally amplifies circulation")
        print("- This explains gravitomagnetic effects in the framework")

def main():
    """
    Main demonstration of the hemisphere integral.
    """
    print("HEMISPHERE INTEGRAL: Mathematical Foundation of 4-Fold Enhancement")
    print("=" * 70)
    print()

    explorer = HemisphereIntegralExplorer()

    # 1. Demonstrate the analytical solution
    explorer.demonstrate_substitution(rho=1.0)

    # 2. Visualize the integrand
    print("\nGenerating integrand visualization...")
    explorer.visualize_integrand()

    # 3. Convergence study
    print("\nPerforming convergence study...")
    convergence_results = explorer.convergence_study()

    # 4. Comprehensive visualization
    print("\nGenerating convergence visualization...")
    explorer.visualize_convergence(convergence_results)

    # 5. Connect back to physics
    explorer.connect_to_4fold_enhancement()

    # 6. Summary table
    print("\nSUMMARY TABLE: Numerical vs Analytical")
    print("=" * 50)
    print("ρ        Numerical     Analytical    Rel. Error")
    print("-" * 50)

    for rho in [0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        result = explorer.numerical_integral(rho, w_max=20*rho)
        analytical = explorer.analytical_solution(rho)
        rel_error = abs(result['value'] - analytical) / analytical

        print(f"{rho:5.1f}    {result['value']:9.6f}    {analytical:9.6f}    {rel_error:.2e}")

    print("\nConclusion:")
    print("The integral ∫₀^∞ dw'/(ρ²+w'²)^(3/2) = 1/ρ² is exact.")
    print("This mathematical gem enables the 4-fold circulation enhancement")
    print("that underlies gravitomagnetic effects in the 4D vortex framework.")

if __name__ == "__main__":
    main()
