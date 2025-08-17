import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import quad
import sympy as sp
import time

# Physical constants
G = 6.67430e-11  # m^3 kg^-1 s^-2
c = 2.99792458e8  # m/s
M_sun = 1.989e30  # kg
R_sun = 6.957e8   # m

class AetherPhotonDeflection:
    """
    Calculates photon deflection using Fermat's principle in the aether framework.
    
    Key concepts:
    - Photons travel along "straight" paths in aether
    - Aether flows toward gravitational sources (vortex sinks) 
    - Aether density varies, creating position-dependent refractive index
    - Path bending comes from varying v_eff(r) with 4-fold enhancement
    
    Critical corrections implemented:
    - 4-fold geometric enhancement from 4D→3D projection (Section 2.7)
    - Correct sign convention: positive deflection = bending toward mass
    - Ray-tracing with Snell's law for stable numerics
    """
    
    def __init__(self, mass=M_sun, include_1pn=True):
        self.mass = mass
        self.include_1pn = include_1pn
        self.rs = 2 * G * mass / c**2  # Schwarzschild radius
        
    def gravitational_potential(self, r):
        """
        Scalar potential from aether framework (Section 4.4):
        Φ(r) = -GM/r + (GM)²/(2c²r²) + O(ε³)
        
        The first term is Newtonian, second is 1 PN correction
        """
        phi_newton = -G * self.mass / r
        
        if self.include_1pn:
            phi_1pn = (G * self.mass)**2 / (2 * c**2 * r**2)
            return phi_newton + phi_1pn
        else:
            return phi_newton
    
    def effective_refractive_index(self, r):
        """
        Aether refractive index from density variation with 4-fold enhancement:
        
        From Section 2.7: 4D→3D projection yields 4-fold geometric enhancement
        This affects the refractive index coupling to the gravitational potential.
        
        n_eff(r) = 1 + 4|Φ|/(2c²) = 1 - 2Φ/c² (since Φ < 0 near masses)
        """
        phi = self.gravitational_potential(r)
        
        # Include the 4-fold enhancement from 4D→3D projection (Section 2.7)
        # Factor of 4 comes from geometric projection of vortex sheets
        phi_over_c2 = -4 * phi / (2 * c**2)  # 4-fold enhancement + sign correction
        
        # This gives n > 1 near masses (light slows down)
        return 1 + phi_over_c2
    
    def path_parameterization(self, t, y_params):
        """
        Parameterize photon path as y(x) where x = t
        y_params: coefficients for polynomial path approximation
        """
        # Use polynomial approximation for smooth path
        x = t
        y = sum(y_params[i] * x**i for i in range(len(y_params)))
        return x, y
    
    def optical_path_integrand(self, x, y, dydx):
        """
        Integrand for optical path length: n_eff(r) * ds
        where ds = sqrt(1 + (dy/dx)²) dx
        """
        r = np.sqrt(x**2 + y**2)
        
        # Avoid singularity at origin
        r = max(r, 0.1 * R_sun)
        
        n_eff = self.effective_refractive_index(r)
        ds = np.sqrt(1 + dydx**2)
        
        return n_eff * ds
    
    def total_optical_path(self, y_params, x_range):
        """
        Calculate total optical path length for given path parameters
        """
        x_start, x_end = x_range
        
        def integrand(x):
            # Calculate y and dy/dx from parameters
            y = sum(y_params[i] * x**i for i in range(len(y_params)))
            dydx = sum(i * y_params[i] * x**(i-1) for i in range(1, len(y_params)))
            
            return self.optical_path_integrand(x, y, dydx)
        
        try:
            path_length, _ = quad(integrand, x_start, x_end, limit=100)
            return path_length
        except:
            return 1e10  # Return large value for invalid paths
    
    def find_photon_path(self, impact_parameter, x_range=None):
        """
        Find photon path using high-precision ray-tracing with Snell's law.
        Enhanced for maximum numerical accuracy.
        """
        if x_range is None:
            x_range = (-20 * R_sun, 20 * R_sun)  # Larger range for better asymptotics
        
        # Initial conditions: photon starts far away, moving horizontally
        x0 = x_range[0]
        y0 = impact_parameter
        
        # Initial direction (toward +x)
        vx0 = 1.0
        vy0 = 0.0
        
        # Initial position and velocity
        initial_state = [x0, y0, vx0, vy0]
        
        def ray_equations(s, state):
            """
            High-precision ray equations in medium with varying refractive index.
            Enhanced with better numerical derivatives and singularity handling.
            """
            x, y, vx, vy = state
            
            r = np.sqrt(x**2 + y**2)
            r = max(r, 0.01 * R_sun)  # Tighter singularity avoidance
            
            n = self.effective_refractive_index(r)
            
            # High-precision gradient calculation using central differences
            # Use smaller step size and higher-order method
            dr = min(0.001 * R_sun, 0.001 * r)  # Adaptive step size
            
            # 5-point stencil for better accuracy
            r_vals = np.array([r - 2*dr, r - dr, r, r + dr, r + 2*dr])
            r_vals = np.maximum(r_vals, 0.01 * R_sun)  # Avoid singularity
            
            n_vals = np.array([self.effective_refractive_index(ri) for ri in r_vals])
            
            # 5-point central difference formula: more accurate than 3-point
            dn_dr = (-n_vals[4] + 8*n_vals[3] - 8*n_vals[1] + n_vals[0]) / (12 * dr)
            
            # ∇n = (dn/dr) * r̂ = (dn/dr) * (x/r, y/r)
            dn_dx = dn_dr * x / r if r > 0 else 0
            dn_dy = dn_dr * y / r if r > 0 else 0
            
            # Normalize velocity with higher precision
            v_mag = np.sqrt(vx**2 + vy**2)
            if v_mag < 1e-15:  # Tighter tolerance
                v_mag = 1.0
            
            dx_ds = vx / v_mag
            dy_ds = vy / v_mag
            
            # Ray equations: d(n*v)/ds = ∇n with higher precision
            # Include second-order terms for better accuracy
            v_dot_grad_n = vx * dn_dx + vy * dn_dy
            
            dvx_ds = dn_dx / n - vx * v_dot_grad_n / (n * v_mag**2)
            dvy_ds = dn_dy / n - vy * v_dot_grad_n / (n * v_mag**2)
            
            return [dx_ds, dy_ds, dvx_ds, dvy_ds]
        
        # High-precision integration settings
        s_span = (0, 3 * abs(x_range[1] - x_range[0]))  # Longer path for better convergence
        
        # More evaluation points for smoother result
        n_points = 50000  # Increased from 10000
        s_eval = np.linspace(0, s_span[1], n_points)
        
        try:
            from scipy.integrate import solve_ivp
            
            # High-precision solver settings
            sol = solve_ivp(ray_equations, s_span, initial_state, 
                          t_eval=s_eval, 
                          method='DOP853',  # 8th-order Runge-Kutta (higher than RK45)
                          rtol=1e-14,       # Much tighter relative tolerance
                          atol=1e-16,       # Much tighter absolute tolerance
                          max_step=R_sun/1000,  # Smaller maximum step size
                          dense_output=True)    # For smooth interpolation
            
            if sol.success:
                return sol.t, sol.y, True
            else:
                print(f"Integration failed: {sol.message}")
                return None, None, False
                
        except Exception as e:
            print(f"High-precision integration failed: {e}")
            return None, None, False
    
    def calculate_deflection_angle(self, s_vals, path_data):
        """
        High-precision deflection angle calculation from ray-traced path.
        Uses sophisticated asymptotic analysis for maximum accuracy.
        """
        if s_vals is None or path_data is None:
            return 0.0
        
        x_vals, y_vals, vx_vals, vy_vals = path_data
        n_points = len(x_vals)
        
        # Use larger asymptotic regions for better accuracy
        n_edge = max(100, n_points // 10)  # Use outer 10% of points
        
        # Find points that are sufficiently far from the mass for asymptotic analysis
        r_vals = np.sqrt(x_vals**2 + y_vals**2)
        
        # Initial asymptotic region: far left, away from Sun
        initial_mask = (r_vals[:n_edge] > 5 * R_sun)
        if np.sum(initial_mask) < 10:  # Fallback if not enough far points
            initial_mask = np.arange(n_edge) < min(50, n_edge)
        
        # Final asymptotic region: far right, away from Sun  
        final_mask = (r_vals[-n_edge:] > 5 * R_sun)
        if np.sum(final_mask) < 10:  # Fallback if not enough far points
            final_mask = np.arange(max(0, n_points-n_edge), n_points) >= max(0, n_points-50)
        
        # Calculate initial direction using weighted average (more weight to farther points)
        initial_indices = np.where(initial_mask)[0] if np.any(initial_mask) else np.arange(min(50, n_edge))
        weights_initial = r_vals[initial_indices]**2  # Weight by r² for better asymptotic behavior
        weights_initial /= np.sum(weights_initial)
        
        initial_vx = np.average(vx_vals[initial_indices], weights=weights_initial)
        initial_vy = np.average(vy_vals[initial_indices], weights=weights_initial)
        
        # Calculate final direction using weighted average
        final_indices = np.where(final_mask)[0] + (n_points - n_edge) if np.any(final_mask) else np.arange(max(0, n_points-50), n_points)
        weights_final = r_vals[final_indices]**2  # Weight by r² for better asymptotic behavior
        weights_final /= np.sum(weights_final)
        
        final_vx = np.average(vx_vals[final_indices], weights=weights_final)
        final_vy = np.average(vy_vals[final_indices], weights=weights_final)
        
        # Calculate angles with higher precision
        initial_angle = np.arctan2(initial_vy, initial_vx)
        final_angle = np.arctan2(final_vy, final_vx)
        
        # Raw deflection angle with careful handling of branch cuts
        deflection_raw = final_angle - initial_angle
        
        # Handle angle wrapping around ±π more carefully
        while deflection_raw > np.pi:
            deflection_raw -= 2 * np.pi
        while deflection_raw < -np.pi:
            deflection_raw += 2 * np.pi
        
        # For gravitational deflection: bending toward mass = positive
        # Our setup has mass at origin, ray passing above (positive y)
        # Bending downward (negative change in angle) should be positive deflection
        deflection_rad = -deflection_raw  # Sign flip for correct convention
        
        # Convert to arcseconds with high precision
        deflection_arcsec = deflection_rad * 206264.806247096  # More precise conversion factor
        
        return deflection_arcsec
    
    def generate_path_points(self, optimal_params, x_range, n_points=1000):
        """
        Generate points along the optimal path for plotting
        """
        x_vals = np.linspace(x_range[0], x_range[1], n_points)
        y_vals = np.array([sum(optimal_params[i] * x**i for i in range(len(optimal_params))) 
                          for x in x_vals])
        
        return x_vals, y_vals
    
    def plot_results(self, s_vals, path_data, impact_parameter, deflection_angle, x_range):
        """
        Create visualization of photon path and aether properties
        """
        if s_vals is None or path_data is None:
            print("No valid path data to plot!")
            return None
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        
        x_vals, y_vals, vx_vals, vy_vals = path_data
        
        # 1. Photon path
        ax1.plot(x_vals/R_sun, y_vals/R_sun, 'b-', linewidth=2, label='Photon path')
        ax1.plot([x_range[0]/R_sun, x_range[1]/R_sun], 
                [impact_parameter/R_sun, impact_parameter/R_sun], 
                'r--', alpha=0.5, label='Straight line')
        
        # Add Sun
        circle = plt.Circle((0, 0), 1, color='orange', alpha=0.7)
        ax1.add_patch(circle)
        
        ax1.set_xlabel('x / R_sun')
        ax1.set_ylabel('y / R_sun')
        ax1.set_title(f'Photon Deflection (δφ = {deflection_angle:.3f}")')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_aspect('equal')
        ax1.set_xlim(x_range[0]/R_sun, x_range[1]/R_sun)
        
        # 2. Gravitational potential
        r_vals = np.linspace(0.5*R_sun, 10*R_sun, 1000)
        phi_vals = [self.gravitational_potential(r) for r in r_vals]
        
        ax2.plot(r_vals/R_sun, np.array(phi_vals)/c**2, 'g-', linewidth=2)
        ax2.set_xlabel('r / R_sun')
        ax2.set_ylabel('Φ/c²')
        ax2.set_title('Gravitational Potential')
        ax2.grid(True, alpha=0.3)
        
        # 3. Refractive index
        n_eff_vals = [self.effective_refractive_index(r) for r in r_vals]
        
        ax3.plot(r_vals/R_sun, n_eff_vals, 'm-', linewidth=2)
        ax3.set_xlabel('r / R_sun')
        ax3.set_ylabel('n_eff(r)')
        ax3.set_title('Aether Refractive Index')
        ax3.grid(True, alpha=0.3)
        
        # 4. Path deviation from straight line
        straight_line = np.full_like(y_vals, impact_parameter)
        deviation = (y_vals - straight_line) / R_sun
        
        ax4.plot(x_vals/R_sun, deviation * 1e6, 'c-', linewidth=2)  # Convert to micro R_sun
        ax4.set_xlabel('x / R_sun')
        ax4.set_ylabel('Path deviation (μR_sun)')
        ax4.set_title('Deviation from Straight Line')
        ax4.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        return fig

def run_deflection_calculation():
    """
    Main function to run high-precision photon deflection calculation
    """
    print("=" * 70)
    print("HIGH-PRECISION PHOTON DEFLECTION IN AETHER FRAMEWORK")
    print("=" * 70)
    
    # Initialize calculator
    calc = AetherPhotonDeflection(mass=M_sun, include_1pn=True)
    
    # Parameters for high precision
    impact_parameter = 1.1 * R_sun  # Just grazing the Sun
    x_range = (-20 * R_sun, 20 * R_sun)  # Larger range for better asymptotics
    
    print(f"Central mass: {M_sun:.6e} kg (Sun)")
    print(f"Impact parameter: {impact_parameter/R_sun:.3f} R_sun")
    print(f"Integration range: {x_range[0]/R_sun:.1f} to {x_range[1]/R_sun:.1f} R_sun")
    print(f"Using 8th-order Runge-Kutta with 50,000 points")
    print(f"Tolerance: rtol=1e-14, atol=1e-16")
    print("\nCalculating photon path using high-precision ray-tracing...")
    
    # Find optimal path using high-precision ray-tracing
    import time
    start_time = time.time()
    
    s_vals, path_data, success = calc.find_photon_path(impact_parameter, x_range)
    
    computation_time = time.time() - start_time
    print(f"Computation completed in {computation_time:.2f} seconds")
    
    if not success or s_vals is None:
        print("ERROR: High-precision ray-tracing failed!")
        return None, None, None
    
    # Calculate deflection
    deflection_angle = calc.calculate_deflection_angle(s_vals, path_data)
    
    # High-precision theoretical predictions
    gr_prediction = 4 * G * M_sun / (c**2 * impact_parameter) * 206264.806247096
    newton_prediction = gr_prediction / 2  # Newtonian result
    
    print("\n" + "=" * 50)
    print("HIGH-PRECISION RESULTS:")
    print("=" * 50)
    print(f"Aether framework deflection: {deflection_angle:.6f} arcseconds")
    print(f"General Relativity prediction: {gr_prediction:.6f} arcseconds")
    print(f"Newtonian prediction: {newton_prediction:.6f} arcseconds")
    print(f"Ratio (Aether/GR): {deflection_angle/gr_prediction:.8f}")
    print(f"Difference from GR: {abs(deflection_angle - gr_prediction):.6f} arcseconds")
    print(f"Relative error: {abs(deflection_angle - gr_prediction)/gr_prediction*100:.4f}%")
    
    # Enhanced physical insights
    print("\n" + "=" * 50)
    print("ENHANCED PHYSICAL INSIGHTS:")
    print("=" * 50)
    
    r_min = impact_parameter
    phi_at_closest = calc.gravitational_potential(r_min)
    n_eff_at_closest = calc.effective_refractive_index(r_min)
    
    print(f"At closest approach (r = {r_min/R_sun:.3f} R_sun):")
    print(f"  Gravitational potential: Φ/c² = {phi_at_closest/c**2:.3e}")
    print(f"  4-fold enhanced coupling: 4|Φ|/c² = {-4*phi_at_closest/c**2:.3e}")
    print(f"  Aether refractive index: n_eff = {n_eff_at_closest:.9f}")
    print(f"  Aether speed: v_eff/c = {1/n_eff_at_closest:.9f}")
    print(f"  Speed reduction: Δv/c = {(1 - 1/n_eff_at_closest)*1e6:.3f} ppm")
    
    # Path analysis
    x_vals, y_vals, vx_vals, vy_vals = path_data
    r_vals = np.sqrt(x_vals**2 + y_vals**2)
    r_min_actual = np.min(r_vals)
    closest_idx = np.argmin(r_vals)
    
    print(f"\nPath analysis:")
    print(f"  Actual closest approach: {r_min_actual/R_sun:.6f} R_sun")
    print(f"  Maximum path deviation: {np.max(np.abs(y_vals - impact_parameter))/R_sun*1e6:.3f} μR_sun")
    print(f"  Integration points used: {len(s_vals):,}")
    
    # Create visualization
    print("\nGenerating high-resolution visualization...")
    fig = calc.plot_results(s_vals, path_data, impact_parameter, deflection_angle, x_range)
    
    return deflection_angle, gr_prediction, fig

if __name__ == "__main__":
    # Run the calculation
    deflection, gr_pred, figure = run_deflection_calculation()
    
    # High-precision analysis with different impact parameters
    print("\n" + "=" * 70)
    print("HIGH-PRECISION ANALYSIS WITH DIFFERENT IMPACT PARAMETERS:")
    print("=" * 70)
    
    calc = AetherPhotonDeflection(mass=M_sun, include_1pn=True)
    
    # Test more impact parameters for better statistics
    b_values = np.array([1.1, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]) * R_sun
    deflections = []
    gr_predictions = []
    computation_times = []
    
    for i, b in enumerate(b_values):
        print(f"\nCalculating {i+1}/{len(b_values)}: b = {b/R_sun:.1f} R_sun...")
        try:
            start_time = time.time()
            s_vals, path_data, success = calc.find_photon_path(b, (-20*R_sun, 20*R_sun))
            comp_time = time.time() - start_time
            
            if success and s_vals is not None:
                defl = calc.calculate_deflection_angle(s_vals, path_data)
                gr_pred = 4 * G * M_sun / (c**2 * b) * 206264.806247096
                
                deflections.append(defl)
                gr_predictions.append(gr_pred)
                computation_times.append(comp_time)
                
                rel_error = abs(defl - gr_pred) / gr_pred * 100
                print(f"  Aether: {defl:.6f}\", GR: {gr_pred:.6f}\", Ratio: {defl/gr_pred:.8f}")
                print(f"  Relative error: {rel_error:.4f}%, Time: {comp_time:.1f}s")
            else:
                print(f"  High-precision ray-tracing failed for b = {b/R_sun:.1f} R_sun")
        except Exception as e:
            print(f"  Error for b = {b/R_sun:.1f} R_sun: {e}")
    
    print(f"\nTotal computation time: {sum(computation_times):.1f} seconds")
    print(f"Average time per calculation: {np.mean(computation_times):.1f} seconds")
    
    # High-precision comparison plot
    if deflections:
        plt.figure(figsize=(15, 8))
        
        valid_b = b_values[:len(deflections)]
        
        plt.subplot(1, 3, 1)
        plt.loglog(valid_b/R_sun, deflections, 'bo-', label='Aether Framework', linewidth=2, markersize=8)
        plt.loglog(valid_b/R_sun, gr_predictions[:len(deflections)], 'r--', 
                  label='General Relativity', linewidth=2)
        plt.xlabel('Impact Parameter (R_sun)', fontsize=12)
        plt.ylabel('Deflection Angle (arcsec)', fontsize=12)
        plt.title('Light Deflection vs Impact Parameter\n(High Precision)', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 3, 2)
        ratios = np.array(deflections) / np.array(gr_predictions[:len(deflections)])
        plt.semilogx(valid_b/R_sun, ratios, 'go-', linewidth=2, markersize=8)
        plt.axhline(y=1, color='r', linestyle='--', alpha=0.7, label='Perfect Agreement')
        plt.xlabel('Impact Parameter (R_sun)', fontsize=12)
        plt.ylabel('Aether/GR Ratio', fontsize=12)
        plt.title('Agreement with General Relativity\n(High Precision)', fontsize=14)
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        # New subplot: Relative error
        plt.subplot(1, 3, 3)
        rel_errors = np.abs(np.array(deflections) - np.array(gr_predictions[:len(deflections)])) / np.array(gr_predictions[:len(deflections)]) * 100
        plt.semilogy(valid_b/R_sun, rel_errors, 'mo-', linewidth=2, markersize=8)
        plt.xlabel('Impact Parameter (R_sun)', fontsize=12)
        plt.ylabel('Relative Error (%)', fontsize=12)
        plt.title('Relative Error vs Impact Parameter\n(High Precision)', fontsize=14)
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
        # Statistical summary
        print(f"\nSTATISTICAL SUMMARY:")
        print(f"Average ratio (Aether/GR): {np.mean(ratios):.8f}")
        print(f"Standard deviation of ratios: {np.std(ratios):.8f}")
        print(f"Average relative error: {np.mean(rel_errors):.4f}%")
        print(f"Maximum relative error: {np.max(rel_errors):.4f}%")
        print(f"Minimum relative error: {np.min(rel_errors):.4f}%")
        
    print("\n" + "=" * 70)
    print("HIGH-PRECISION CALCULATION COMPLETE!")
    print("=" * 70)
