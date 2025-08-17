import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import quad
import sympy as sp

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
        Find photon path using ray-tracing with Snell's law in varying medium.
        This is more stable than the optimization approach.
        """
        if x_range is None:
            x_range = (-10 * R_sun, 10 * R_sun)
        
        # Use ray-tracing approach: integrate the ray equation
        # dn/dr * cos(θ) = d/ds(n * sin(θ)) where θ is angle from radial
        
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
            Ray equations in medium with varying refractive index
            dx/ds = vx/|v|, dy/ds = vy/|v|
            d(n*vx)/ds = n_x, d(n*vy)/ds = n_y (where n_x = ∂n/∂x)
            """
            x, y, vx, vy = state
            
            r = np.sqrt(x**2 + y**2)
            r = max(r, 0.1 * R_sun)  # Avoid singularity
            
            n = self.effective_refractive_index(r)
            
            # Calculate gradient of refractive index
            dr = 0.01 * R_sun  # Small step for numerical derivative
            r_plus = r + dr
            r_minus = max(r - dr, 0.1 * R_sun)
            
            dn_dr = (self.effective_refractive_index(r_plus) - 
                    self.effective_refractive_index(r_minus)) / (2 * dr)
            
            # ∇n = (dn/dr) * r̂ = (dn/dr) * (x/r, y/r)
            dn_dx = dn_dr * x / r if r > 0 else 0
            dn_dy = dn_dr * y / r if r > 0 else 0
            
            # Normalize velocity
            v_mag = np.sqrt(vx**2 + vy**2)
            if v_mag < 1e-10:
                v_mag = 1.0
            
            dx_ds = vx / v_mag
            dy_ds = vy / v_mag
            
            # Ray equations: d(n*v)/ds = ∇n
            dvx_ds = dn_dx / n - vx * (vx * dn_dx + vy * dn_dy) / (n * v_mag**2)
            dvy_ds = dn_dy / n - vy * (vx * dn_dx + vy * dn_dy) / (n * v_mag**2)
            
            return [dx_ds, dy_ds, dvx_ds, dvy_ds]
        
        # Integrate the ray
        s_span = (0, 2 * abs(x_range[1] - x_range[0]))  # Path length parameter
        s_eval = np.linspace(0, s_span[1], 10000)
        
        try:
            from scipy.integrate import solve_ivp
            
            sol = solve_ivp(ray_equations, s_span, initial_state, 
                          t_eval=s_eval, method='RK45', 
                          rtol=1e-10, atol=1e-12, max_step=R_sun/100)
            
            if sol.success:
                return sol.t, sol.y, True
            else:
                return None, None, False
                
        except Exception as e:
            print(f"Integration failed: {e}")
            return None, None, False
    
    def calculate_deflection_angle(self, s_vals, path_data):
        """
        Calculate deflection angle from ray-traced path.
        
        Deflection convention: positive = bending toward the mass
        For grazing ray above Sun, bending downward = positive deflection
        """
        if s_vals is None or path_data is None:
            return 0.0
        
        x_vals, y_vals, vx_vals, vy_vals = path_data
        
        # Find asymptotic directions (far from the mass)
        # Use points far from origin where curvature is minimal
        n_points = len(x_vals)
        n_edge = max(10, n_points // 20)  # Use outer 5% of points
        
        # Initial direction (average over first n_edge points)
        initial_vx = np.mean(vx_vals[:n_edge])
        initial_vy = np.mean(vy_vals[:n_edge])
        initial_angle = np.arctan2(initial_vy, initial_vx)
        
        # Final direction (average over last n_edge points)  
        final_vx = np.mean(vx_vals[-n_edge:])
        final_vy = np.mean(vy_vals[-n_edge:])
        final_angle = np.arctan2(final_vy, final_vx)
        
        # Raw deflection angle
        deflection_raw = final_angle - initial_angle
        
        # Handle angle wrapping around ±π
        if deflection_raw > np.pi:
            deflection_raw -= 2 * np.pi
        elif deflection_raw < -np.pi:
            deflection_raw += 2 * np.pi
        
        # For gravitational deflection: bending toward mass = positive
        # Our setup has mass at origin, ray passing above (positive y)
        # Bending downward (negative change in angle) should be positive deflection
        deflection_rad = -deflection_raw  # Sign flip for correct convention
        
        # Convert to arcseconds
        deflection_arcsec = deflection_rad * 206265
        
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
    Main function to run photon deflection calculation
    """
    print("=" * 60)
    print("PHOTON DEFLECTION IN AETHER FRAMEWORK")
    print("=" * 60)
    
    # Initialize calculator
    calc = AetherPhotonDeflection(mass=M_sun, include_1pn=True)
    
    # Parameters
    impact_parameter = 1.1 * R_sun  # Just grazing the Sun
    x_range = (-10 * R_sun, 10 * R_sun)
    
    print(f"Central mass: {M_sun:.3e} kg (Sun)")
    print(f"Impact parameter: {impact_parameter/R_sun:.2f} R_sun")
    print(f"Integration range: {x_range[0]/R_sun:.1f} to {x_range[1]/R_sun:.1f} R_sun")
    print("\nCalculating photon path using ray-tracing in varying aether...")
    
    # Find optimal path using ray-tracing
    s_vals, path_data, success = calc.find_photon_path(impact_parameter, x_range)
    
    if not success or s_vals is None:
        print("ERROR: Ray-tracing failed!")
        return None, None, None
    
    # Calculate deflection
    deflection_angle = calc.calculate_deflection_angle(s_vals, path_data)
    
    # Theoretical predictions
    gr_prediction = 4 * G * M_sun / (c**2 * impact_parameter) * 206265  # GR result in arcsec
    newton_prediction = gr_prediction / 2  # Newtonian result
    
    print("\n" + "=" * 40)
    print("RESULTS:")
    print("=" * 40)
    print(f"Aether framework deflection: {deflection_angle:.3f} arcseconds")
    print(f"General Relativity prediction: {gr_prediction:.3f} arcseconds")
    print(f"Newtonian prediction: {newton_prediction:.3f} arcseconds")
    print(f"Ratio (Aether/GR): {deflection_angle/gr_prediction:.4f}")
    print(f"Difference from GR: {abs(deflection_angle - gr_prediction):.3f} arcseconds")
    
    # Physical insights
    print("\n" + "=" * 40)
    print("PHYSICAL INSIGHTS:")
    print("=" * 40)
    
    r_min = impact_parameter
    phi_at_closest = calc.gravitational_potential(r_min)
    n_eff_at_closest = calc.effective_refractive_index(r_min)
    
    print(f"At closest approach (r = {r_min/R_sun:.2f} R_sun):")
    print(f"  Gravitational potential: Φ/c² = {phi_at_closest/c**2:.2e}")
    print(f"  4-fold enhanced coupling: 4|Φ|/c² = {-4*phi_at_closest/c**2:.2e}")
    print(f"  Aether refractive index: n_eff = {n_eff_at_closest:.6f}")
    print(f"  Aether speed: v_eff/c = {1/n_eff_at_closest:.6f}")
    print(f"  Speed reduction: Δv/c = {(1 - 1/n_eff_at_closest)*1e6:.2f} ppm")
    
    # Create visualization
    print("\nGenerating visualization...")
    fig = calc.plot_results(s_vals, path_data, impact_parameter, deflection_angle, x_range)
    
    return deflection_angle, gr_prediction, fig

if __name__ == "__main__":
    # Run the calculation
    deflection, gr_pred, figure = run_deflection_calculation()
    
    # Additional analysis with different impact parameters
    print("\n" + "=" * 60)
    print("ANALYSIS WITH DIFFERENT IMPACT PARAMETERS:")
    print("=" * 60)
    
    calc = AetherPhotonDeflection(mass=M_sun, include_1pn=True)
    
    b_values = np.array([1.1, 1.5, 2.0, 3.0, 5.0]) * R_sun
    deflections = []
    gr_predictions = []
    
    for b in b_values:
        print(f"\nCalculating for b = {b/R_sun:.1f} R_sun...")
        try:
            s_vals, path_data, success = calc.find_photon_path(b, (-10*R_sun, 10*R_sun))
            if success and s_vals is not None:
                defl = calc.calculate_deflection_angle(s_vals, path_data)
                gr_pred = 4 * G * M_sun / (c**2 * b) * 206265
                
                deflections.append(defl)
                gr_predictions.append(gr_pred)
                
                print(f"  Aether: {defl:.3f}\", GR: {gr_pred:.3f}\", Ratio: {defl/gr_pred:.4f}")
            else:
                print(f"  Ray-tracing failed for b = {b/R_sun:.1f} R_sun")
        except Exception as e:
            print(f"  Error for b = {b/R_sun:.1f} R_sun: {e}")
    
    # Plot comparison
    if deflections:
        plt.figure(figsize=(10, 6))
        
        valid_b = b_values[:len(deflections)]
        
        plt.subplot(1, 2, 1)
        plt.loglog(valid_b/R_sun, deflections, 'bo-', label='Aether Framework', linewidth=2)
        plt.loglog(valid_b/R_sun, gr_predictions[:len(deflections)], 'r--', 
                  label='General Relativity', linewidth=2)
        plt.xlabel('Impact Parameter (R_sun)')
        plt.ylabel('Deflection Angle (arcsec)')
        plt.title('Light Deflection vs Impact Parameter')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.subplot(1, 2, 2)
        ratios = np.array(deflections) / np.array(gr_predictions[:len(deflections)])
        plt.semilogx(valid_b/R_sun, ratios, 'go-', linewidth=2)
        plt.axhline(y=1, color='r', linestyle='--', alpha=0.7, label='Perfect Agreement')
        plt.xlabel('Impact Parameter (R_sun)')
        plt.ylabel('Aether/GR Ratio')
        plt.title('Agreement with General Relativity')
        plt.legend()
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
        
    print("\n" + "=" * 60)
    print("CALCULATION COMPLETE!")
    print("=" * 60)
