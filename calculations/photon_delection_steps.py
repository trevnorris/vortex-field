#!/usr/bin/env python3
"""
Incremental Aether Framework Implementation
==========================================

Starting from the WORKING version (99.6% agreement) and adding components incrementally:

BASELINE (Working):
✓ Scalar potential Φ(r) with 1 PN corrections
✓ 4-fold geometric enhancement  
✓ Refractive index n(r) = 1 + 4|Φ|/(2c²)
✓ Simple ray tracing with Snell's law
Result: 1.584 arcsec (99.61% of GR)

INCREMENTAL ADDITIONS:
□ Step 1: Add vector potential A(r,θ) from solar rotation
□ Step 2: Add "suck" flow effects
□ Step 3: Add "swirl" flow effects
□ Step 4: Add complete flow coupling

Each step is tested independently to see the effect.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

# Constants
c = 2.99792458e8      # Speed of light (m/s)
G = 6.67430e-11       # Gravitational constant
M_sun = 1.98847e30    # Solar mass (kg)
R_sun = 6.96342e8     # Solar radius (m)
J_sun = 1.92e42       # Solar angular momentum (kg⋅m²/s)

class IncrementalAetherFramework:
    """
    Start with working version, add components incrementally
    """
    
    def __init__(self, mass=M_sun, spin=J_sun, include_1pn=True):
        self.mass = mass
        self.spin = spin
        self.include_1pn = include_1pn
        
        # Control flags for incremental testing
        self.include_vector_potential = False
        self.include_suck_flow = False  
        self.include_swirl_flow = False
        self.include_flow_coupling = False
        
        # Derived parameters
        self.rs = 2 * G * mass / c**2
        
        print(f"Incremental Aether Framework initialized:")
        print(f"  Mass: {self.mass:.2e} kg")
        print(f"  Spin: {self.spin:.2e} kg⋅m²/s")
        print(f"  Schwarzschild radius: {self.rs:.1f} m")
        print(f"  Include 1 PN: {self.include_1pn}")
        
    def set_step(self, step_name):
        """Set which components to include for testing"""
        # Reset all flags
        self.include_vector_potential = False
        self.include_suck_flow = False
        self.include_swirl_flow = False
        self.include_flow_coupling = False
        
        if step_name == "baseline":
            # Working version - just refractive index
            pass  # All flags stay False
            
        elif step_name == "step1_vector":
            # Add vector potential (but don't use it in ray equations yet)
            self.include_vector_potential = True
            
        elif step_name == "step2_suck":
            # Add vector potential + suck flow
            self.include_vector_potential = True
            self.include_suck_flow = True
            
        elif step_name == "step3_swirl":
            # Add vector potential + suck flow + swirl flow
            self.include_vector_potential = True
            self.include_suck_flow = True
            self.include_swirl_flow = True
            
        elif step_name == "step4_complete":
            # Everything
            self.include_vector_potential = True
            self.include_suck_flow = True
            self.include_swirl_flow = True
            self.include_flow_coupling = True
            
        print(f"\nStep set to: {step_name}")
        print(f"  Vector potential: {self.include_vector_potential}")
        print(f"  Suck flow: {self.include_suck_flow}")
        print(f"  Swirl flow: {self.include_swirl_flow}")
        print(f"  Flow coupling: {self.include_flow_coupling}")
    
    def scalar_potential(self, r):
        """
        WORKING VERSION: Scalar potential with 1 PN corrections
        This is exactly what gave us 99.6% agreement
        """
        if np.any(r <= 0):
            return np.full_like(r, -np.inf)
        
        # Newtonian potential
        phi_0 = -G * self.mass / r
        
        # 1 PN correction
        if self.include_1pn:
            phi_2 = (G * self.mass)**2 / (2 * c**2 * r**2)
            return phi_0 + phi_2
        else:
            return phi_0
    
    def refractive_index(self, r):
        """
        WORKING VERSION: 4-fold enhanced refractive index
        This is the key formula that worked!
        """
        phi = self.scalar_potential(r)
        
        # 4-fold enhancement from 4D→3D projection (Section 2.7)
        phi_over_c2 = -4 * phi / (2 * c**2)
        
        return 1 + phi_over_c2
    
    def vector_potential(self, x, y):
        """
        STEP 1: Add vector potential from solar rotation
        A(r,θ) = (2G/c²)(J × r)/r³
        """
        if not self.include_vector_potential or self.spin == 0:
            return np.array([0.0, 0.0])
        
        r = np.sqrt(x**2 + y**2)
        if r <= 0:
            return np.array([0.0, 0.0])
        
        # For rotation about z-axis: A_φ = (2G J)/(c² r²)
        # Convert to Cartesian
        A_phi = 2 * G * self.spin / (c**2 * r**2)
        
        cos_phi = x / r
        sin_phi = y / r
        
        A_x = -A_phi * sin_phi
        A_y = A_phi * cos_phi
        
        return np.array([A_x, A_y])
    
    def suck_flow_velocity(self, x, y):
        """
        STEP 2: Add "suck" flow v_suck = -∇Φ
        """
        if not self.include_suck_flow:
            return np.array([0.0, 0.0])
        
        r = np.sqrt(x**2 + y**2)
        if r <= 1e-10:
            return np.array([0.0, 0.0])
        
        # Calculate -dΦ/dr numerically
        dr = 1e-4
        r_plus = r + dr
        r_minus = max(r - dr, 1e-10)
        
        phi_plus = self.scalar_potential(r_plus)
        phi_minus = self.scalar_potential(r_minus)
        
        dPhi_dr = (phi_plus - phi_minus) / (2 * dr)
        v_suck_r = -dPhi_dr
        
        # Convert to Cartesian
        v_suck_x = v_suck_r * x / r
        v_suck_y = v_suck_r * y / r
        
        return np.array([v_suck_x, v_suck_y])
    
    def swirl_flow_velocity(self, x, y):
        """
        STEP 3: Add "swirl" flow v_swirl = ∇×A
        """
        if not self.include_swirl_flow or not self.include_vector_potential:
            return np.array([0.0, 0.0])
        
        r = np.sqrt(x**2 + y**2)
        if r <= 1e-10:
            return np.array([0.0, 0.0])
        
        # Frame-dragging creates tangential velocity
        # For A_φ = 2GJ/(c²r²), the tangential velocity is approximately
        omega_fd = G * self.spin / (c**2 * r**3)  # Frame-dragging angular velocity
        
        # Tangential velocity components
        v_swirl_x = -omega_fd * y
        v_swirl_y = omega_fd * x
        
        return np.array([v_swirl_x, v_swirl_y])
    
    def total_flow_velocity(self, x, y):
        """
        Combine suck + swirl flows
        """
        v_suck = self.suck_flow_velocity(x, y)
        v_swirl = self.swirl_flow_velocity(x, y)
        return v_suck + v_swirl
    
    def ray_equations_incremental(self, s, state):
        """
        Ray equations using EXACT formulation from successful version
        
        State: [x, y, vx, vy] - position and velocity (not wavevector!)
        Parameter: s - path length (not time!)
        """
        x, y, vx, vy = state
        
        r = np.sqrt(x**2 + y**2)
        r = max(r, 0.1 * R_sun)  # Same safety as successful version
        
        n = self.refractive_index(r)
        
        # Calculate gradient of refractive index (same as successful version)
        dr = 0.01 * R_sun  # Same step size as successful version
        r_plus = r + dr
        r_minus = max(r - dr, 0.1 * R_sun)
        
        dn_dr = (self.refractive_index(r_plus) - 
                self.refractive_index(r_minus)) / (2 * dr)
        
        # ∇n = (dn/dr) * r̂ = (dn/dr) * (x/r, y/r)
        dn_dx = dn_dr * x / r if r > 0 else 0
        dn_dy = dn_dr * y / r if r > 0 else 0
        
        # Normalize velocity
        v_mag = np.sqrt(vx**2 + vy**2)
        if v_mag < 1e-10:
            v_mag = 1.0
        
        dx_ds = vx / v_mag
        dy_ds = vy / v_mag
        
        # === FLOW EFFECTS (for later steps) ===
        if self.include_flow_coupling:
            # STEP 4: Add flow to velocity
            v_flow = self.total_flow_velocity(x, y)
            # Modify the velocity with flow
            effective_vx = vx + v_flow[0] * v_mag / c  # Scale appropriately
            effective_vy = vy + v_flow[1] * v_mag / c
        else:
            effective_vx = vx
            effective_vy = vy
        
        # Ray equations: d(n*v)/ds = ∇n (EXACT from successful version)
        dvx_ds = dn_dx / n - effective_vx * (effective_vx * dn_dx + effective_vy * dn_dy) / (n * v_mag**2)
        dvy_ds = dn_dy / n - effective_vy * (effective_vx * dn_dx + effective_vy * dn_dy) / (n * v_mag**2)
        
        return [dx_ds, dy_ds, dvx_ds, dvy_ds]
    
    def trace_photon_incremental(self, impact_parameter, x_range=None):
        """
        Trace photon using EXACT setup from successful version
        """
        if x_range is None:
            x_range = (-10 * R_sun, 10 * R_sun)
        
        # Initial conditions (EXACT from successful version)
        x0 = x_range[0]
        y0 = impact_parameter
        
        # Initial direction (traveling in +x direction)
        vx0 = 1.0  # Normalized velocity component
        vy0 = 0.0
        
        initial_state = [x0, y0, vx0, vy0]
        
        # Path length parameter (not time!)
        total_distance = 2 * abs(x_range[1] - x_range[0])
        s_span = (0, total_distance)
        
        print(f"    Tracing with path length: {total_distance/R_sun:.1f} R_sun")
        
        try:
            solution = solve_ivp(
                self.ray_equations_incremental,
                s_span,  # Path length span
                initial_state,
                method='DOP853',
                rtol=1e-10,
                atol=1e-12,
                max_step=total_distance/1000,  # Fine steps
                dense_output=True
            )
            
            if solution.success:
                print(f"    Integration: SUCCESS")
                return solution
            else:
                print(f"    Integration FAILED: {solution.message}")
                return None
                
        except Exception as e:
            print(f"    Integration ERROR: {e}")
            return None
    
    def calculate_deflection_incremental(self, solution):
        """
        Calculate deflection using EXACT method from successful version
        """
        if solution is None or not solution.success:
            return 0.0
        
        x_vals = solution.y[0]
        y_vals = solution.y[1]
        vx_vals = solution.y[2]  # Velocity components (not wavevector!)
        vy_vals = solution.y[3]
        
        n_points = len(x_vals)
        if n_points < 100:
            return 0.0
        
        # Use asymptotic directions (EXACT from successful version)
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
        # Sign flip for correct convention (EXACT from successful version)
        deflection_rad = -deflection_raw  
        
        return deflection_rad
    
    def test_step(self, step_name, impact_parameter=1.1*R_sun):
        """
        Test a specific step and return results
        """
        print(f"\n{'='*50}")
        print(f"TESTING: {step_name.upper()}")
        print(f"{'='*50}")
        
        # Set the step
        self.set_step(step_name)
        
        # Trace photon
        print(f"Tracing photon with impact parameter {impact_parameter/R_sun:.2f} R_sun...")
        solution = self.trace_photon_incremental(impact_parameter)
        
        if solution is None or not solution.success:
            print(f"❌ Integration FAILED for {step_name}")
            return None, 0.0, 0.0
        
        # Calculate deflection
        deflection_rad = self.calculate_deflection_incremental(solution)
        deflection_arcsec = deflection_rad * 180 * 3600 / np.pi
        
        # GR prediction
        gr_prediction_rad = 4 * G * self.mass / (c**2 * impact_parameter)
        gr_prediction_arcsec = gr_prediction_rad * 180 * 3600 / np.pi
        
        # Results
        ratio = deflection_rad / gr_prediction_rad if gr_prediction_rad != 0 else 0
        
        print(f"Results for {step_name}:")
        print(f"  Deflection: {deflection_arcsec:.3f} arcsec")
        print(f"  GR prediction: {gr_prediction_arcsec:.3f} arcsec")
        print(f"  Ratio: {ratio:.4f}")
        print(f"  Agreement: {ratio*100:.2f}%")
        
        return solution, deflection_arcsec, gr_prediction_arcsec
    
    def run_incremental_test(self):
        """
        Run all steps incrementally and compare results
        """
        print("="*70)
        print("INCREMENTAL AETHER FRAMEWORK TEST")
        print("="*70)
        
        steps = [
            ("baseline", "Refractive index only (working version)"),
            ("step1_vector", "Add vector potential A(r,θ)"),
            ("step2_suck", "Add suck flow v_suck = -∇Φ"),
            ("step3_swirl", "Add swirl flow v_swirl = ∇×A"),
            ("step4_complete", "Add flow coupling to ray equations")
        ]
        
        results = []
        
        for step_name, description in steps:
            print(f"\n{'-'*50}")
            print(f"Testing: {description}")
            print(f"{'-'*50}")
            
            solution, deflection, gr_pred = self.test_step(step_name)
            
            if solution is not None:
                results.append((step_name, description, deflection, gr_pred))
                
                # Check if this step broke things
                if deflection < 0.1:  # Essentially zero
                    print(f"⚠️  WARNING: {step_name} gave near-zero deflection!")
                    print("   This step may have broken the working version.")
                    break
                elif abs(deflection - 1.584) > 0.5:  # Major deviation from working
                    print(f"⚠️  WARNING: {step_name} deviated significantly from baseline")
                    print(f"   Expected ~1.584, got {deflection:.3f}")
            else:
                print(f"❌ {step_name} failed completely")
                break
        
        # Summary
        print(f"\n{'='*70}")
        print("INCREMENTAL TEST SUMMARY")
        print(f"{'='*70}")
        
        if results:
            print(f"{'Step':<15} {'Description':<35} {'Deflection':<12} {'Agreement':<12}")
            print(f"{'-'*75}")
            
            for step_name, description, deflection, gr_pred in results:
                agreement = deflection / gr_pred * 100 if gr_pred != 0 else 0
                print(f"{step_name:<15} {description[:30]:<35} {deflection:<12.3f} {agreement:<12.2f}%")
        
        return results

def main():
    """Run incremental test"""
    framework = IncrementalAetherFramework()
    results = framework.run_incremental_test()
    
    if results:
        print(f"\n{'='*50}")
        print("ANALYSIS:")
        print(f"{'='*50}")
        
        baseline_deflection = results[0][2]  # First result should be baseline
        
        print(f"Baseline (working version): {baseline_deflection:.3f} arcsec")
        
        for i, (step_name, description, deflection, gr_pred) in enumerate(results[1:], 1):
            change = deflection - baseline_deflection
            print(f"Step {i}: {deflection:.3f} arcsec (change: {change:+.3f})")
            
            if abs(change) > 0.1:
                if change > 0:
                    print(f"   → Step {i} IMPROVED agreement")
                else:
                    print(f"   → Step {i} DEGRADED agreement")
            else:
                print(f"   → Step {i} had minimal effect")

if __name__ == "__main__":
    main()
