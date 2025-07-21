import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, ifft, fftfreq

# Parameters based on de Nova et al. (2020) for stable horizon
N = 2048  # High resolution
L = 100.0  # Larger box to avoid boundaries
dx = L / N
x = np.linspace(-L/2, L/2, N)
dt = 0.0002  # Balanced dt
t_max = 100.0  # Longer for steady state
steps = int(t_max / dt)

g = 1.0  # Cubic for approximation (lit uses quintic, but cubic stable here)
V0 = 0.5  # Potential height for step
k0 = 0.6  # Momentum shift for subsonic initial, left-moving (positive k0 for v negative)
rho0 = 1.0

# Initial: Uniform density with phase ramp for leftward flow
psi = np.sqrt(rho0) * np.exp(1j * k0 * x)  # v = -k0 (leftward, since phase gradient positive for negative v in convention)

# Step potential V(x) = V0 * Theta(x) for horizon at x=0
V = V0 * (x > 0).astype(float)

# Soft absorbers at both ends to damp waves/reflections
absorb_width = L/10
absorb_strength = 0.05
absorb_left = -1j * absorb_strength * (1 / (1 + np.exp((x + L/2 - absorb_width)/2.0)))  # Left
absorb_right = -1j * absorb_strength * (1 / (1 + np.exp(-(x - L/2 + absorb_width)/2.0)))  # Right
absorb = absorb_left + absorb_right

# Kinetic
k = 2 * np.pi * fftfreq(N, dx)
K = k**2 / 2

def split_step_gp(psi, dt, g, V, absorb):
    psi = np.exp(-1j * dt/2 * (g * np.abs(psi)**2 + V) + absorb * dt/2) * psi
    psi_k = fft(psi)
    psi_k = np.exp(-1j * dt * K) * psi_k
    psi = ifft(psi_k)
    psi = np.exp(-1j * dt/2 * (g * np.abs(psi)**2 + V) + absorb * dt/2) * psi
    return psi

# Simulate real time evolution
for step in range(steps):
    psi = split_step_gp(psi, dt, g, V, absorb)

# Final diagnostics
rho_final = np.abs(psi)**2
c_s_final = np.sqrt(g * rho_final)
phase_final = np.unwrap(np.angle(psi))
v_final = -np.gradient(phase_final) / dx  # Negative for leftward convention

# Horizon: Find where |v| > c_s on left (x<0, supersonic downstream)
mach = np.abs(v_final) / c_s_final
horizon_idx = np.where((x < 0) & (mach > 1))[0]
if len(horizon_idx) > 0:
    horizon_x = x[horizon_idx[-1]]  # Rightmost on left side (near step)
    print(f"Horizon approx at x={horizon_x:.2f}, where v first exceeds c_s")
else:
    print("No horizon formed")

dip = np.min(rho_final[x < 0]) if len(horizon_idx) > 0 else np.min(rho_final)
print(f"Min density near horizon: {dip:.3f} (rarefaction)")

# Plot (uncomment locally)
plt.figure()
plt.plot(x, rho_final, label='Density ρ')
plt.plot(x, np.abs(v_final), label='|v| (flow speed)')
plt.plot(x, c_s_final, label='c_s')
if len(horizon_idx) > 0:
    plt.axvline(horizon_x, color='r', ls='--', label='Horizon')
plt.xlim(-50, 50)
plt.legend()
plt.show()
