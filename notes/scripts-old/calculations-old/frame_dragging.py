import math

# Constants
G = 6.67430e-11
c = 2.99792458e8
J_earth = 5.86e33  # angular momentum
r_gpb = 7.02e6     # orbit radius

# Gravitomagnetic field magnitude (simplified for polar orbit average)
B_g = G * J_earth / (c**2 * r_gpb**3)

# Precession rate (mas/yr)
Omega_rad = 0.5 * B_g  # factor 0.5 for gyro precession
mas_per_yr = Omega_rad * (180 * 3600 * 1e3 / math.pi) * (365.25 * 86400)  # rad/s to mas/yr

print(f"Frame-dragging precession: {mas_per_yr:.1f} mas/yr")
