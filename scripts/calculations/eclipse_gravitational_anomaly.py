import math
from sympy import symbols, series, simplify, sqrt

# Step 1: Symbolic verification with SymPy
# x = R/d (radius to distance ratio)
G, M, d, x = symbols('G M d x', positive=True)

g_point = G * M / d**2
g_disk = 2 * G * M / (x**2 * d**2) * (1 - 1 / sqrt(1 + x**2))
delta_g_sym = g_point - g_disk
expanded = series(delta_g_sym, x, 0, 5)
print("Symbolic expansion of Delta g:")
print(simplify(expanded))
print("\n")  # Line break for readability

# Step 2: Numerical calculation
# Physical constants
G_num = 6.6743e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M_num = 1.9885e30   # Mass of the Sun (kg)
d_num = 1.496e11    # Average Earth-Sun distance (m)
R_num = 6.957e8     # Radius of the Sun (m)

# Surface mass density for disk model (kg/m^2)
sigma = M_num / (math.pi * R_num**2)

# Point source gravity (m/s^2)
g_point_num = G_num * M_num / d_num**2

# Disk source gravity on axis (m/s^2)
g_disk_num = 2 * math.pi * G_num * sigma * (1 - d_num / math.sqrt(d_num**2 + R_num**2))

# Anomaly as absolute difference (m/s^2)
delta_g_num = abs(g_disk_num - g_point_num)

# Output numerical results
print("Delta g (m/s^2):", delta_g_num)
print("Delta g (μGal):", delta_g_num * 1e8)  # Conversion: 1 μGal = 10^{-8} m/s^2
