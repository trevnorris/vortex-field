import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import argrelextrema
import math

# Constants
G = 6.67430e-11
c = 2.99792458e8
M = 1.9885e30  # Sun mass
mu = G * M
a = 5.7909e10  # semi-major axis
e = 0.20563    # eccentricity
period = 7.60055e6  # seconds

# Specific angular momentum h
h = math.sqrt(mu * a * (1 - e**2))

# ODE in polar: y = [r, theta, dr/dt, dtheta/dt]
def deriv(t, y):
    r, theta, dr, dtheta = y
    r2 = r**2
    r3 = r * r2
    r4 = r2 * r2

    # Newtonian
    a_r_n = -mu / r2
    cent = r * dtheta**2

    # PN correction for precession (negative sign)
    pn_corr = -3 * mu * h**2 / (c**2 * r4)

    ddr = a_r_n + cent + pn_corr
    ddtheta = -2 * dr * dtheta / r

    return [dr, dtheta, ddr, ddtheta]

# Initial at perihelion
r0 = a * (1 - e)
dr0 = 0
dtheta0 = h / r0**2
theta0 = 0
y0 = [r0, theta0, dr0, dtheta0]

# Time for 100 orbits for better accumulation
t_span = (0, 100 * period)
t_eval = np.linspace(0, t_span[1], 1000000)  # 1 million points for 100 orbits (~10k per orbit)

# Integrate
sol = solve_ivp(deriv, t_span, y0, method='RK45', t_eval=t_eval, rtol=1e-12, atol=1e-12)

# Find rough perihelia: local minima in r
r = sol.y[0]
idx = argrelextrema(r, np.less)[0]

# Filter to ensure they are perihelion (r close to r_min)
r_min_expected = a * (1 - e)
idx = idx[np.abs(r[idx] - r_min_expected) < 0.05 * a]  # Tighter filter

# Refine each min with quadratic fit
theta_peri = []
for i in idx:
    # Take 21 points around i
    s = slice(max(0, i-10), min(len(t_eval), i+11))
    t_local = t_eval[s]
    r_local = r[s]
    # Polyfit quadratic
    poly = np.polyfit(t_local, r_local, 2)
    if poly[0] > 0:  # Parabola opens up for min
        t_min = -poly[1] / (2 * poly[0])
        # Interp theta at t_min
        theta_min = np.interp(t_min, t_eval, sol.y[1])
        theta_peri.append(theta_min)

# Convert to array
theta_peri = np.array(theta_peri)

# Precession per orbit (rad)
dtheta = np.diff(theta_peri)
precession_rad = dtheta - 2 * math.pi

# Average precession per orbit in arcsec
avg_precession = np.mean(precession_rad) * (180 * 3600 / math.pi)

# Per century
orbits_century = (365.25 * 100 * 86400) / period
advance = avg_precession * orbits_century

print(f"Number of perihelia found: {len(theta_peri)}")
print(f"Average precession per orbit: {avg_precession:.6f} arcsec")
print(f"Perihelion advance: {advance:.2f} arcsec/century")
