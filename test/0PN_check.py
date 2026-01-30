from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

# =========================
#  Mercury orbit (2-body) with Newtonian gravity
# =========================

# Physical constants (SI)
G = 6.67430e-11                 # m^3 kg^-1 s^-2
M_sun = 1.98847e30              # kg
mu = G * M_sun                  # Sun gravitational parameter, m^3/s^2

# Mercury approximate orbital elements (can tweak)
a = 5.7909227e10                # semi-major axis [m] (~0.387 AU)
e = 0.205630                    # eccentricity

# Perihelion initial condition in orbital plane (2D)
rp = a * (1 - e)                # perihelion distance
vp = np.sqrt(mu * (1 + e) / (a * (1 - e)))  # speed at perihelion

# Initial state: position on +x axis, velocity along +y axis
y0 = np.array([rp, 0.0, 0.0, vp], dtype=float)

def two_body(t, y):
    x, y_pos, vx, vy = y
    r2 = x*x + y_pos*y_pos
    r = np.sqrt(r2)
    ax = -mu * x / (r**3)
    ay = -mu * y_pos / (r**3)
    return np.array([vx, vy, ax, ay], dtype=float)

# Orbital period (Kepler 3rd law): T = 2π sqrt(a^3/mu)
T = 2 * np.pi * np.sqrt(a**3 / mu)

# Simulate N orbits
N_orbits = 2
t0, tf = 0.0, N_orbits * T
t_eval = np.linspace(t0, tf, 4000)

sol = solve_ivp(
    two_body,
    [t0, tf],
    y0,
    method="RK45",
    t_eval=t_eval,
    rtol=1e-10,
    atol=1e-12,
)

# =========================
#  Plot orbit and diagnostics
# =========================

x = sol.y[0]
y_pos = sol.y[1]
vx = sol.y[2]
vy = sol.y[3]
r = np.sqrt(x**2 + y_pos**2)
v2 = vx**2 + vy**2

# Specific orbital energy (should be constant in 2-body)
E = 0.5 * v2 - mu / r
E0 = E[0]

# Plot orbit (xy)
plt.figure()
plt.plot(0, 0, "o", markersize=8, label="Sun")
plt.plot(x, y_pos, label="Mercury (2-body)")
plt.gca().set_aspect("equal", "box")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.grid(True)
plt.legend()
plt.title("Mercury orbit around the Sun (Newtonian 2-body)")

# Plot radius vs time
plt.figure()
plt.plot(sol.t / 86400, r)
plt.xlabel("time [days]")
plt.ylabel("r [m]")
plt.grid(True)
plt.title("Distance from Sun vs time")

# Plot relative energy error
plt.figure()
plt.plot(sol.t / 86400, (E - E0) / abs(E0))
plt.xlabel("time [days]")
plt.ylabel("relative energy error")
plt.grid(True)
plt.title("Energy conservation check (should stay near 0)")

plt.show()
