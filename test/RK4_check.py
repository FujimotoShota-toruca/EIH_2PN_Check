from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

def f(t, y):
    return -2 * y

sol = solve_ivp(f, [0, 5], [1], method='RK45', t_eval=np.linspace(0, 5, 100))

plt.plot(sol.t, sol.y[0], label='solve_ivp (RK45)')
plt.plot(sol.t, np.exp(-2 * sol.t), 'r--', label='Exact solution')
plt.legend()
plt.grid(True)
plt.show()