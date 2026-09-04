import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Model parameters
mu_max = 0.40   # Maximum specific growth rate [1/h]
K_s = 0.50      # Monod half-saturation constant [g/L]
Y_xs = 0.50     # Biomass yield on substrate [g biomass / g substrate]
k_d = 0.01      # Biomass decay coefficient [1/h]


def monod_growth(S, mu_max, K_s):
    mu = mu_max * S / (K_s + S)
    return mu


def batch_model(t, y):
    X, S = y

    mu = monod_growth(S, mu_max, K_s)

    dXdt = (mu - k_d) * X
    dSdt = -(mu * X) / Y_xs

    return [dXdt, dSdt]


# Initial conditions
X0 = 0.10
S0 = 20.0

y0 = [X0, S0]


# Simulation time
t_start = 0.0
t_end = 30.0

t_eval = np.linspace(t_start, t_end, 500)


# Solve the differential equations
solution = solve_ivp(
    batch_model,
    [t_start, t_end],
    y0,
    t_eval=t_eval
)


# Extract results
time = solution.t
X = solution.y[0]
S = solution.y[1]


# Plot results
plt.plot(time, X, label="Biomass, X")
plt.plot(time, S, label="Substrate, S")

plt.xlabel("Time [h]")
plt.ylabel("Concentration [g/L]")
plt.title("Batch Bioreactor")

plt.legend()
plt.grid()
plt.show()