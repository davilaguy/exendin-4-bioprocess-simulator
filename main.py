import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ============================================================
# MODEL PARAMETERS
# ============================================================

mu_max = 0.40   # Maximum specific growth rate [1/h]
K_s = 0.50      # Monod half-saturation constant [g/L]
Y_xs = 0.50     # Biomass yield on substrate [g biomass / g substrate]
k_d = 0.01      # Biomass decay coefficient [1/h]

F = 0.10        # Feed flow rate [L/h]
S_f = 200.0     # Substrate concentration in feed [g/L]


# ============================================================
# KINETICS
# ============================================================

def monod_growth(S, mu_max, K_s):
    mu = mu_max * S / (K_s + S)
    return mu


# ============================================================
# BATCH MODEL
# ============================================================

def batch_model(t, y):
    X, S = y

    mu = monod_growth(S, mu_max, K_s)

    dXdt = (mu - k_d) * X
    dSdt = -(mu * X) / Y_xs

    return [dXdt, dSdt]


# ============================================================
# FED-BATCH MODEL
# ============================================================

def fed_batch_model(t, y):
    X, S, V = y

    # Current microbial growth rate
    mu = monod_growth(S, mu_max, K_s)

    # Current dilution rate
    D = F / V

    # Dynamic mass balances
    dXdt = (mu - k_d - D) * X

    dSdt = (
        D * (S_f - S)
        - (mu * X) / Y_xs
    )

    dVdt = F

    return [dXdt, dSdt, dVdt]


# ============================================================
# INITIAL CONDITIONS
# ============================================================

X0 = 0.10      # Initial biomass concentration [g/L]
S0 = 20.0      # Initial substrate concentration [g/L]
V0 = 2.0       # Initial reactor volume [L]

y0 = [X0, S0, V0]


# ============================================================
# SIMULATION SETTINGS
# ============================================================

t_start = 0.0
t_end = 30.0

t_eval = np.linspace(t_start, t_end, 500)


# ============================================================
# SOLVE DIFFERENTIAL EQUATIONS
# ============================================================

solution = solve_ivp(
    fed_batch_model,
    [t_start, t_end],
    y0,
    t_eval=t_eval
)


# ============================================================
# EXTRACT RESULTS
# ============================================================

time = solution.t

X = solution.y[0]
S = solution.y[1]
V = solution.y[2]


# ============================================================
# CALCULATE TOTAL MASSES
# ============================================================

# Concentration [g/L] * volume [L] = total mass [g]
biomass_mass = X * V
substrate_mass = S * V

initial_biomass_mass = X0 * V0
initial_substrate_mass = S0 * V0

# Substrate added through the feed
feed_substrate_mass = F * S_f * t_end

# Total substrate ever supplied to the reactor
total_substrate_supplied = (
    initial_substrate_mass
    + feed_substrate_mass
)

# Substrate consumed by the cells
substrate_consumed = (
    total_substrate_supplied
    - substrate_mass[-1]
)


# ============================================================
# VALIDATION
# ============================================================

# Because F is constant:
# V(t) = V0 + F*t
expected_final_volume = V0 + F * t_end

volume_error = (
    V[-1] - expected_final_volume
)


# ============================================================
# PRINT PROCESS SUMMARY
# ============================================================

print("\n--- FED-BATCH PROCESS SUMMARY ---")

print("Solver successful:", solution.success)

print("\nConcentrations:")
print("Final biomass concentration:", X[-1], "g/L")
print("Final substrate concentration:", S[-1], "g/L")

print("\nReactor:")
print("Final volume:", V[-1], "L")
print("Expected final volume:", expected_final_volume, "L")
print("Volume error:", volume_error, "L")

print("\nMass balance:")
print("Initial biomass mass:", initial_biomass_mass, "g")
print("Final biomass mass:", biomass_mass[-1], "g")
print("Initial substrate mass:", initial_substrate_mass, "g")
print("Feed substrate mass:", feed_substrate_mass, "g")
print("Total substrate supplied:", total_substrate_supplied, "g")
print("Final substrate mass:", substrate_mass[-1], "g")
print("Substrate consumed:", substrate_consumed, "g")


# ============================================================
# PLOT CONCENTRATIONS
# ============================================================

plt.figure()

plt.plot(time, X, label="Biomass, X")
plt.plot(time, S, label="Substrate, S")

plt.xlabel("Time [h]")
plt.ylabel("Concentration [g/L]")
plt.title("Fed-Batch Bioreactor")

plt.legend()
plt.grid()

plt.savefig(
    "Images/fed_batch_concentrations.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# PLOT REACTOR VOLUME
# ============================================================

plt.figure()

plt.plot(time, V)

plt.xlabel("Time [h]")
plt.ylabel("Reactor Volume [L]")
plt.title("Fed-Batch Reactor Volume")

plt.grid()

plt.savefig(
    "Images/fed_batch_volume.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()