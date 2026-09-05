import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# ============================================================
# MODEL PARAMETERS
# ============================================================

mu_max = 0.40   # Maximum specific growth rate [1/h]
K_s = 0.50      # Monod half-saturation constant [g/L]
Y_xs = 0.50     # Biomass yield on substrate [g biomass/g substrate]
k_d = 0.01      # Biomass decay coefficient [1/h]

S_f = 200.0     # Substrate concentration in feed [g/L]

# Feed strategy
t_feed = 12.0       # Time at which feeding begins [h]
mu_target = 0.20    # Desired specific growth rate during feeding [1/h]


# ============================================================
# KINETICS
# ============================================================

def monod_growth(S, mu_max, K_s):
    mu = mu_max * S / (K_s + S)
    return mu


# ============================================================
# FEED STRATEGY
# ============================================================

def feed_rate(t, X, V):

    # Batch phase: no feed
    if t < t_feed:
        return 0.0

    # Model-based feed required to support mu_target
    F = (mu_target * X * V) / (Y_xs * S_f)

    return F


# ============================================================
# FED-BATCH MODEL
# ============================================================

def fed_batch_model(t, y):
    X, S, V = y

    # Specific microbial growth rate
    mu = monod_growth(S, mu_max, K_s)

    # Current feed rate
    F = feed_rate(t, X, V)

    # Current dilution rate
    D = F / V

    # Dynamic balances
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
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-9
)


# ============================================================
# EXTRACT RESULTS
# ============================================================

time = solution.t

X = solution.y[0]
S = solution.y[1]
V = solution.y[2]


# ============================================================
# RECONSTRUCT FEED AND DILUTION PROFILES
# ============================================================

F_profile = np.array([
    feed_rate(t, x, v)
    for t, x, v in zip(time, X, V)
])

D_profile = F_profile / V


# ============================================================
# CALCULATE TOTAL MASSES
# ============================================================

biomass_mass = X * V
substrate_mass = S * V

initial_biomass_mass = X0 * V0
initial_substrate_mass = S0 * V0

# Feed is now time-varying, so integrate F(t)*Sf over time
feed_substrate_mass = np.trapezoid(
    F_profile * S_f,
    time
)

total_substrate_supplied = (
    initial_substrate_mass
    + feed_substrate_mass
)

substrate_consumed = (
    total_substrate_supplied
    - substrate_mass[-1]
)


# ============================================================
# PROCESS SUMMARY
# ============================================================

print("\n--- FED-BATCH PROCESS SUMMARY ---")

print("Solver successful:", solution.success)

print("\nConcentrations:")
print(f"Final biomass concentration: {X[-1]:.2f} g/L")
print(f"Final substrate concentration: {S[-1]:.3f} g/L")

print("\nReactor:")
print(f"Final volume: {V[-1]:.2f} L")
print(f"Final feed rate: {F_profile[-1]:.3f} L/h")

print("\nMass balance:")
print(f"Initial biomass mass: {initial_biomass_mass:.2f} g")
print(f"Final biomass mass: {biomass_mass[-1]:.2f} g")
print(f"Initial substrate mass: {initial_substrate_mass:.2f} g")
print(f"Feed substrate mass: {feed_substrate_mass:.2f} g")
print(f"Total substrate supplied: {total_substrate_supplied:.2f} g")
print(f"Final substrate mass: {substrate_mass[-1]:.2f} g")
print(f"Substrate consumed: {substrate_consumed:.2f} g")


# ============================================================
# PLOT CONCENTRATIONS
# ============================================================

plt.figure()

plt.plot(time, X, label="Biomass, X")
plt.plot(time, S, label="Substrate, S")

plt.axvline(
    t_feed,
    linestyle="--",
    label="Feed start"
)

plt.xlabel("Time [h]")
plt.ylabel("Concentration [g/L]")
plt.title("Model-Based Fed-Batch Bioreactor")

plt.legend()
plt.grid()

plt.savefig(
    "Images/model_based_concentrations.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# PLOT FEED RATE
# ============================================================

plt.figure()

plt.plot(time, F_profile)

plt.axvline(
    t_feed,
    linestyle="--"
)

plt.xlabel("Time [h]")
plt.ylabel("Feed Rate [L/h]")
plt.title("Model-Based Feed Profile")

plt.grid()

plt.savefig(
    "Images/model_based_feed.png",
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
plt.title("Reactor Volume")

plt.grid()

plt.savefig(
    "Images/model_based_volume.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()