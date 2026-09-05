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


# ============================================================
# FEED STRATEGY PARAMETERS
# ============================================================

t_feed = 12.0       # Feed start time [h]
mu_target = 0.20    # Target specific growth rate [1/h]


# ============================================================
# PRODUCT FORMATION PARAMETERS
# ============================================================

t_induction = 18.0  # Peptide-production induction time [h]

alpha = 0.20        # Growth-associated coefficient [mg product/g biomass]
beta = 0.05         # Non-growth-associated production [mg/(g biomass*h)]

k_p = 0.01          # Product degradation coefficient [1/h]


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

    # Initial batch phase
    if t < t_feed:
        return 0.0

    # Model-based feed designed around target growth rate
    F = (mu_target * X * V) / (Y_xs * S_f)

    return F


# ============================================================
# PRODUCT FORMATION
# ============================================================

def product_rate(t, mu):

    # No recombinant peptide production before induction
    if t < t_induction:
        return 0.0

    # Luedeking-Piret-type production model
    q_p = alpha * mu + beta

    return q_p


# ============================================================
# FED-BATCH BIOREACTOR MODEL
# ============================================================

def fed_batch_model(t, y):

    X, S, P, V = y

    # Microbial growth
    mu = monod_growth(S, mu_max, K_s)

    # Feed rate
    F = feed_rate(t, X, V)

    # Dilution rate
    D = F / V

    # Specific product formation rate
    q_p = product_rate(t, mu)

    # Biomass balance
    dXdt = (mu - k_d - D) * X

    # Substrate balance
    dSdt = (
        D * (S_f - S)
        - (mu * X) / Y_xs
    )

    # Therapeutic peptide balance
    dPdt = (
        q_p * X
        - D * P
        - k_p * P
    )

    # Reactor volume balance
    dVdt = F

    return [dXdt, dSdt, dPdt, dVdt]


# ============================================================
# INITIAL CONDITIONS
# ============================================================

X0 = 0.10      # Initial biomass concentration [g/L]
S0 = 20.0      # Initial substrate concentration [g/L]
P0 = 0.0       # Initial product concentration [mg/L]
V0 = 2.0       # Initial reactor volume [L]

y0 = [X0, S0, P0, V0]


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
P = solution.y[2]
V = solution.y[3]


# ============================================================
# RECONSTRUCT OPERATING PROFILES
# ============================================================

F_profile = np.array([
    feed_rate(t, x, v)
    for t, x, v in zip(time, X, V)
])

D_profile = F_profile / V


mu_profile = np.array([
    monod_growth(s, mu_max, K_s)
    for s in S
])


q_p_profile = np.array([
    product_rate(t, mu)
    for t, mu in zip(time, mu_profile)
])


# ============================================================
# MASS CALCULATIONS
# ============================================================

biomass_mass = X * V
substrate_mass = S * V

# P is mg/L and V is L, therefore product mass is mg
product_mass = P * V

initial_biomass_mass = X0 * V0
initial_substrate_mass = S0 * V0

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
# PRODUCTIVITY METRICS
# ============================================================

final_product_mass = product_mass[-1]

overall_productivity = (
    final_product_mass / t_end
)

production_time = (
    t_end - t_induction
)

production_phase_productivity = (
    final_product_mass / production_time
)


# ============================================================
# PROCESS SUMMARY
# ============================================================

print("\n--- VENOM-DERIVED PEPTIDE BIOPROCESS ---")

print("Solver successful:", solution.success)

print("\nFinal reactor state:")
print(f"Biomass concentration: {X[-1]:.2f} g/L")
print(f"Substrate concentration: {S[-1]:.3f} g/L")
print(f"Peptide concentration: {P[-1]:.2f} mg/L")
print(f"Reactor volume: {V[-1]:.2f} L")

print("\nBiomass and substrate:")
print(f"Final biomass mass: {biomass_mass[-1]:.2f} g")
print(f"Substrate supplied: {total_substrate_supplied:.2f} g")
print(f"Substrate consumed: {substrate_consumed:.2f} g")

print("\nTherapeutic peptide:")
print(f"Induction time: {t_induction:.1f} h")
print(f"Final peptide mass: {final_product_mass:.2f} mg")
print(f"Overall productivity: {overall_productivity:.2f} mg/h")
print(
    f"Production-phase productivity: "
    f"{production_phase_productivity:.2f} mg/h"
)


# ============================================================
# BIOMASS AND SUBSTRATE PLOT
# ============================================================

plt.figure()

plt.plot(time, X, label="Biomass, X")
plt.plot(time, S, label="Substrate, S")

plt.axvline(
    t_feed,
    linestyle="--",
    label="Feed start"
)

plt.axvline(
    t_induction,
    linestyle=":",
    label="Induction"
)

plt.xlabel("Time [h]")
plt.ylabel("Concentration [g/L]")
plt.title("Fed-Batch Biomass and Substrate")

plt.legend()
plt.grid()

plt.savefig(
    "Images/biomass_substrate_induction.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# THERAPEUTIC PEPTIDE PLOT
# ============================================================

plt.figure()

plt.plot(time, P)

plt.axvline(
    t_induction,
    linestyle="--",
    label="Induction"
)

plt.xlabel("Time [h]")
plt.ylabel("Peptide Concentration [mg/L]")
plt.title("Recombinant Therapeutic Peptide Production")

plt.legend()
plt.grid()

plt.savefig(
    "Images/therapeutic_peptide_production.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# FEED PROFILE
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