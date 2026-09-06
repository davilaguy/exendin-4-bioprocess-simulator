import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Literature-informed growth parameters

mu_g_max = 0.177      # Maximum growth rate on glycerol [1/h]
mu_m_max = 0.070      # Maximum growth rate on methanol [1/h]

Y_xg = 0.40           # Biomass yield on glycerol [g biomass/g glycerol]
Y_xm = 0.55           # Biomass yield on methanol [g biomass/g methanol]

G_feed = 1260.0       # Glycerol concentration in feed [g/L]
M_feed = 792.0        # Methanol concentration in feed [g/L]

k_d = 0.01            # Biomass decay coefficient [1/h]


# Substrate kinetic parameters

K_g = 0.50            # Glycerol half-saturation constant [g/L]

K_m = 0.10            # Methanol kinetic constant [g/L]
M_inhibition_peak = 3.65   # Approximate methanol inhibition region [g/L]

K_i_m = M_inhibition_peak**2 / K_m


# Process schedule

t_glycerol_feed_start = 32.0   # Start glycerol fed-batch [h]
t_glycerol_feed_end = 46.0     # Stop glycerol feed [h]

t_methanol_start = 48.0        # Start methanol induction [h]

mu_g_target = 0.14             # Target growth rate during glycerol feeding [1/h]
mu_m_target = 0.015            # Target growth rate during methanol feeding [1/h]


# Exendin-4 production parameters

Y_p_m = 4.70          # Product formation calibration factor [mg product/g methanol]
K_g_repression = 0.10 # Glycerol repression parameter [g/L]
k_p = 0.005           # Product degradation coefficient [1/h]


# Oxygen-transfer parameters

C_star = 0.20         # Saturated dissolved oxygen concentration [mmol/L]
K_o = 0.01            # Oxygen half-saturation constant [mmol/L]

qO2_g_max = 2.20      # Max oxygen uptake on glycerol [mmol O2/(g biomass*h)]
qO2_m_max = 1.70      # Max oxygen uptake on methanol [mmol O2/(g biomass*h)]

kla = 450.0           # Volumetric oxygen-transfer coefficient [1/h]


# Create output directory if it does not already exist

os.makedirs("Images", exist_ok=True)


# Glycerol growth kinetics

def glycerol_growth(G):

    G = max(G, 0.0)

    mu_g = mu_g_max * G / (K_g + G)

    return mu_g


# Methanol growth kinetics with substrate inhibition

def methanol_growth(M):

    M = max(M, 0.0)

    if M == 0.0:
        return 0.0

    kinetic_term = M / (
        K_m
        + M
        + M**2 / K_i_m
    )

    maximum_term = M_inhibition_peak / (
        K_m
        + M_inhibition_peak
        + M_inhibition_peak**2 / K_i_m
    )

    mu_m = mu_m_max * kinetic_term / maximum_term

    return mu_m


# Glycerol feed strategy

def glycerol_feed_rate(t, X, V):

    if t < t_glycerol_feed_start:
        return 0.0

    if t >= t_glycerol_feed_end:
        return 0.0

    F_g = (
        mu_g_target
        * X
        * V
        / (Y_xg * G_feed)
    )

    return F_g


# Methanol feed strategy

def methanol_feed_rate(t, X, V):

    if t < t_methanol_start:
        return 0.0

    F_m = (
        mu_m_target
        * X
        * V
        / (Y_xm * M_feed)
    )

    return F_m


# Exendin-4 production model

def product_rate(t, G, q_m):

    if t < t_methanol_start:
        return 0.0

    G = max(G, 0.0)

    glycerol_repression = 1.0 / (
        1.0 + G / K_g_repression
    )

    q_p = (
        Y_p_m
        * q_m
        * glycerol_repression
    )

    return q_p


# Complete bioreactor model

def bioreactor_model(t, y):

    X, G, M, P, CL, V = y

    # Substrate-dependent growth rates before oxygen limitation

    mu_g_raw = glycerol_growth(G)
    mu_m_raw = methanol_growth(M)

    # Oxygen limitation

    CL_physical = max(CL, 0.0)

    oxygen_limitation = (
        CL_physical
        / (K_o + CL_physical)
    )

    mu_g = mu_g_raw * oxygen_limitation
    mu_m = mu_m_raw * oxygen_limitation

    # Feed rates

    F_g = glycerol_feed_rate(t, X, V)
    F_m = methanol_feed_rate(t, X, V)

    F_total = F_g + F_m

    # Dilution rate

    D = F_total / V

    # Specific substrate uptake rates

    q_g = mu_g / Y_xg
    q_m = mu_m / Y_xm

    # Product formation rate

    q_p = product_rate(t, G, q_m)

    # Oxygen uptake rate

    glycerol_oxygen_fraction = 0.0
    methanol_oxygen_fraction = 0.0

    if mu_g_max > 0.0:
        glycerol_oxygen_fraction = mu_g / mu_g_max

    if mu_m_max > 0.0:
        methanol_oxygen_fraction = mu_m / mu_m_max

    OUR = X * (
        qO2_g_max * glycerol_oxygen_fraction
        + qO2_m_max * methanol_oxygen_fraction
    )

    # Oxygen transfer rate

    OTR = kla * (
        C_star - CL
    )

    # Biomass balance

    dXdt = (
        mu_g
        + mu_m
        - k_d
        - D
    ) * X

    # Glycerol balance

    dGdt = (
        (F_g / V) * (G_feed - G)
        - (F_m / V) * G
        - q_g * X
    )

    # Methanol balance

    dMdt = (
        (F_m / V) * (M_feed - M)
        - (F_g / V) * M
        - q_m * X
    )

    # Exendin-4 balance

    dPdt = (
        q_p * X
        - D * P
        - k_p * P
    )

    # Dissolved oxygen balance

    dCLdt = (
        OTR
        - OUR
    )

    # Reactor volume balance

    dVdt = F_total

    return [
        dXdt,
        dGdt,
        dMdt,
        dPdt,
        dCLdt,
        dVdt
    ]


# Initial conditions

X0 = 0.10      # Biomass concentration [g/L]
G0 = 40.0      # Glycerol concentration [g/L]
M0 = 0.0       # Methanol concentration [g/L]
P0 = 0.0       # Exendin-4 concentration [mg/L]

CL0 = C_star   # Initially oxygen-saturated broth [mmol/L]

V0 = 2.0       # Initial working volume [L]

y0 = [
    X0,
    G0,
    M0,
    P0,
    CL0,
    V0
]


# Simulation settings

t_start = 0.0
t_end = 120.0

t_eval = np.linspace(
    t_start,
    t_end,
    1200
)


# Solve differential equations

solution = solve_ivp(
    bioreactor_model,
    [t_start, t_end],
    y0,
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-9,
    max_step=0.1
)


# Extract state trajectories

time = solution.t

X = solution.y[0]
G = solution.y[1]
M = solution.y[2]
P = solution.y[3]
CL = solution.y[4]
V = solution.y[5]


# Convert dissolved oxygen concentration to percent saturation

DO = (
    100.0
    * np.maximum(CL, 0.0)
    / C_star
)


# Reconstruct feed profiles

F_g_profile = np.array([
    glycerol_feed_rate(t, x, v)
    for t, x, v in zip(time, X, V)
])

F_m_profile = np.array([
    methanol_feed_rate(t, x, v)
    for t, x, v in zip(time, X, V)
])


# Reconstruct kinetic profiles

mu_g_raw_profile = np.array([
    glycerol_growth(g)
    for g in G
])

mu_m_raw_profile = np.array([
    methanol_growth(m)
    for m in M
])

oxygen_limitation_profile = (
    np.maximum(CL, 0.0)
    / (
        K_o
        + np.maximum(CL, 0.0)
    )
)

mu_g_profile = (
    mu_g_raw_profile
    * oxygen_limitation_profile
)

mu_m_profile = (
    mu_m_raw_profile
    * oxygen_limitation_profile
)


# Reconstruct oxygen uptake and transfer profiles

glycerol_oxygen_fraction_profile = (
    mu_g_profile / mu_g_max
)

methanol_oxygen_fraction_profile = (
    mu_m_profile / mu_m_max
)

OUR_profile = X * (
    qO2_g_max
    * glycerol_oxygen_fraction_profile
    +
    qO2_m_max
    * methanol_oxygen_fraction_profile
)

OTR_profile = kla * (
    C_star - CL
)


# Mass calculations

biomass_mass = X * V
glycerol_mass = G * V
methanol_mass = M * V
product_mass = P * V

glycerol_fed = np.trapezoid(
    F_g_profile * G_feed,
    time
)

methanol_fed = np.trapezoid(
    F_m_profile * M_feed,
    time
)


# Process summary calculations

minimum_DO = np.min(DO)

minimum_DO_index = np.argmin(DO)
minimum_DO_time = time[minimum_DO_index]

maximum_OUR = np.max(OUR_profile)
maximum_OTR = np.max(OTR_profile)


# Print process summary

print("\nPichia pastoris Exendin-4 process")

print("\nSolver:")
print("Successful:", solution.success)

print("\nFinal reactor state:")
print(f"Biomass: {X[-1]:.2f} g/L")
print(f"Glycerol: {G[-1]:.3f} g/L")
print(f"Methanol: {M[-1]:.3f} g/L")
print(f"Exendin-4 fusion protein: {P[-1]:.2f} mg/L")
print(f"Dissolved oxygen: {DO[-1]:.1f} %")
print(f"Volume: {V[-1]:.2f} L")

print("\nProcess totals:")
print(f"Final biomass mass: {biomass_mass[-1]:.2f} g")
print(f"Final product mass: {product_mass[-1]:.2f} mg")
print(f"Glycerol added: {glycerol_fed:.2f} g")
print(f"Methanol added: {methanol_fed:.2f} g")

print("\nMethanol:")
print(f"Maximum methanol concentration: {np.max(M):.3f} g/L")

print("\nOxygen:")
print(f"Minimum dissolved oxygen: {minimum_DO:.1f} %")
print(f"Minimum DO time: {minimum_DO_time:.1f} h")
print(f"Maximum OUR: {maximum_OUR:.2f} mmol/(L*h)")
print(f"Maximum OTR: {maximum_OTR:.2f} mmol/(L*h)")


# Plot biomass and glycerol

plt.figure()

plt.plot(
    time,
    X,
    label="Biomass"
)

plt.plot(
    time,
    G,
    label="Glycerol"
)

plt.axvline(
    t_glycerol_feed_start,
    linestyle="--",
    label="Glycerol feed start"
)

plt.axvline(
    t_methanol_start,
    linestyle=":",
    label="Methanol induction"
)

plt.xlabel("Time [h]")
plt.ylabel("Concentration [g/L]")
plt.title("Pichia pastoris Growth")

plt.legend()
plt.grid()

plt.savefig(
    "Images/pichia_growth_with_oxygen.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Plot methanol concentration

plt.figure()

plt.plot(
    time,
    M
)

plt.axvline(
    t_methanol_start,
    linestyle="--",
    label="Methanol induction"
)

plt.axhline(
    M_inhibition_peak,
    linestyle=":",
    label="Methanol inhibition region"
)

plt.xlabel("Time [h]")
plt.ylabel("Methanol [g/L]")
plt.title("Methanol During AOX1 Induction")

plt.legend()
plt.grid()

plt.savefig(
    "Images/methanol_with_oxygen.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Plot Exendin-4 production

plt.figure()

plt.plot(
    time,
    P
)

plt.axvline(
    t_methanol_start,
    linestyle="--",
    label="Methanol induction"
)

plt.xlabel("Time [h]")
plt.ylabel("Fusion Protein [mg/L]")
plt.title("Recombinant Exendin-4 Production")

plt.legend()
plt.grid()

plt.savefig(
    "Images/exendin4_with_oxygen.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Plot dissolved oxygen

plt.figure()

plt.plot(
    time,
    DO
)

plt.axhline(
    30.0,
    linestyle="--",
    label="30% DO target"
)

plt.axvline(
    t_methanol_start,
    linestyle=":",
    label="Methanol induction"
)

plt.xlabel("Time [h]")
plt.ylabel("Dissolved Oxygen [% saturation]")
plt.title("Dissolved Oxygen Without Control")

plt.legend()
plt.grid()

plt.savefig(
    "Images/uncontrolled_DO.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Plot oxygen demand versus oxygen transfer

plt.figure()

plt.plot(
    time,
    OUR_profile,
    label="OUR"
)

plt.plot(
    time,
    OTR_profile,
    label="OTR"
)

plt.xlabel("Time [h]")
plt.ylabel("Oxygen Rate [mmol/(L*h)]")
plt.title("Oxygen Uptake and Transfer")

plt.legend()
plt.grid()

plt.savefig(
    "Images/oxygen_transfer_vs_demand.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Plot feed strategy

plt.figure()

plt.plot(
    time,
    F_g_profile,
    label="Glycerol feed"
)

plt.plot(
    time,
    F_m_profile,
    label="Methanol feed"
)

plt.xlabel("Time [h]")
plt.ylabel("Feed Rate [L/h]")
plt.title("Fed-Batch Feed Strategy")

plt.legend()
plt.grid()

plt.savefig(
    "Images/feed_strategy_with_oxygen.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()