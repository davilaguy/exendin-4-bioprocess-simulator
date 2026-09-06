import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Literature-informed growth parameters

mu_g_max = 0.177      # Maximum growth rate on glycerol [1/h]
mu_m_max = 0.070      # Representative Mut+ maximum growth rate on methanol [1/h]

Y_xg = 0.40           # Biomass yield on glycerol [g biomass/g glycerol]
Y_xm = 0.55           # Biomass yield on methanol [g biomass/g methanol]

G_feed = 1260.0       # Glycerol feed concentration [g/L]
M_feed = 792.0        # Approximate pure methanol concentration [g/L]

k_d = 0.01            # Biomass decay coefficient [1/h]


# Provisional kinetic parameters

K_g = 0.50            # Glycerol half-saturation constant [g/L]

K_m = 0.10            # Methanol kinetic constant [g/L]
M_inhibition_peak = 3.65   # Methanol concentration where inhibition becomes important [g/L]

K_i_m = M_inhibition_peak**2 / K_m

# Oxygen transfer

C_star = 0.20          # Oxygen saturation concentration [mmol/L]

qO2_g = 2.2            # Oxygen demand on glycerol [mmol O2/(g*h)]
qO2_m = 1.7            # Oxygen demand on methanol [mmol O2/(g*h)]

kla = 450.0            # Volumetric oxygen-transfer coefficient [1/h]
K_o = 0.01

# Process schedule

t_glycerol_feed_start = 32.0
t_glycerol_feed_end = 46.0

t_methanol_start = 48.0

mu_g_target = 0.14
mu_m_target = 0.015


# Product parameters

Y_p_m = 4.70          # Product yield factor [mg product/g methanol]
K_g_repression = 0.10 # Glycerol repression parameter [g/L]
k_p = 0.005           # Product degradation coefficient [1/h]


# Growth kinetics

def glycerol_growth(G):

    G = max(G, 0.0)

    mu_g = mu_g_max * G / (K_g + G)

    return mu_g


def methanol_growth(M):

    M = max(M, 0.0)

    if M == 0:
        return 0.0

    kinetic_term = M / (
        K_m
        + M
        + M**2 / K_i_m
    )

    max_kinetic_term = M_inhibition_peak / (
        K_m
        + M_inhibition_peak
        + M_inhibition_peak**2 / K_i_m
    )

    mu_m = mu_m_max * kinetic_term / max_kinetic_term

    return mu_m


# Feed strategies

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


# Exendin-4 production

def product_rate(t, G, mu_m):

    if t < t_methanol_start:
        return 0.0

    G = max(G, 0.0)

    q_m = mu_m / Y_xm

    glycerol_repression = 1.0 / (
        1.0 + G / K_g_repression
    )

    q_p = (
        Y_p_m
        * q_m
        * glycerol_repression
    )

    return q_p


# Bioreactor model

def bioreactor_model(t, y):

    X, G, M, P, CL, V = y

    mu_g = glycerol_growth(G)
    mu_m = methanol_growth(M)

    F_g = glycerol_feed_rate(t, X, V)
    F_m = methanol_feed_rate(t, X, V)

    F_total = F_g + F_m
    D = F_total / V

    q_g = mu_g / Y_xg
    q_m = mu_m / Y_xm

    q_p = product_rate(t, G, mu_m)

    OUR = X * (
    qO2_g * mu_g / mu_g_max
    + qO2_m * mu_m / mu_m_max
    )

    OTR = kla * (C_star - CL)
    dCLdt = OTR - OUR

    dXdt = (
        mu_g
        + mu_m
        - k_d
        - D
    ) * X

    dGdt = (
        (F_g / V) * (G_feed - G)
        - (F_m / V) * G
        - q_g * X
    )

    dMdt = (
        (F_m / V) * (M_feed - M)
        - (F_g / V) * M
        - q_m * X
    )

    dPdt = (
        q_p * X
        - D * P
        - k_p * P
    )

    dVdt = F_total

    return [
        dXdt,
        dGdt,
        dMdt,
        dPdt,
        dVdt
    ]


# Initial conditions

X0 = 0.10      # Biomass [g/L]
G0 = 40.0      # Glycerol [g/L]
M0 = 0.0       # Methanol [g/L]
P0 = 0.0       # Exendin-4 fusion protein [mg/L]
V0 = 2.0       # Working volume [L]
CL0 = C_star

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


# Solve model

solution = solve_ivp(
    bioreactor_model,
    [t_start, t_end],
    y0,
    t_eval=t_eval,
    rtol=1e-7,
    atol=1e-9,
    max_step=0.1
)


# Extract results

time = solution.t

X = solution.y[0]
G = solution.y[1]
M = solution.y[2]
P = solution.y[3]
V = solution.y[4]


# Reconstruct feed and kinetic profiles

F_g_profile = np.array([
    glycerol_feed_rate(t, x, v)
    for t, x, v in zip(time, X, V)
])

F_m_profile = np.array([
    methanol_feed_rate(t, x, v)
    for t, x, v in zip(time, X, V)
])

mu_g_profile = np.array([
    glycerol_growth(g)
    for g in G
])

mu_m_profile = np.array([
    methanol_growth(m)
    for m in M
])


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


# Process summary

print("\nPichia pastoris Exendin-4 process")

print("\nSolver:")
print("Successful:", solution.success)

print("\nFinal reactor state:")
print(f"Biomass: {X[-1]:.2f} g/L")
print(f"Glycerol: {G[-1]:.3f} g/L")
print(f"Methanol: {M[-1]:.3f} g/L")
print(f"Exendin-4 fusion protein: {P[-1]:.2f} mg/L")
print(f"Volume: {V[-1]:.2f} L")

print("\nProcess totals:")
print(f"Final biomass mass: {biomass_mass[-1]:.2f} g")
print(f"Final product mass: {product_mass[-1]:.2f} mg")
print(f"Glycerol added: {glycerol_fed:.2f} g")
print(f"Methanol added: {methanol_fed:.2f} g")

print("\nMethanol:")
print(f"Maximum methanol concentration: {np.max(M):.3f} g/L")


# Biomass and glycerol

plt.figure()

plt.plot(time, X, label="Biomass")
plt.plot(time, G, label="Glycerol")

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
    "Images/pichia_growth.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Methanol concentration

plt.figure()

plt.plot(time, M)

plt.axvline(
    t_methanol_start,
    linestyle="--",
    label="Methanol induction"
)

plt.axhline(
    M_inhibition_peak,
    linestyle=":",
    label="Inhibition threshold"
)

plt.xlabel("Time [h]")
plt.ylabel("Methanol [g/L]")
plt.title("Methanol During AOX1 Induction")

plt.legend()
plt.grid()

plt.savefig(
    "Images/methanol_profile.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Exendin-4 production

plt.figure()

plt.plot(time, P)

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
    "Images/exendin4_production.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Feed rates

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
    "Images/feed_strategy.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()