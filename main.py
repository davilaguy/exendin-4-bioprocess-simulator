import os

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Growth parameters

mu_g_max = 0.177
mu_m_max = 0.070

Y_xg = 0.40
Y_xm = 0.55

G_feed = 1260.0
M_feed = 792.0

k_d = 0.01


# Substrate kinetics

K_g = 0.50

K_m = 0.10
M_inhibition_peak = 3.65

K_i_m = M_inhibition_peak**2 / K_m


# Process timing

t_glycerol_feed_start = 32.0
t_glycerol_feed_end = 46.0
t_methanol_start = 48.0

mu_g_target = 0.14
mu_m_target = 0.015


# Exendin-4 production

Y_p_m = 4.70
K_g_repression = 0.10
k_p = 0.005


# Oxygen model

C_star = 0.20
K_o = 0.01

qO2_g_max = 2.20
qO2_m_max = 1.70


# Agitation and oxygen transfer

rpm_ref = 600.0
kla_ref = 450.0

kla_rpm_exponent = 1.20

rpm_min = 300.0
rpm_max = 1200.0


# PI dissolved-oxygen controller

DO_setpoint = 30.0

rpm_base = 300.0

Kp = 150.0
Ki = 10.0

integral_limit = 1000.0


# Initial conditions

X0 = 0.10
G0 = 40.0
M0 = 0.0
P0 = 0.0
CL0 = C_star
V0 = 2.0

integral0 = 0.0


# Simulation settings

t_start = 0.0
t_end = 120.0

t_eval = np.linspace(
    t_start,
    t_end,
    2401
)


# Output folder

os.makedirs(
    "Images",
    exist_ok=True
)


# Glycerol growth

def glycerol_growth(G):

    G = max(G, 0.0)

    return (
        mu_g_max
        * G
        / (K_g + G)
    )


# Methanol growth

def methanol_growth(M):

    M = max(M, 0.0)

    if M == 0.0:
        return 0.0

    kinetic_term = (
        M
        / (
            K_m
            + M
            + M**2 / K_i_m
        )
    )

    maximum_term = (
        M_inhibition_peak
        / (
            K_m
            + M_inhibition_peak
            + M_inhibition_peak**2 / K_i_m
        )
    )

    return (
        mu_m_max
        * kinetic_term
        / maximum_term
    )


# Glycerol feed

def glycerol_feed_rate(t, X, V):

    if t < t_glycerol_feed_start:
        return 0.0

    if t >= t_glycerol_feed_end:
        return 0.0

    return (
        mu_g_target
        * X
        * V
        / (Y_xg * G_feed)
    )


# Methanol feed

def methanol_feed_rate(t, X, V):

    if t < t_methanol_start:
        return 0.0

    return (
        mu_m_target
        * X
        * V
        / (Y_xm * M_feed)
    )


# Exendin-4 formation

def product_rate(t, G, q_m):

    if t < t_methanol_start:
        return 0.0

    G = max(G, 0.0)

    glycerol_repression = (
        1.0
        / (
            1.0
            + G / K_g_repression
        )
    )

    return (
        Y_p_m
        * q_m
        * glycerol_repression
    )


# Agitation to kLa relationship

def kla_from_rpm(rpm):

    return (
        kla_ref
        * (
            rpm / rpm_ref
        )**kla_rpm_exponent
    )


# PI controller

def controller_output(DO, integral_error):

    error = (
        DO_setpoint - DO
    )

    rpm_command = (
        rpm_base
        + Kp * error
        + Ki * integral_error
    )

    rpm = np.clip(
        rpm_command,
        rpm_min,
        rpm_max
    )

    return (
        error,
        rpm,
        rpm_command
    )


# Complete reactor model

def bioreactor_model(t, y):

    X, G, M, P, CL, V, integral_error = y

    # Current dissolved oxygen

    DO = (
        100.0
        * max(CL, 0.0)
        / C_star
    )

    # PI control

    error, rpm, rpm_command = controller_output(
        DO,
        integral_error
    )

    # Anti-windup

    at_lower_limit = (
        rpm <= rpm_min
    )

    at_upper_limit = (
        rpm >= rpm_max
    )

    pushing_lower = (
        error < 0.0
    )

    pushing_upper = (
        error > 0.0
    )

    if (
        at_lower_limit
        and pushing_lower
    ):
        dintegraldt = 0.0

    elif (
        at_upper_limit
        and pushing_upper
    ):
        dintegraldt = 0.0

    elif (
        integral_error >= integral_limit
        and error > 0.0
    ):
        dintegraldt = 0.0

    elif (
        integral_error <= -integral_limit
        and error < 0.0
    ):
        dintegraldt = 0.0

    else:
        dintegraldt = error

    # Substrate-dependent growth

    mu_g_raw = glycerol_growth(G)
    mu_m_raw = methanol_growth(M)

    # Oxygen limitation

    CL_physical = max(
        CL,
        0.0
    )

    oxygen_limitation = (
        CL_physical
        / (
            K_o
            + CL_physical
        )
    )

    mu_g = (
        mu_g_raw
        * oxygen_limitation
    )

    mu_m = (
        mu_m_raw
        * oxygen_limitation
    )

    # Feed rates

    F_g = glycerol_feed_rate(
        t,
        X,
        V
    )

    F_m = methanol_feed_rate(
        t,
        X,
        V
    )

    F_total = (
        F_g + F_m
    )

    D = (
        F_total / V
    )

    # Substrate uptake

    q_g = (
        mu_g / Y_xg
    )

    q_m = (
        mu_m / Y_xm
    )

    # Product formation

    q_p = product_rate(
        t,
        G,
        q_m
    )

    # Oxygen uptake

    OUR = X * (
        qO2_g_max
        * mu_g
        / mu_g_max
        +
        qO2_m_max
        * mu_m
        / mu_m_max
    )

    # Oxygen transfer

    kla = kla_from_rpm(
        rpm
    )

    OTR = (
        kla
        * (
            C_star - CL
        )
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
        (F_g / V)
        * (G_feed - G)
        -
        (F_m / V)
        * G
        -
        q_g * X
    )

    # Methanol balance

    dMdt = (
        (F_m / V)
        * (M_feed - M)
        -
        (F_g / V)
        * M
        -
        q_m * X
    )

    # Exendin-4 balance

    dPdt = (
        q_p * X
        - D * P
        - k_p * P
    )

    # Dissolved oxygen balance

    dCLdt = (
        OTR - OUR
    )

    # Volume balance

    dVdt = (
        F_total
    )

    return [
        dXdt,
        dGdt,
        dMdt,
        dPdt,
        dCLdt,
        dVdt,
        dintegraldt
    ]


# Initial state

y0 = [
    X0,
    G0,
    M0,
    P0,
    CL0,
    V0,
    integral0
]


# Solve the process

solution = solve_ivp(
    bioreactor_model,
    [t_start, t_end],
    y0,
    t_eval=t_eval,
    method="BDF",
    rtol=1e-7,
    atol=1e-9,
    max_step=0.05
)


# Extract states

time = solution.t

X = solution.y[0]
G = solution.y[1]
M = solution.y[2]
P = solution.y[3]
CL = solution.y[4]
V = solution.y[5]

integral_error = solution.y[6]


# Dissolved oxygen

DO = (
    100.0
    * np.maximum(
        CL,
        0.0
    )
    / C_star
)


# Reconstruct controller output

error_profile = (
    DO_setpoint - DO
)

rpm_command_profile = (
    rpm_base
    + Kp * error_profile
    + Ki * integral_error
)

rpm_profile = np.clip(
    rpm_command_profile,
    rpm_min,
    rpm_max
)

kla_profile = np.array([
    kla_from_rpm(rpm)
    for rpm in rpm_profile
])


# Reconstruct feed profiles

F_g_profile = np.array([
    glycerol_feed_rate(
        t,
        x,
        v
    )
    for t, x, v
    in zip(
        time,
        X,
        V
    )
])

F_m_profile = np.array([
    methanol_feed_rate(
        t,
        x,
        v
    )
    for t, x, v
    in zip(
        time,
        X,
        V
    )
])


# Reconstruct growth rates

mu_g_raw_profile = np.array([
    glycerol_growth(g)
    for g in G
])

mu_m_raw_profile = np.array([
    methanol_growth(m)
    for m in M
])

oxygen_limitation_profile = (
    np.maximum(
        CL,
        0.0
    )
    / (
        K_o
        + np.maximum(
            CL,
            0.0
        )
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


# Reconstruct oxygen rates

OUR_profile = X * (
    qO2_g_max
    * mu_g_profile
    / mu_g_max
    +
    qO2_m_max
    * mu_m_profile
    / mu_m_max
)

OTR_profile = (
    kla_profile
    * (
        C_star - CL
    )
)


# Mass calculations

biomass_mass = (
    X * V
)

product_mass = (
    P * V
)

glycerol_fed = np.trapezoid(
    F_g_profile * G_feed,
    time
)

methanol_fed = np.trapezoid(
    F_m_profile * M_feed,
    time
)


# Process summary

print(
    "\nPichia pastoris Exendin-4 "
    "PI-controlled process"
)

print(
    "\nSolver successful:",
    solution.success
)

print("\nFinal reactor state:")

print(
    f"Biomass: "
    f"{X[-1]:.2f} g/L"
)

print(
    f"Glycerol: "
    f"{G[-1]:.3f} g/L"
)

print(
    f"Methanol: "
    f"{M[-1]:.3f} g/L"
)

print(
    f"Exendin-4 fusion protein: "
    f"{P[-1]:.2f} mg/L"
)

print(
    f"Dissolved oxygen: "
    f"{DO[-1]:.1f} %"
)

print(
    f"Agitation: "
    f"{rpm_profile[-1]:.0f} rpm"
)

print(
    f"Volume: "
    f"{V[-1]:.2f} L"
)


print("\nController performance:")

print(
    f"Minimum DO: "
    f"{np.min(DO):.1f} %"
)

print(
    f"Maximum RPM: "
    f"{np.max(rpm_profile):.0f} rpm"
)

print(
    f"Minimum RPM: "
    f"{np.min(rpm_profile):.0f} rpm"
)


print("\nProcess totals:")

print(
    f"Final biomass mass: "
    f"{biomass_mass[-1]:.2f} g"
)

print(
    f"Final product mass: "
    f"{product_mass[-1]:.2f} mg"
)

print(
    f"Glycerol added: "
    f"{glycerol_fed:.2f} g"
)

print(
    f"Methanol added: "
    f"{methanol_fed:.2f} g"
)


# Dissolved oxygen control

plt.figure()

plt.plot(
    time,
    DO,
    label="Dissolved oxygen"
)

plt.axhline(
    DO_setpoint,
    linestyle="--",
    label="DO setpoint"
)

plt.axvline(
    t_methanol_start,
    linestyle=":",
    label="Methanol induction"
)

plt.xlabel("Time [h]")
plt.ylabel("DO [% saturation]")
plt.title("PI Dissolved-Oxygen Control")

plt.legend()
plt.grid()

plt.savefig(
    "Images/PI_DO_control.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Agitation response

plt.figure()

plt.plot(
    time,
    rpm_profile
)

plt.axhline(
    rpm_max,
    linestyle="--",
    label="Maximum RPM"
)

plt.xlabel("Time [h]")
plt.ylabel("Agitation [rpm]")
plt.title("PI Agitation Response")

plt.legend()
plt.grid()

plt.savefig(
    "Images/PI_agitation.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Biomass and glycerol

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
    "Images/PI_growth.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# Exendin-4 production

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
    "Images/PI_exendin4.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# OUR and OTR

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
plt.title("Oxygen Demand and Transfer")

plt.legend()
plt.grid()

plt.savefig(
    "Images/PI_oxygen_rates.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()