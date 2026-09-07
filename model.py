import numpy as np
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


# Initial reactor conditions

X0 = 0.10
M0 = 0.0
P0 = 0.0
CL0 = C_star
V0 = 2.0

integral0 = 0.0


# Simulation time

t_start = 0.0
t_end = 120.0


def glycerol_growth(G):

    G = max(G, 0.0)

    return (
        mu_g_max
        * G
        / (K_g + G)
    )


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


def methanol_feed_rate(
    t,
    X,
    V,
    methanol_start
):

    if t < methanol_start:
        return 0.0

    return (
        mu_m_target
        * X
        * V
        / (Y_xm * M_feed)
    )


def product_rate(
    t,
    G,
    q_m,
    methanol_start
):

    if t < methanol_start:
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


def kla_from_rpm(rpm):

    return (
        kla_ref
        * (
            rpm / rpm_ref
        )**kla_rpm_exponent
    )


def run_simulation(
    do_setpoint=30.0,
    initial_glycerol=40.0,
    methanol_start=48.0,
    maximum_rpm=1200.0,
    kp=150.0,
    ki=10.0
):

    rpm_base = rpm_min
    integral_limit = 1000.0

    y0 = [
        X0,
        initial_glycerol,
        M0,
        P0,
        CL0,
        V0,
        integral0
    ]

    def controller_output(
        DO,
        integral_error
    ):

        error = (
            do_setpoint - DO
        )

        rpm_command = (
            rpm_base
            + kp * error
            + ki * integral_error
        )

        rpm = np.clip(
            rpm_command,
            rpm_min,
            maximum_rpm
        )

        return (
            error,
            rpm
        )

    def bioreactor_model(t, y):

        X, G, M, P, CL, V, integral_error = y

        DO = (
            100.0
            * max(CL, 0.0)
            / C_star
        )

        error, rpm = controller_output(
            DO,
            integral_error
        )

        at_lower_limit = (
            rpm <= rpm_min
        )

        at_upper_limit = (
            rpm >= maximum_rpm
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

        mu_g_raw = glycerol_growth(G)
        mu_m_raw = methanol_growth(M)

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

        F_g = glycerol_feed_rate(
            t,
            X,
            V
        )

        F_m = methanol_feed_rate(
            t,
            X,
            V,
            methanol_start
        )

        F_total = (
            F_g + F_m
        )

        D = (
            F_total / V
        )

        q_g = (
            mu_g / Y_xg
        )

        q_m = (
            mu_m / Y_xm
        )

        q_p = product_rate(
            t,
            G,
            q_m,
            methanol_start
        )

        OUR = X * (
            qO2_g_max
            * mu_g
            / mu_g_max
            +
            qO2_m_max
            * mu_m
            / mu_m_max
        )

        kla = kla_from_rpm(
            rpm
        )

        OTR = (
            kla
            * (
                C_star - CL
            )
        )

        dXdt = (
            mu_g
            + mu_m
            - k_d
            - D
        ) * X

        dGdt = (
            (F_g / V)
            * (G_feed - G)
            -
            (F_m / V)
            * G
            -
            q_g * X
        )

        dMdt = (
            (F_m / V)
            * (M_feed - M)
            -
            (F_g / V)
            * M
            -
            q_m * X
        )

        dPdt = (
            q_p * X
            - D * P
            - k_p * P
        )

        dCLdt = (
            OTR - OUR
        )

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

    t_eval = np.linspace(
        t_start,
        t_end,
        2401
    )

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

    time = solution.t

    X = solution.y[0]
    G = solution.y[1]
    M = solution.y[2]
    P = solution.y[3]
    CL = solution.y[4]
    V = solution.y[5]

    integral_error = (
        solution.y[6]
    )

    DO = (
        100.0
        * np.maximum(
            CL,
            0.0
        )
        / C_star
    )

    error_profile = (
        do_setpoint - DO
    )

    rpm_command_profile = (
        rpm_base
        + kp * error_profile
        + ki * integral_error
    )

    rpm_profile = np.clip(
        rpm_command_profile,
        rpm_min,
        maximum_rpm
    )

    kla_profile = np.array([
        kla_from_rpm(rpm)
        for rpm in rpm_profile
    ])

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
            v,
            methanol_start
        )
        for t, x, v
        in zip(
            time,
            X,
            V
        )
    ])

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

    summary = {
        "solver_success": solution.success,
        "final_biomass": X[-1],
        "final_glycerol": G[-1],
        "final_methanol": M[-1],
        "final_product": P[-1],
        "final_DO": DO[-1],
        "final_rpm": rpm_profile[-1],
        "final_volume": V[-1],
        "minimum_DO": np.min(DO),
        "maximum_rpm": np.max(rpm_profile),
        "final_biomass_mass": biomass_mass[-1],
        "final_product_mass": product_mass[-1],
        "glycerol_added": glycerol_fed,
        "methanol_added": methanol_fed
    }

    return {
        "time": time,
        "X": X,
        "G": G,
        "M": M,
        "P": P,
        "DO": DO,
        "V": V,
        "rpm": rpm_profile,
        "kla": kla_profile,
        "OUR": OUR_profile,
        "OTR": OTR_profile,
        "F_g": F_g_profile,
        "F_m": F_m_profile,
        "summary": summary
    }