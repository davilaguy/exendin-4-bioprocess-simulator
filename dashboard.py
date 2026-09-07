import matplotlib.pyplot as plt
import streamlit as st

from model import (
    run_simulation,
    t_glycerol_feed_start
)


st.set_page_config(
    page_title="Exendin-4 Bioprocess Simulator",
    page_icon="🧬",
    layout="wide"
)


st.title(
    "Venom-to-Therapeutic"
)

st.subheader(
    "Exendin-4 Bioprocess Simulator"
)

st.caption(
    "Dynamic fed-batch Pichia pastoris process model with "
    "glycerol growth, methanol-induced recombinant expression, "
    "oxygen transfer, and closed-loop dissolved-oxygen control."
)


with st.sidebar:

    st.header(
        "Process Settings"
    )

    DO_setpoint = st.slider(
        "DO setpoint [%]",
        min_value=20.0,
        max_value=50.0,
        value=30.0,
        step=1.0
    )

    initial_glycerol = st.slider(
        "Initial glycerol [g/L]",
        min_value=20.0,
        max_value=60.0,
        value=40.0,
        step=2.0
    )

    methanol_start = st.slider(
        "Methanol induction time [h]",
        min_value=47.0,
        max_value=60.0,
        value=48.0,
        step=1.0
    )

    maximum_rpm = st.slider(
        "Maximum agitation [rpm]",
        min_value=800,
        max_value=1500,
        value=1200,
        step=50
    )

    st.divider()

    st.header(
        "PI Controller"
    )

    Kp = st.slider(
        "Proportional gain",
        min_value=50.0,
        max_value=250.0,
        value=150.0,
        step=10.0
    )

    Ki = st.slider(
        "Integral gain",
        min_value=0.0,
        max_value=25.0,
        value=10.0,
        step=1.0
    )

    st.caption(
        "Changing PI gains can produce poor or "
        "unstable controller performance."
    )


with st.spinner(
    "Simulating bioreactor..."
):

    results = run_simulation(
        do_setpoint=DO_setpoint,
        initial_glycerol=initial_glycerol,
        methanol_start=methanol_start,
        maximum_rpm=maximum_rpm,
        kp=Kp,
        ki=Ki
    )


summary = results["summary"]

time = results["time"]
X = results["X"]
G = results["G"]
M = results["M"]
P = results["P"]
DO = results["DO"]
V = results["V"]
rpm = results["rpm"]

OUR = results["OUR"]
OTR = results["OTR"]

F_g = results["F_g"]
F_m = results["F_m"]


st.divider()

st.subheader(
    "Final Reactor State"
)


col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "Biomass",
    f"{summary['final_biomass']:.1f} g/L"
)

col2.metric(
    "Exendin-4",
    f"{summary['final_product']:.1f} mg/L"
)

col3.metric(
    "Dissolved Oxygen",
    f"{summary['final_DO']:.1f} %"
)

col4.metric(
    "Agitation",
    f"{summary['final_rpm']:.0f} rpm"
)


col5, col6, col7, col8 = st.columns(4)

col5.metric(
    "Glycerol",
    f"{max(summary['final_glycerol'], 0):.3f} g/L"
)

col6.metric(
    "Methanol",
    f"{max(summary['final_methanol'], 0):.3f} g/L"
)

col7.metric(
    "Volume",
    f"{summary['final_volume']:.2f} L"
)

col8.metric(
    "Minimum DO",
    f"{summary['minimum_DO']:.1f} %"
)


st.divider()

st.subheader(
    "Process Dynamics"
)


left, right = st.columns(2)


with left:

    fig_growth, ax_growth = plt.subplots()

    ax_growth.plot(
        time,
        X,
        label="Biomass"
    )

    ax_growth.plot(
        time,
        G,
        label="Glycerol"
    )

    ax_growth.axvline(
        t_glycerol_feed_start,
        linestyle="--",
        label="Glycerol feed"
    )

    ax_growth.axvline(
        methanol_start,
        linestyle=":",
        label="Methanol induction"
    )

    ax_growth.set_xlabel(
        "Time [h]"
    )

    ax_growth.set_ylabel(
        "Concentration [g/L]"
    )

    ax_growth.set_title(
        "Biomass and Glycerol"
    )

    ax_growth.grid()

    ax_growth.legend()

    st.pyplot(
        fig_growth
    )

    plt.close(
        fig_growth
    )


with right:

    fig_product, ax_product = plt.subplots()

    ax_product.plot(
        time,
        P
    )

    ax_product.axvline(
        methanol_start,
        linestyle="--",
        label="Methanol induction"
    )

    ax_product.set_xlabel(
        "Time [h]"
    )

    ax_product.set_ylabel(
        "Exendin-4 Fusion Protein [mg/L]"
    )

    ax_product.set_title(
        "Recombinant Exendin-4 Production"
    )

    ax_product.grid()

    ax_product.legend()

    st.pyplot(
        fig_product
    )

    plt.close(
        fig_product
    )


left, right = st.columns(2)


with left:

    fig_DO, ax_DO = plt.subplots()

    ax_DO.plot(
        time,
        DO,
        label="Dissolved oxygen"
    )

    ax_DO.axhline(
        DO_setpoint,
        linestyle="--",
        label="DO setpoint"
    )

    ax_DO.axvline(
        methanol_start,
        linestyle=":",
        label="Methanol induction"
    )

    ax_DO.set_xlabel(
        "Time [h]"
    )

    ax_DO.set_ylabel(
        "DO [% saturation]"
    )

    ax_DO.set_title(
        "Dissolved-Oxygen Control"
    )

    ax_DO.grid()

    ax_DO.legend()

    st.pyplot(
        fig_DO
    )

    plt.close(
        fig_DO
    )


with right:

    fig_rpm, ax_rpm = plt.subplots()

    ax_rpm.plot(
        time,
        rpm
    )

    ax_rpm.axhline(
        maximum_rpm,
        linestyle="--",
        label="Maximum RPM"
    )

    ax_rpm.set_xlabel(
        "Time [h]"
    )

    ax_rpm.set_ylabel(
        "Agitation [rpm]"
    )

    ax_rpm.set_title(
        "PI Controller Agitation Response"
    )

    ax_rpm.grid()

    ax_rpm.legend()

    st.pyplot(
        fig_rpm
    )

    plt.close(
        fig_rpm
    )


st.divider()

st.subheader(
    "Engineering Diagnostics"
)


left, right = st.columns(2)


with left:

    fig_oxygen, ax_oxygen = plt.subplots()

    ax_oxygen.plot(
        time,
        OUR,
        label="OUR"
    )

    ax_oxygen.plot(
        time,
        OTR,
        label="OTR"
    )

    ax_oxygen.set_xlabel(
        "Time [h]"
    )

    ax_oxygen.set_ylabel(
        "Oxygen Rate [mmol/(L·h)]"
    )

    ax_oxygen.set_title(
        "Oxygen Demand vs Transfer"
    )

    ax_oxygen.grid()

    ax_oxygen.legend()

    st.pyplot(
        fig_oxygen
    )

    plt.close(
        fig_oxygen
    )


with right:

    fig_feed, ax_feed = plt.subplots()

    ax_feed.plot(
        time,
        F_g,
        label="Glycerol feed"
    )

    ax_feed.plot(
        time,
        F_m,
        label="Methanol feed"
    )

    ax_feed.set_xlabel(
        "Time [h]"
    )

    ax_feed.set_ylabel(
        "Feed Rate [L/h]"
    )

    ax_feed.set_title(
        "Feed Strategy"
    )

    ax_feed.grid()

    ax_feed.legend()

    st.pyplot(
        fig_feed
    )

    plt.close(
        fig_feed
    )


st.divider()

st.subheader(
    "Process Totals"
)


col1, col2, col3, col4 = st.columns(4)

col1.metric(
    "Biomass Mass",
    f"{summary['final_biomass_mass']:.1f} g"
)

col2.metric(
    "Exendin-4 Mass",
    f"{summary['final_product_mass']:.0f} mg"
)

col3.metric(
    "Glycerol Added",
    f"{summary['glycerol_added']:.1f} g"
)

col4.metric(
    "Methanol Added",
    f"{summary['methanol_added']:.1f} g"
)


st.caption(
    "This model is a literature-informed engineering "
    "simulation and includes calibrated and assumed parameters. "
    "It is not a validated pharmaceutical manufacturing model."
)