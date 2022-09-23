"""
Conceptual Aircraft Design Tool
(for PRJ-22 and AP-701 courses)

Cap. Eng. Ney Rafael Secco (ney@ita.br)
Aircraft Design Department
Aeronautics Institute of Technology

08-2022

The code uses several historical regression from
aircraft design books to make a quick initial
sizing procedure.

Generally, the user should call only the 'analyze'
function from this module.
"""

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# CONSTANTS
ft2m = 0.3048
kt2ms = 0.514444
lb2N = 4.44822
nm2m = 1852.0
gravity = 9.81
gamma_air = 1.4
R_air = 287

# ========================================
# MAIN FUNCTION


def analyze(
    airplane=None,
    print_log=False,  # Plot results on the terminal screen
    plot=False,  # Generate 3D plot of the aircraft
    W0_guess=None,  # Guess for MTOW [N]
    T0_guess=None,  # Guess for Takeoff total thrust [N]
):
    """
    This is the main function that should be used for aircraft analysis.
    """

    # Load standard airplane if none is provided
    if airplane is None:
        airplane = standard_airplane()

    # Use an average wing loading for transports
    # to estime W0_guess and T0_guess if none are provided
    if W0_guess is None:
        W0_guess = 5e3 * airplane["S_w"]
    if T0_guess is None:
        T0_guess = 0.3 * W0_guess

    ### ADD CODE FROM SECTION 4.3 HERE ###

    # Generate geometry
    geometry(airplane)

    # Converge MTOW and Takeoff Thrust
    thrust_matching(W0_guess, T0_guess, airplane)

    # Balance analysis
    balance(airplane)

    # Landing gear design
    landing_gear(airplane)

    if print_log:
        print("We [kg]", airplane["We"] / gravity)
        print("Wf [kg]", airplane["Wf"] / gravity)
        print("W0 [kg]", airplane["W0"] / gravity)
        print("T0 [kg]", airplane["T0"] / gravity)
        print("T0/W0", airplane["T0"] / airplane["W0"])
        print("W0/S", airplane["W0"] / airplane["S_w"])
        print("deltaS_wlan", airplane["deltaS_wlan"])
        print("xcg_fwd", airplane["xcg_fwd"])
        print("xcg_aft", airplane["xcg_aft"])
        print("xnp", airplane["xnp"])
        print("SM_fwd", airplane["SM_fwd"])
        print("SM_aft", airplane["SM_aft"])
        print("b_tank_b_w", airplane["b_tank_b_w"])
        print("CLv", airplane["CLv"])

        if airplane["frac_nlg_fwd"] is not None:
            print("frac_nlg_fwd", airplane["frac_nlg_fwd"])
            print("frac_nlg_aft", airplane["frac_nlg_aft"])
            print("alpha_tipback [deg]", airplane["alpha_tipback"] * 180.0 / np.pi)
            print(
                "alpha_tailstrike [deg]", airplane["alpha_tailstrike"] * 180.0 / np.pi
            )
            print("phi_overturn [deg]", airplane["phi_overturn"] * 180.0 / np.pi)

    if plot:
        plot3d(airplane)

    return airplane


# ========================================
# DISCIPLINE MODULES


def geometry(airplane):

    # Unpack dictionary
    S_w = airplane["S_w"]
    AR_w = airplane["AR_w"]
    taper_w = airplane["taper_w"]
    sweep_w = airplane["sweep_w"]
    dihedral_w = airplane["dihedral_w"]
    xr_w = airplane["xr_w"]
    zr_w = airplane["zr_w"]
    Cht = airplane["Cht"]
    AR_h = airplane["AR_h"]
    taper_h = airplane["taper_h"]
    sweep_h = airplane["sweep_h"]
    dihedral_h = airplane["dihedral_h"]
    Lc_h = airplane["Lc_h"]
    zr_h = airplane["zr_h"]
    Cvt = airplane["Cvt"]
    AR_v = airplane["AR_v"]
    taper_v = airplane["taper_v"]
    sweep_v = airplane["sweep_v"]
    Lb_v = airplane["Lb_v"]
    zr_v = airplane["zr_v"]

    ### ADD CODE FROM SECTION 3.1 HERE ###

    # WING

    b_w = np.sqrt(AR_w * S_w)
    cr_w = (2 * S_w) / (b_w * (1 + taper_w))
    ct_w = taper_w * cr_w
    yt_w = b_w / 2
    xt_w = xr_w + yt_w * np.tan(sweep_w) + (cr_w - ct_w) / 4
    zt_w = zr_w + yt_w * np.tan(dihedral_w)
    cm_w = (2 * cr_w / 3) * ((1 + taper_w + taper_w**2) / (1 + taper_w))
    ym_w = (b_w / 6) * ((1 + 2 * taper_w) / (1 + taper_w))
    xm_w = xr_w + ym_w * np.tan(sweep_w) + (cr_w - cm_w) / 4
    zm_w = zr_w + ym_w * np.tan(dihedral_w)

    # HORIZONTAL TAIL

    L_h = Lc_h * cm_w
    S_h = (S_w * cm_w / L_h) * Cht
    b_h = np.sqrt(AR_h * S_h)
    cr_h = (2 * S_h) / (b_h * (1 + taper_h))
    ct_h = taper_h * cr_h
    cm_h = (2 * cr_h / 3) * (1 + taper_h + taper_h**2) / (1 + taper_h)
    xm_h = xm_w + L_h + (cm_w - cm_h) / 4
    ym_h = (b_h / 6) * (1 + 2 * taper_h) / (1 + taper_h)
    zm_h = zr_h + ym_h * np.tan(dihedral_h)
    xr_h = xm_h - ym_h * np.tan(sweep_h) + (cm_h - cr_h) / 4
    yt_h = b_h / 2
    xt_h = xr_h + yt_h * np.tan(sweep_h) + (cr_h - ct_h) / 4
    zt_h = zr_h + yt_h * np.tan(dihedral_h)

    # VERTICAL TAIL

    L_v = Lb_v * b_w
    S_v = (S_w * b_w / L_v) * Cvt
    b_v = np.sqrt(AR_v * S_v)
    cr_v = (2 * S_v) / (b_v * (1 + taper_v))
    ct_v = taper_v * cr_v
    cm_v = (2 * cr_v / 3) * (1 + taper_v + taper_v**2) / (1 + taper_v)
    xm_v = xm_w + L_v + (cm_w - cm_v) / 4
    zm_v = zr_v + (b_v / 3) * (1 + 2 * taper_v) / (1 + taper_v)
    xr_v = xm_v - (zm_v - zr_v) * np.tan(sweep_v) + (cm_v - cr_v) / 4
    zt_v = zr_v + b_v
    xt_v = xr_v + (zt_v - zr_v) * np.tan(sweep_v) + (cr_v - ct_v) / 4

    # Update dictionary with new results
    airplane["b_w"] = b_w
    airplane["cr_w"] = cr_w
    airplane["xt_w"] = xt_w
    airplane["yt_w"] = yt_w
    airplane["zt_w"] = zt_w
    airplane["ct_w"] = ct_w
    airplane["xm_w"] = xm_w
    airplane["ym_w"] = ym_w
    airplane["zm_w"] = zm_w
    airplane["cm_w"] = cm_w
    airplane["S_h"] = S_h
    airplane["b_h"] = b_h
    airplane["xr_h"] = xr_h
    airplane["cr_h"] = cr_h
    airplane["xt_h"] = xt_h
    airplane["yt_h"] = yt_h
    airplane["zt_h"] = zt_h
    airplane["ct_h"] = ct_h
    airplane["xm_h"] = xm_h
    airplane["ym_h"] = ym_h
    airplane["zm_h"] = zm_h
    airplane["cm_h"] = cm_h
    airplane["S_v"] = S_v
    airplane["b_v"] = b_v
    airplane["xr_v"] = xr_v
    airplane["cr_v"] = cr_v
    airplane["xt_v"] = xt_v
    airplane["zt_v"] = zt_v
    airplane["ct_v"] = ct_v
    airplane["xm_v"] = xm_v
    airplane["zm_v"] = zm_v
    airplane["cm_v"] = cm_v

    # All variables are stored in the dictionary.
    # There is no need to return anything
    return None


# ----------------------------------------


def aerodynamics(
    airplane,
    Mach,
    altitude,
    CL,
    W0_guess,
    n_engines_failed=0,
    highlift_config="clean",
    lg_down=0,
    h_ground=0,
    method=2,
):
    """
    Mach: float -> Freestream Mach number.

    altitude: float -> Flight altitude [meters].

    CL: float -> Lift coefficient

    W0_guess: float -> Latest MTOW estimate [N]

    n_engines_failed: integer -> number of engines failed. Windmilling drag is
                                 added here. This number should be less than the
                                 total number of engines.

    highlift_config: 'clean', 'takeoff', or 'landing' -> Configuration of high-lift devices

    lg_down: 0 or 1 -> 0 for retraced landing gear or 1 for extended landing gear

    h_ground: float -> Distance between wing and the ground for ground effect [m].
                       Use 0 for no ground effect.

    method: 1 or 2 -> Method 1 applies a single friction coefficient
                      to the entire wetted area of the aircraft (based on Howe).
                      Method 2 is more refined since it computes friction and
                      form factors for each component.
    """

    # Wetted areas from Torenbeek's Appendix B

    # Unpacking dictionary
    S_w = airplane["S_w"]
    AR_w = airplane["AR_w"]
    cr_w = airplane["cr_w"]
    ct_w = airplane["ct_w"]
    taper_w = airplane["taper_w"]
    sweep_w = airplane["sweep_w"]
    tcr_w = airplane["tcr_w"]
    tct_w = airplane["tct_w"]
    b_w = airplane["b_w"]
    cm_w = airplane["cm_w"]

    clmax_w = airplane["clmax_w"]

    S_h = airplane["S_h"]
    cr_h = airplane["cr_h"]
    ct_h = airplane["ct_h"]
    taper_h = airplane["taper_h"]
    sweep_h = airplane["sweep_h"]
    tcr_h = airplane["tcr_h"]
    tct_h = airplane["tct_h"]
    b_h = airplane["b_h"]
    cm_h = airplane["cm_h"]

    S_v = airplane["S_v"]
    cr_v = airplane["cr_v"]
    ct_v = airplane["ct_v"]
    taper_v = airplane["taper_v"]
    sweep_v = airplane["sweep_v"]
    tcr_v = airplane["tcr_v"]
    tct_v = airplane["tct_v"]
    b_v = airplane["b_v"]
    cm_v = airplane["cm_v"]

    L_f = airplane["L_f"]
    D_f = airplane["D_f"]

    L_n = airplane["L_n"]
    D_n = airplane["D_n"]

    x_nlg = airplane["x_nlg"]  # This is only used to check if we have LG

    n_engines = airplane["n_engines"]
    n_engines_under_wing = airplane["n_engines_under_wing"]

    flap_type = airplane["flap_type"]
    c_flap_c_wing = airplane["c_flap_c_wing"]
    b_flap_b_wing = airplane["b_flap_b_wing"]

    slat_type = airplane["slat_type"]
    c_slat_c_wing = airplane["c_slat_c_wing"]
    b_slat_b_wing = airplane["b_slat_b_wing"]
    k_exc_drag = airplane["k_exc_drag"]

    # Default rugosity value (smooth paint from Raymer Tab 12.5)
    rugosity = 0.634e-5

    ### WING

    # Average t/c
    tc_avg = 0.5 * (tcr_w + tct_w)

    # Compute the wing planform area hidden by the fuselage
    D_f_b_wing = D_f / b_w
    S_hid_S_wing = D_f_b_wing * (2 - D_f_b_wing * (1 - taper_w)) / (1 + taper_w)

    # Exposed Area
    Sexp = S_w * (1 - S_hid_S_wing)

    # Wetted Area
    tau = tcr_w / tct_w
    Swet_w = 2 * Sexp * (1 + 0.25 * tcr_w * (1 + tau * taper_w) / (1 + taper_w))

    # Friction coefficient
    Cf_w = Cf_calc(Mach, altitude, length=cm_w, rugosity=rugosity, k_lam=0.05)

    # Form factor
    FF_w = FF_surface(Mach, tcr_w, tct_w, sweep_w, b_w, cr_w, ct_w, cm_w)

    # Interference factor
    Q_w = 1.0

    # Drag coefficient
    CD0_w = Cf_w * FF_w * Q_w * Swet_w / S_w

    ### HORIZONTAL TAIL

    # Exposed Area
    Sexp = S_h

    # Wetted Area
    tau = tcr_h / tct_h
    Swet_h = 2 * Sexp * (1 + 0.25 * tcr_h * (1 + tau * taper_h) / (1 + taper_h))

    # Friction coefficient
    Cf_h = Cf_calc(Mach, altitude, length=cm_h, rugosity=rugosity, k_lam=0.05)

    # Form factor
    FF_h = FF_surface(Mach, tcr_h, tct_h, sweep_h, b_h, cr_h, ct_h, cm_h)

    # Interference factor
    Q_h = 1.05

    # Drag coefficient
    CD0_h = Cf_h * FF_h * Q_h * Swet_h / S_w

    ### VERTICAL TAIL

    # Exposed Area
    Sexp = S_v

    # Wetted Area
    tau = tcr_v / tct_v
    Swet_v = 2 * Sexp * (1 + 0.25 * tcr_v * (1 + tau * taper_v) / (1 + taper_v))

    # Friction coefficient
    Cf_v = Cf_calc(Mach, altitude, length=cm_v, rugosity=rugosity, k_lam=0.05)

    # Form factor
    FF_v = FF_surface(Mach, tcr_v, tct_v, sweep_v, 2 * b_v, cr_v, ct_v, cm_v)

    # Interference factor
    Q_v = 1.05

    # Drag coefficient
    CD0_v = Cf_v * FF_v * Q_v * Swet_v / S_w

    ### FUSELAGE

    # Wetted area
    lambda_fus = L_f / D_f
    Swet_f = (
        np.pi
        * D_f
        * L_f
        * (1 - 2 / lambda_fus) ** (2.0 / 3.0)
        * (1 + 1 / lambda_fus**2)
    )

    # Friction coefficient
    Cf_f = Cf_calc(Mach, altitude, length=L_f, rugosity=rugosity, k_lam=0.05)

    # Form factor
    FF_f = 1 + 60 / lambda_fus**3 + lambda_fus / 400

    # Interference factor
    Q_f = 1.0

    # Drag coefficient
    CD0_f = Cf_f * FF_f * Q_f * Swet_f / S_w

    ### NACELLE

    # Wetted area (where we take the number of nacelles into account)
    Swet_n = n_engines * np.pi * D_n * L_n

    # Friction coefficient
    Cf_n = Cf_calc(Mach, altitude, length=L_n, rugosity=rugosity, k_lam=0.05)

    # Form factor
    lambda_n = L_n / D_n
    FF_n = 1 + 0.35 / lambda_n

    # Interference factor
    Q_n = 1.2

    # Drag coefficient
    CD0_n = Cf_n * FF_n * Q_n * Swet_n / S_w

    ### VISCOUS DRAG

    # Total wetted area
    Swet = Swet_w + Swet_h + Swet_v + Swet_f + Swet_n

    if method == 1:

        # Wetted area ratio
        Sr = Swet / S_w

        # t/c correction
        tau = (Sr - 2) / Sr + 1.9 / Sr * (1 + 0.526 * (4 * tc_avg) ** 3)

        # Other parameters for jet aircraft
        Af = 0.93
        clam = 0.05
        Tf = 1.1

        # Friction coefficient (Howe Eq 6.13)
        Cfe = (
            0.005
            * (1 - 2 * clam / Sr)
            * tau
            * (
                1
                - 0.2 * Mach
                + 0.12 * (Mach * np.sqrt(np.cos(sweep_w)) / (Af - tc_avg)) ** 20
            )
            * Tf
            * S_w ** (-0.1)
        )

        # Viscous drag
        CD0 = Cfe * Swet / S_w

    elif method == 2:

        # Add all drag coefficients
        CD0 = CD0_w + CD0_h + CD0_v + CD0_f + CD0_n

    ### INDUCED

    # Oswald Factor (Howe Eq 6.14)
    f_taper = 0.005 * (1 + 1.5 * (taper_w - 0.6) ** 2)
    e = (
        1
        / (1 + 0.12 * Mach**6)
        / (
            1
            + (0.142 + AR_w * (10 * tc_avg) ** 0.33 * f_taper) / np.cos(sweep_w) ** 2
            + 0.1 * (3 * n_engines_under_wing + 1) / (4 + AR_w) ** 0.8
        )
    )

    # Induced drag term
    K = 1 / np.pi / AR_w / e

    ### GROUND EFFECT
    if h_ground > 0:
        aux = 33 * (h_ground / b_w) ** 1.5
        Kge = aux / (1 + aux)  # Raymer Eq. 12.61
        K = K * Kge

    # Induced drag
    CDind_clean = K * CL**2

    ### WAVE DRAG (Korn Equation)

    if Mach > 0.4:

        # I adjusted the supercritical airfoil factor from 0.95 to 0.91
        Mach_dd = (
            0.91 / np.cos(sweep_w)
            - tc_avg / np.cos(sweep_w) ** 2
            - CL / 10 / np.cos(sweep_w) ** 3
        )
        Mach_crit = Mach_dd - (0.1 / 80) ** (1 / 3)

        if Mach > Mach_crit:
            CDwave = 20 * (Mach - Mach_crit) ** 4
        else:
            CDwave = 0.0

    else:
        CDwave = 0.0

    ### HIGH LIFT

    ### Clean wing CLmax (Raymer Eq. 5.7)
    CLmax_clean = 0.9 * clmax_w * np.cos(sweep_w)

    ### Flaps deflection
    if flap_type is not None:

        # Compute flapped area (Raymer Fig. 12.21)
        # Here we consider a trapezoidal wing.
        # The second term is the subtraction of the wing area inside the fuselage
        S_flap_S_wing = (b_flap_b_wing * (2 - b_flap_b_wing * (1 - taper_w))) / (
            1 + taper_w
        ) - S_hid_S_wing

        # Sweep at flap hinge line
        sweep_flap = geo_change_sweep(
            0.25, 1 - c_flap_c_wing, sweep_w, b_w / 2, cr_w, ct_w
        )

        # Take coefficients depending on the flap type
        # dclamx - Raymer
        # Fflap - Raymer Eq. 12.61
        # flap_def - Max deflection for Torenbeek flap chart
        if flap_type == "plain":
            dclmax = 0.9
            Fflap = 0.0144
            flap_def = {
                "clean": 0 * np.pi / 180,
                "takeoff": 20 * np.pi / 180,
                "landing": 60 * np.pi / 180,
            }
        elif flap_type == "slotted":
            dclmax = 1.3
            Fflap = 0.0074
            flap_def = {
                "clean": 0 * np.pi / 180,
                "takeoff": 20 * np.pi / 180,
                "landing": 40 * np.pi / 180,
            }
        elif flap_type == "fowler":
            dclmax = 1.3 * (
                1 + c_flap_c_wing
            )  # Assumed that flap extends until trailing edge
            Fflap = 0.0074
            flap_def = {
                "clean": 0 * np.pi / 180,
                "takeoff": 15 * np.pi / 180,
                "landing": 40 * np.pi / 180,
            }
        elif flap_type == "double slotted":
            dclmax = 1.6 * (
                1 + c_flap_c_wing
            )  # Assumed that flap extends until trailing edge
            Fflap = 0.0074
            flap_def = {
                "clean": 0 * np.pi / 180,
                "takeoff": 20 * np.pi / 180,
                "landing": 50 * np.pi / 180,
            }
        elif flap_type == "triple slotted":
            dclmax = 1.9 * (
                1 + c_flap_c_wing
            )  # Assumed that flap extends until trailing edge
            Fflap = 0.0074
            flap_def = {
                "clean": 0 * np.pi / 180,
                "takeoff": 20 * np.pi / 180,
                "landing": 40 * np.pi / 180,
            }

        # Factor to adjust lift contribution for intermediate deflections.
        # Based on Torenbeek's table.
        if highlift_config == "clean":
            lift_factor = 0.0
        elif highlift_config == "takeoff":
            lift_factor = 0.60
        elif highlift_config == "landing":
            lift_factor = 1.0

        deltaCLmax_flap = (
            0.9 * dclmax * S_flap_S_wing * np.cos(sweep_flap) * lift_factor
        )  # Raymer Eq 12.21

        # Parasite drag contribution (Raymer Eq 12.61)
        # Added a max to avoid negative values if flap is not deflected
        CD0_flap = max(
            0.0,
            Fflap
            * c_flap_c_wing
            * S_flap_S_wing
            * (flap_def[highlift_config] * 180 / np.pi - 10),
        )

    else:
        CD0_flap = 0.0
        deltaCLmax_flap = 0.0

    ### Slats deflection
    if slat_type is not None:

        CD0_slat = 0.0

        # Compute flapped area (Raymer Fig. 12.21)
        # Here we consider a trapezoidal wing.
        # The second term is the subtraction of the wing area inside the fuselage
        D_f_b_wing = D_f / b_w
        S_slat_S_wing = (
            b_slat_b_wing * (2 - b_slat_b_wing * (1 - taper_w))
            - D_f_b_wing * (2 - D_f_b_wing * (1 - taper_w))
        ) / (1 + taper_w)

        sweep_slat = geo_change_sweep(0.25, c_slat_c_wing, sweep_w, b_w / 2, cr_w, ct_w)

        if slat_type == "fixed":
            dclmax = 0.2
        elif slat_type == "flap":
            dclmax = 0.3
        elif slat_type == "kruger":
            dclmax = 0.3
        elif slat_type == "slat":
            dclmax = 0.4 * (1 + c_slat_c_wing)

        # Factor to adjust lift contribution for intermediate deflections.
        # Based on Torenbeek's table.
        if highlift_config == "clean":
            lift_factor = 0.0
        elif highlift_config == "takeoff":
            lift_factor = 0.60
        elif highlift_config == "landing":
            lift_factor = 1.0

        deltaCLmax_slat = (
            0.9 * dclmax * S_slat_S_wing * np.cos(sweep_slat) * lift_factor
        )  # Raymer Eq 12.21

    else:
        CD0_slat = 0.0
        deltaCLmax_slat = 0.0

    # Maximum lift
    CLmax = CLmax_clean + deltaCLmax_flap + deltaCLmax_slat

    # Induced drag contribution due to high-lift devices (Raymer Eq. 12.62)
    # This term seemed too prohibitive, penalizing climb gradients. So I took the
    # average of the values proposed by Raymer.
    CDind_flap = (0.22 * (deltaCLmax_flap + deltaCLmax_slat)) ** 2 * np.cos(sweep_w)

    # Update total induced drag
    CDind = CDind_clean + CDind_flap

    ### Landing gear
    if x_nlg is not None:  # Check if we have a LG

        # (ESDU)
        if flap_type is not None:
            lg_factor = (
                0.57 - 0.26 * flap_def[highlift_config] / flap_def["landing"]
            ) * 1e-3
        else:
            lg_factor = 0.57e-3
        CD0_lg = lg_down * lg_factor * (W0_guess / gravity) ** 0.785 / S_w

    else:
        CD0_lg = 0.0

    ### Windmill engine
    # Vn_V = 0.42
    # CDwdm = (0.0785*D_n**2 + 1/(1 + 0.16*Mach**2)*np.pi/2*D_n**2*Vn_V*(1-Vn_V))/S_w
    # CD0_wdm = n_engines_failed*CDwdm
    CD0_wdm = n_engines_failed * 0.3 * np.pi / 4 * D_n**2 / S_w  # Raymer Eq 12.40

    # Add all parasite drag values found so far
    CD0 = CD0 + CD0_flap + CD0_slat + CD0_lg + CD0_wdm

    ### Excrescence
    CD0_exc = CD0 * k_exc_drag / (1 - k_exc_drag)
    CD0 = CD0 + CD0_exc

    # Total drag
    CD = CD0 + CDind + CDwave

    # Create a drag breakdown dictionary
    dragDict = {
        "CD0_lg": CD0_lg,
        "CD0_wdm": CD0_wdm,
        "CD0_exc": CD0_exc,
        "CD0_flap": CD0_flap,
        "CD0_slat": CD0_slat,
        "CD0": CD0,
        "CDind_clean": CDind_clean,
        "CDind_flap": CDind_flap,
        "CDind": CDind,
        "CDwave": CDwave,
        "CLmax_clean": CLmax_clean,
        "deltaCLmax_flap": deltaCLmax_flap,
        "deltaCLmax_slat": deltaCLmax_slat,
        "CLmax": CLmax,
        "K": K,
        "Swet": Swet,
    }

    if method == 2:
        dragDict["CD0_w"] = CD0_w
        dragDict["CD0_h"] = CD0_h
        dragDict["CD0_v"] = CD0_v
        dragDict["CD0_f"] = CD0_f
        dragDict["CD0_n"] = CD0_n

    # Update dictionary
    airplane["Swet_f"] = Swet_f

    return CD, CLmax, dragDict


# ----------------------------------------


def engineTSFC(Mach, altitude, airplane):
    """
    This function computes the engine thrust-specific fuel
    consumption and thrust correction factor compared to
    static sea-level conditions. The user has to define the
    engine parameters in a 'engine' dictionary within
    the airplane dictionary. The engine model must be
    identified by the 'model' field of the engine dictionary.
    The following engine models are available:

    Howe TSFC turbofan model:
    requires the bypass ratio. An optional sea-level TSFC
    could also be provided. Otherwise, standard parameters
    are used.
    airplane['engine'] = {'model': 'howe turbofan',
                          'BPR': 3.04,
                          'Cbase': 0.7/3600} # Could also be None

    Thermodynamic cycle turbojet:
    This model uses a simplified thermodynamic model of
    turbofans to estimate maximum thrust and TSFC

    airplane['engine'] = {'model': 'thermo turbojet'
                          'data': dictionary (check turbojet_model function)}

    The user can also leave a 'weight' field in the dictionary
    to replace the weight estimation.
    """

    # Get a reference to the engine dictionary
    engine = airplane["engine"]

    # Check which model was given
    if engine["model"].lower() == "howe turbofan":

        # Unpack dictionary
        BPR = engine["BPR"]
        Cbase = engine["Cbase"]  # This is sea-level static TSFC

        # Atmospheric conditions at cruise altitude
        T, p, rho, mi = atmosphere(altitude)

        # Density ratio
        sigma = rho / 1.225

        # Base TSFC
        if Cbase is None:  # User did not provide a base TSFC
            if BPR < 4.0:
                Cbase = 0.85 / 3600
            else:
                Cbase = 0.70 / 3600

        else:

            # Correct Cbase so that the equation gives the desired static TSFC at sea-level
            Cbase = Cbase / (1 - 0.15 * BPR**0.65)

        # Howe Eq 3.12a
        C = (
            Cbase
            * (1 - 0.15 * BPR**0.65)
            * (1 + 0.28 * (1 + 0.063 * BPR**2) * Mach)
            * sigma**0.08
        )

        # Cruise traction correction for takeoff conditions by Scholz
        # (https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_5_PreliminarySizing.pdf)
        kT = (0.0013 * BPR - 0.0397) * altitude / 1000.0 - 0.0248 * BPR + 0.7125

    elif engine["model"].lower() == "thermo turbojet":

        C, F = turbojet_model(Mach, altitude, engine["data"])

        # Check if maximum sea-level thrust was already computed
        if "T0_eng" not in engine:

            _, F0 = turbojet_model(0.0, 0.0, engine["data"])
            airplane["engine"]["T0_eng"] = F0

        # Compute thrust correction factor
        kT = F / engine["T0_eng"]

    return C, kT


# ----------------------------------------


def empty_weight(W0_guess, T0_guess, airplane):

    # Unpack dictionary
    S_w = airplane["S_w"]
    AR_w = airplane["AR_w"]
    taper_w = airplane["taper_w"]
    sweep_w = airplane["sweep_w"]
    xm_w = airplane["xm_w"]
    cm_w = airplane["cm_w"]
    tcr_w = airplane["tcr_w"]
    S_h = airplane["S_h"]
    xm_h = airplane["xm_h"]
    cm_h = airplane["cm_h"]
    S_v = airplane["S_v"]
    xm_v = airplane["xm_v"]
    cm_v = airplane["cm_v"]
    L_f = airplane["L_f"]
    Swet_f = airplane["Swet_f"]
    n_engines = airplane["n_engines"]
    x_n = airplane["x_n"]
    L_n = airplane["L_n"]
    x_nlg = airplane["x_nlg"]
    x_mlg = airplane["x_mlg"]
    altitude_cruise = airplane["altitude_cruise"]
    Mach_cruise = airplane["Mach_cruise"]

    # Select appropriate parameters for weight regression

    if airplane["type"] == "transport":

        # Wing weight (Raymer Eq 15.25)
        # I increased the AR_w exponent from 0.5 to 0.55 to make it more sensitive.
        # Otherwise, the optimum would be around AR = 12, which may be too optimistic.
        Nz = 1.5 * 2.5  # Ultimate load factor
        Scsw = 0.15 * S_w  # Area of control surfaces
        W_w = (
            0.0051
            * (W0_guess * Nz / lb2N) ** 0.557
            * (S_w / ft2m**2) ** 0.649
            * AR_w**0.55
            * tcr_w ** (-0.4)
            * (1 + taper_w) ** 0.1
            / np.cos(sweep_w)
            * (Scsw / ft2m**2) ** 0.1
            * lb2N
        )
        xcg_w = xm_w + 0.4 * cm_w

        # Surface densities for remaining components (kg/m2) - Raymer Tab 15.2
        W_h_dens = 27
        W_v_dens = 27
        W_f_dens = 24
        W_lg_fact = 0.043
        W_eng_fact = 1.3
        W_allelse_fact = 0.17

    elif airplane["type"] == "fighter":

        # Wing weight (Raymer Eq 15.1)
        Nz = 1.5 * 6.0  # Ultimate load factor
        Scsw = 0.15 * S_w  # Area of control surfaces
        W_w = (
            0.0103
            * (W0_guess * Nz / lb2N) ** 0.5
            * (S_w / ft2m**2) ** 0.622
            * AR_w**0.785
            * tcr_w ** (-0.4)
            * (1 + taper_w) ** 0.05
            / np.cos(sweep_w)
            * (Scsw / ft2m**2) ** 0.04
            * lb2N
        )
        xcg_w = xm_w + 0.4 * cm_w

        # Surface densities for remaining components (kg/m2) - Raymer Tab 15.2
        W_h_dens = 20
        W_v_dens = 26
        W_f_dens = 23
        W_lg_fact = 0.033
        W_eng_fact = 1.3
        W_allelse_fact = 0.17

    elif airplane["type"] == "general":

        # Cruise dynamic pressure
        T, _, rho, _ = atmosphere(altitude_cruise)
        a_cruise = np.sqrt(gamma_air * R_air * T)
        v_cruise = Mach_cruise * a_cruise
        q_cruise = 0.5 * rho * v_cruise**2

        # Wing weight (Raymer Eq 15.46)
        Nz = 1.5 * 4.4  # Ultimate load factor
        W_w = (
            0.036
            * (W0_guess * Nz / lb2N) ** 0.49
            * (S_w / ft2m**2) ** 0.758
            * (AR_w / np.cos(sweep_w) ** 2) ** 0.6
            * (100 * tcr_w / np.cos(sweep_w)) ** (-0.3)
            * (taper_w) ** 0.04
            * q_cruise**0.006
            * lb2N
        )
        xcg_w = xm_w + 0.4 * cm_w

        # Surface densities for remaining components (kg/m2) - Raymer Tab 15.2
        W_h_dens = 10
        W_v_dens = 10
        W_f_dens = 7
        W_lg_fact = 0.057
        W_eng_fact = 1.4
        W_allelse_fact = 0.1

    # Use Raymer Tab 15.2 for the remaining components

    W_h = S_h * gravity * W_h_dens
    xcg_h = xm_h + 0.4 * cm_h

    W_v = S_v * gravity * W_v_dens
    xcg_v = xm_v + 0.4 * cm_v

    W_f = Swet_f * gravity * W_f_dens
    xcg_f = 0.45 * L_f

    # Check if LG is active
    if x_nlg is not None:

        W_nlg = 0.15 * W0_guess * W_lg_fact
        xcg_nlg = x_nlg

        W_mlg = 0.85 * W0_guess * W_lg_fact
        xcg_mlg = x_mlg

    else:

        W_nlg = 0.0
        xcg_nlg = 0.0

        W_mlg = 0.0
        xcg_mlg = 0.0

    # Engine weight
    T_eng = T0_guess / n_engines

    if "weight" in airplane["engine"]:

        # The user already gave the engine weight
        W_eng = airplane["engine"]["weight"]

    elif "turbofan" in airplane["engine"]["model"]:
        BPR = airplane["engine"]["BPR"]

        # Turbofan weight (Raymer Eq. 10.4)
        W_eng = gravity * 14.7 * (T_eng / 1000.0) ** 1.1 * np.exp(-0.045 * BPR)

    W_eng_installed = n_engines * W_eng * W_eng_fact
    xcg_eng = x_n + 0.5 * L_n

    # All else weight
    W_allelse = W_allelse_fact * W0_guess
    xcg_allelse = 0.45 * L_f

    # Empty weight
    We = W_w + W_h + W_v + W_f + W_nlg + W_mlg + W_eng_installed + W_allelse

    # Empty weight CG
    xcg_e = (
        W_w * xcg_w
        + W_h * xcg_h
        + W_v * xcg_v
        + W_f * xcg_f
        + W_nlg * xcg_nlg
        + W_mlg * xcg_mlg
        + W_eng_installed * xcg_eng
        + W_allelse * xcg_allelse
    ) / We

    # Update dictionary
    airplane["W_w"] = W_w
    airplane["W_h"] = W_h
    airplane["W_v"] = W_v
    airplane["W_f"] = W_f
    airplane["W_nlg"] = W_nlg
    airplane["W_mlg"] = W_mlg
    airplane["W_eng"] = W_eng_installed
    airplane["W_allelse"] = W_allelse

    return We, xcg_e


# ----------------------------------------


def fuel_weight(W0_guess, airplane):

    # Unpacking dictionary
    S_w = airplane["S_w"]
    altitude_cruise = airplane["altitude_cruise"]
    Mach_cruise = airplane["Mach_cruise"]
    range_cruise = airplane["range_cruise"]
    loiter_time = airplane["loiter_time"]
    altitude_altcruise = airplane["altitude_altcruise"]
    Mach_altcruise = airplane["Mach_altcruise"]
    range_altcruise = airplane["range_altcruise"]

    # Get engine TSFC
    C_cruise, _ = engineTSFC(Mach_cruise, altitude_cruise, airplane)
    C_altcruise, _ = engineTSFC(Mach_altcruise, altitude_altcruise, airplane)

    # Initialize product of all phases
    Mf = 1.0

    ### Start and warm-up
    Mf = Mf * 0.99

    ### Taxi
    Mf = Mf * 0.99

    ### Take-off
    Mf = Mf * 0.995

    ### Climb
    Mf = Mf * 0.98

    ### Cruise

    # Store weight fraction at beginning of the cruise
    Mf_cruise = Mf

    # Atmospheric conditions at cruise altitude
    T, p, rho, mi = atmosphere(altitude_cruise)

    # Cruise speed
    a_cruise = np.sqrt(gamma_air * R_air * T)
    v_cruise = Mach_cruise * a_cruise

    # Cruise CL
    CL = 2.0 * W0_guess * Mf / rho / S_w / v_cruise**2

    # Cruise C
    CD, CLmax, dragDict = aerodynamics(
        airplane,
        Mach_cruise,
        altitude_cruise,
        CL,
        W0_guess,
        n_engines_failed=0,
        highlift_config="clean",
        lg_down=0,
        h_ground=0,
    )

    Mf = Mf * np.exp(-range_cruise * C_cruise / v_cruise * CD / CL)

    ### Loiter

    # Loiter at max L/D #ESCREVER
    # For now, we take the cruise CD0 and K (which do not take
    # wave drag into account) to estimate L/Dmax
    LDmax = 0.5 / np.sqrt(dragDict["CD0"] * dragDict["K"])

    # Factor to fuel comsumption (Correction based on Raymer Tab. 3.3)
    C_loiter = C_cruise - 0.1 / 3600.0

    Mf = Mf * np.exp(-loiter_time * C_loiter / LDmax)

    ### Descent
    Mf = Mf * 0.99

    ### Cruise 2

    # Atmospheric conditions at cruise altitude
    T, p, rho, mi = atmosphere(altitude_altcruise)

    # Cruise speed
    a_altcruise = np.sqrt(gamma_air * R_air * T)
    v_altcruise = Mach_altcruise * a_altcruise

    # Cruise CL
    CL = 2.0 * W0_guess * Mf / rho / S_w / v_altcruise**2

    # Cruise CD
    CD, CLmax, dragDict = aerodynamics(
        airplane,
        Mach_altcruise,
        altitude_altcruise,
        CL,
        W0_guess,
        n_engines_failed=0,
        highlift_config="clean",
        lg_down=0,
        h_ground=0,
    )

    Mf = Mf * np.exp(-range_altcruise * C_altcruise / v_altcruise * CD / CL)

    ### Landing and Taxi
    Mf = Mf * 0.992

    ### Fuel weight (Raymer Eq 3.13)
    Wf = 1.06 * (1 - Mf) * W0_guess

    return Wf, Mf_cruise


# ----------------------------------------


def weight(W0_guess, T0_guess, airplane):

    # Unpacking dictionary
    W_payload = airplane["W_payload"]
    W_crew = airplane["W_crew"]

    # Set iterator
    delta = 1000

    while abs(delta) > 10:

        # We need to call fuel_weight first since it
        # calls the aerodynamics module to get Swet_f used by
        # the empty weight function
        Wf, Mf_cruise = fuel_weight(W0_guess, airplane)

        We, xcg_e = empty_weight(W0_guess, T0_guess, airplane)

        W0 = We + Wf + W_payload + W_crew

        delta = W0 - W0_guess

        W0_guess = W0

    return W0, We, Wf, Mf_cruise, xcg_e


# ----------------------------------------


def performance(W0, Mf_cruise, airplane):

    """
    This function computes the required thrust and wing areas
    required to meet takeoff, landing, climb, and cruise requirements.

    OUTPUTS:
    T0: real -> Total thrust required to meet all mission phases
    deltaS_wlan: real -> Wing area margin for landing. This value should be positive
                         for a feasible landing.
    """

    # Unpacking dictionary
    S_w = airplane["S_w"]
    n_engines = airplane["n_engines"]
    h_ground = airplane["h_ground"]
    altitude_takeoff = airplane["altitude_takeoff"]
    distance_takeoff = airplane["distance_takeoff"]
    altitude_landing = airplane["altitude_landing"]
    distance_landing = airplane["distance_landing"]
    MLW_frac = airplane["MLW_frac"]
    altitude_cruise = airplane["altitude_cruise"]
    Mach_cruise = airplane["Mach_cruise"]

    ### TAKEOFF

    # Compute air density at takeoff altitude
    T, p, rho, mi = atmosphere(altitude_takeoff)

    # density ratio
    sigma = rho / 1.225

    # Takeoff aerodynamics
    _, CLmaxTO, _ = aerodynamics(
        airplane,
        Mach=0.2,
        altitude=altitude_takeoff,
        CL=0.5,
        W0_guess=W0,
        n_engines_failed=0,
        highlift_config="takeoff",
        lg_down=0,
        h_ground=h_ground,
    )

    T0W0 = 0.2387 / sigma / CLmaxTO / distance_takeoff * W0 / S_w

    T0_to = T0W0 * W0

    ### LANDING

    # Compute air density at landing altitude
    T, p, rho, mi = atmosphere(altitude_landing)

    # Landing aerodynamics
    _, CLmaxLD, _ = aerodynamics(
        airplane,
        Mach=0.2,
        altitude=altitude_landing,
        CL=0.5,
        W0_guess=W0,
        n_engines_failed=0,
        highlift_config="landing",
        lg_down=0,
        h_ground=h_ground,
    )

    # Landing Field Length (Roskam)
    # We removed this since it was favoring the landing (SFL)
    # distance_landing = distance_landing/0.6

    # Approach speed (Roskam adapted to SI)
    Va = 1.701 * np.sqrt(distance_landing)

    # Required stall speed
    Vs = Va / 1.3

    # Required wing area
    S_wlan = 2 * W0 * MLW_frac / rho / Vs**2 / CLmaxLD

    # Compute wing area excess with respect to the landing requirement.
    # The aircraft should have deltaS_wlan >= 0 to satisfy landing.
    deltaS_wlan = S_w - S_wlan

    ### CRUISE

    # Compute air density at cruise altitude
    T, p, rho, mi = atmosphere(altitude_cruise)

    # Cruise speed
    a_cruise = np.sqrt(gamma_air * R_air * T)
    v_cruise = Mach_cruise * a_cruise

    # Cruise CL
    CL = 2.0 * W0 * Mf_cruise / rho / S_w / v_cruise**2

    # Cruise CD # ESCREVER
    CD, _, _ = aerodynamics(
        airplane,
        Mach=Mach_cruise,
        altitude=altitude_cruise,
        CL=CL,
        W0_guess=W0,
        n_engines_failed=0,
        highlift_config="clean",
        lg_down=0,
        h_ground=0,
    )

    # Cruise traction
    T = 0.5 * rho * v_cruise**2 * S_w * CD

    # Cruise correction factor
    _, kT = engineTSFC(Mach_cruise, altitude_cruise, airplane)

    # Corrected thrust
    T0_cruise = T / kT

    ### CLIMB

    # Define standard function for climb analysis
    def climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    ):

        """
        We need a guess for CLmax just to get an approximate drag polar for
        speed computation. We will get the correct CLmax from the aerodynamics module

        kT: Thrust decay factor (e.g. use 0.94 for maximum continuous thrust)
        """

        # Compute air density
        T, p, rho, mi = atmosphere(altitude)

        # Compute stall speed
        Vs = np.sqrt(2 * W0 * Mf / rho / S_w / CLmax_guess)

        # Compute climb speed
        Vclimb = Vs * Ks

        # Compute sound speed and Mach number
        a = np.sqrt(gamma_air * R_air * T)
        Mach = Vclimb / a

        # Dummy run to get CLmax # ESCREVER
        _, CLmax, _ = aerodynamics(
            airplane,
            Mach=Mach,
            altitude=altitude,
            CL=0.5,
            W0_guess=W0,
            n_engines_failed=n_engines_failed,
            highlift_config=highlift_config,
            lg_down=lg_down,
            h_ground=h_ground_climb,
        )

        # Get climb CL
        CL = CLmax / Ks**2

        # Get corresponding CD # ESCREVER
        CD, _, _ = aerodynamics(
            airplane,
            Mach=Mach,
            altitude=altitude,
            CL=CL,
            W0_guess=W0,
            n_engines_failed=n_engines_failed,
            highlift_config=highlift_config,
            lg_down=lg_down,
            h_ground=h_ground_climb,
        )

        # Check number of failed engines
        if n_engines_failed >= n_engines:
            print(
                "Warning: number of engines failed is equal or greater than the number of engines"
            )
            print("We will force n_engines_failed = n_engines-1")
            n_engines_failed = n_engines - 1

        # Compute T/W
        TW = n_engines / (n_engines - n_engines_failed) * (grad + CD / CL)

        # Compute required traction
        T0 = TW * W0 * Mf / kT

        return T0

    # FAR 25.111
    grad = 0.012
    Ks = 1.2
    altitude = altitude_takeoff
    CLmax_guess = CLmaxTO
    lg_down = 0
    h_ground_climb = h_ground
    highlift_config = "takeoff"
    n_engines_failed = 1
    Mf = 1.0
    kT = 1.0
    T0_1 = climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    )

    # FAR 25.121a
    grad = 0.0
    Ks = 1.1
    altitude = altitude_takeoff
    CLmax_guess = CLmaxTO
    lg_down = 1
    h_ground_climb = h_ground
    highlift_config = "takeoff"
    n_engines_failed = 1
    Mf = 1.0
    kT = 1.0
    T0_2 = climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    )

    # FAR 25.121b
    grad = 0.024
    Ks = 1.2
    altitude = altitude_takeoff
    CLmax_guess = CLmaxTO
    lg_down = 0
    h_ground_climb = 0
    highlift_config = "takeoff"
    n_engines_failed = 1
    Mf = 1.0
    kT = 1.0
    T0_3 = climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    )

    # FAR 25.121c
    grad = 0.012
    Ks = 1.25
    altitude = altitude_takeoff
    CLmax_guess = CLmaxTO
    lg_down = 0
    h_ground_climb = 0
    highlift_config = "clean"
    n_engines_failed = 1
    Mf = 1.0
    kT = 0.94
    T0_4 = climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    )

    # FAR 25.119
    grad = 0.032
    Ks = 1.30
    altitude = altitude_landing
    CLmax_guess = CLmaxLD
    lg_down = 1
    h_ground_climb = 0
    highlift_config = "landing"
    n_engines_failed = 0
    Mf = MLW_frac
    kT = 1.0
    T0_5 = climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    )

    # FAR 25.121d
    grad = 0.021
    Ks = 1.40
    altitude = altitude_landing
    CLmax_guess = CLmaxLD
    lg_down = 1
    h_ground_climb = 0
    highlift_config = "takeoff"  # Assumed approach == landing
    n_engines_failed = 1
    Mf = MLW_frac
    kT = 1.0
    T0_6 = climb_analysis(
        grad,
        Ks,
        altitude,
        CLmax_guess,
        lg_down,
        h_ground_climb,
        highlift_config,
        n_engines_failed,
        Mf,
        kT,
    )

    # RATO TAKEOFF

    if "altitude_rato_takeoff" in airplane.keys():
        """
        RATO takeoff formulation from
        SILVA, Daniel Marques. Modelo de Desempenho para Avaliar o Perfil de Missão de Um
        AAM Aplicável ao Projeto Conceitual da Aeronave. 2021. Trabalho de Conclusão
        de Curso (Graduação) { Instituto Tecnológico de Aeronáutica, São José dos Campos.

        This will be assigned as a wing loading requirement that could replace
        the landing requirement if it is more critical
        """

        # Unpack optional arguments
        altitude_rato = airplane["altitude_rato_takeoff"]
        ramp_angle_rato = airplane["ramp_angle_rato_takeoff"]
        c_rato = airplane["c_rato_takeoff"]
        prop_mass_ratio_rato = airplane["prop_mass_ratio_rato_takeoff"]
        t_burn_rato = airplane["t_burn_rato_takeoff"]
        Ks = airplane["Ks_rato_takeoff"]

        # Compute air density at takeoff altitude
        T, p, rho, mi = atmosphere(altitude_rato)

        # density ratio
        sigma = rho / 1.225

        # Takeoff aerodynamics
        _, CLmaxTO, _ = aerodynamics(
            airplane,
            Mach=0.2,
            altitude=altitude_rato,
            CL=0.5,
            W0_guess=W0,
            n_engines_failed=0,
            highlift_config="takeoff",
            lg_down=0,
            h_ground=0,
        )

        # Compute required wing loading to reach desired stall speed
        # factor at the end of the burn
        WS_req = (
            0.5
            * rho
            * CLmaxTO
            / Ks
            * (
                c_rato * np.log(1 - prop_mass_ratio_rato)
                + gravity * t_burn_rato * np.sin(ramp_angle_rato)
            )
            ** 2
        )

        deltaS_rato = S_w - W0 / WS_req
        deltaS_wlan = min(deltaS_wlan, deltaS_rato)

        # Check climb gradient at the end of the burn
        # This will replace last climb requirement
        grad = np.tan(ramp_angle_rato)
        Ks = Ks
        altitude = altitude_rato
        CLmax_guess = CLmaxTO
        lg_down = 1
        h_ground_climb = 0
        highlift_config = "takeoff"
        n_engines_failed = 0
        Mf = 1.0
        kT = 1.0
        T0_6 = climb_analysis(
            grad,
            Ks,
            altitude,
            CLmax_guess,
            lg_down,
            h_ground_climb,
            highlift_config,
            n_engines_failed,
            Mf,
            kT,
        )

    # Get the maximum required thrust with a 5% margin
    T0vec = [T0_to, T0_cruise, T0_1, T0_2, T0_3, T0_4, T0_5, T0_6]
    T0 = 1.05 * max(T0vec)

    return T0, T0vec, deltaS_wlan, CLmaxTO


# ----------------------------------------


def thrust_matching(W0_guess, T0_guess, airplane):

    # Set iterator
    delta = 1000

    # Loop to adjust T0
    while abs(delta) > 10:

        W0, We, Wf, Mf_cruise, xcg_e = weight(W0_guess, T0_guess, airplane)

        T0, T0vec, deltaS_wlan, CLmaxTO = performance(W0, Mf_cruise, airplane)

        # Compute change with respect to previous iteration
        delta = T0 - T0_guess

        # Update guesses for the next iteration
        T0_guess = T0
        W0_guess = W0

    # Update dictionary
    airplane["W0"] = W0
    airplane["We"] = We
    airplane["Wf"] = Wf
    airplane["xcg_e"] = xcg_e
    airplane["T0"] = T0
    airplane["T0vec"] = T0vec
    airplane["deltaS_wlan"] = deltaS_wlan
    airplane["CLmaxTO"] = CLmaxTO

    # Return
    return None


# ----------------------------------------


def balance(airplane):

    # Unpack dictionary
    W0 = airplane["W0"]
    W_payload = airplane["W_payload"]
    xcg_payload = airplane["xcg_payload"]
    W_crew = airplane["W_crew"]
    xcg_crew = airplane["xcg_crew"]
    We = airplane["We"]
    xcg_e = airplane["xcg_e"]
    Wf = airplane["Wf"]
    Mach_cruise = airplane["Mach_cruise"]
    S_w = airplane["S_w"]
    AR_w = airplane["AR_w"]
    sweep_w = airplane["sweep_w"]
    b_w = airplane["b_w"]
    xr_w = airplane["xr_w"]
    cr_w = airplane["cr_w"]
    ct_w = airplane["ct_w"]
    xm_w = airplane["xm_w"]
    cm_w = airplane["cm_w"]
    tcr_w = airplane["tcr_w"]
    tct_w = airplane["tct_w"]
    c_tank_c_w = airplane["c_tank_c_w"]
    x_tank_c_w = airplane["x_tank_c_w"]
    S_h = airplane["S_h"]
    AR_h = airplane["AR_h"]
    sweep_h = airplane["sweep_h"]
    b_h = airplane["b_h"]
    cr_h = airplane["cr_h"]
    ct_h = airplane["ct_h"]
    xm_h = airplane["xm_h"]
    cm_h = airplane["cm_h"]
    eta_h = airplane["eta_h"]
    Cvt = airplane["Cvt"]
    L_f = airplane["L_f"]
    D_f = airplane["D_f"]
    y_n = airplane["y_n"]
    T0 = airplane["T0"]
    n_engines = airplane["n_engines"]
    CLmaxTO = airplane["CLmaxTO"]
    rho_f = airplane["rho_f"]

    ### TANK CG
    """
    We will compute the centroid of a trapezoidal tank.
    The expressions below where derived with the obelisk volume
    and centroid, assuming that both root and tip sections
    have the same c_tank_c_w and t/c
    """

    # Required fuel volume
    Vf = Wf / rho_f / gravity

    # Average wing thickness
    tc_w = 0.5 * (tcr_w + tct_w)

    # Find the span fraction that should be occupied by the fuel tank
    b_tank_b_w = (
        3.0 * Vf / c_tank_c_w / tc_w / (cr_w**2 + ct_w**2 + cr_w * ct_w) / b_w
    )

    # Find the lateral distance of the fuel tank centroid to the symmetry plane
    ycg_f = (
        b_tank_b_w
        * b_w
        / 8
        * (cr_w**2 + 2 * cr_w * ct_w + 3 * ct_w**2)
        / (cr_w**2 + cr_w * ct_w + ct_w**2)
    )

    # Sweep at the tank center line
    sweep_tank = geo_change_sweep(
        0.25, x_tank_c_w + 0.5 * c_tank_c_w, sweep_w, b_w / 2, cr_w, ct_w
    )

    # Longitudinal position of the tank CG
    xcg_f = xr_w + cr_w * (x_tank_c_w + 0.5 * c_tank_c_w) + ycg_f * np.tan(sweep_tank)

    ### CG RANGE

    # Empty airplane
    xcg_1 = xcg_e

    # Crew
    xcg_2 = (We * xcg_e + W_crew * xcg_crew) / (We + W_crew)

    # Payload and crew
    xcg_3 = (We * xcg_e + W_payload * xcg_payload + W_crew * xcg_crew) / (
        We + W_payload + W_crew
    )

    # Fuel and crew
    xcg_4 = (We * xcg_e + Wf * xcg_f + W_crew * xcg_crew) / (We + Wf + W_crew)

    # Payload, crew, and fuel (full airplane)
    xcg_5 = (We * xcg_e + Wf * xcg_f + W_payload * xcg_payload + W_crew * xcg_crew) / W0

    # Find CG range
    xcg_list = [xcg_1, xcg_2, xcg_3, xcg_4, xcg_5]
    xcg_fwd = min(xcg_list)
    xcg_aft = max(xcg_list)

    # We do not need to consider the static margin for the empty case
    # So we compute the flight CG range
    xcg_list = [xcg_2, xcg_3, xcg_4, xcg_5]
    xcg_fwd_flight = min(xcg_list)
    xcg_aft_flight = max(xcg_list)

    ### NEUTRAL POINT

    # Wing lift slope (Raymer Eq 12.6)
    sweep_maxt_w = geo_change_sweep(
        0.25, 0.40, sweep_w, b_w / 2, cr_w, ct_w
    )  # Sweep at max. thickness
    beta2 = 1 - Mach_cruise**2
    CLa_w = (
        2
        * np.pi
        * AR_w
        / (
            2
            + np.sqrt(
                4
                + AR_w**2
                * beta2
                / 0.95**2
                * (1 + np.tan(sweep_maxt_w) ** 2 / beta2)
            )
        )
        * 0.98
    )

    # Wing aerodynamic center at 25% mac
    xac_w = xm_w + 0.25 * cm_w

    # HT lift slope (Raymer Eq 12.6)
    sweep_maxt_h = geo_change_sweep(
        0.25, 0.40, sweep_h, b_h / 2, cr_h, ct_h
    )  # Sweep at max. thickness
    CLa_h = (
        2
        * np.pi
        * AR_h
        / (
            2
            + np.sqrt(
                4
                + AR_h**2
                * beta2
                / 0.95**2
                * (1 + np.tan(sweep_maxt_h) ** 2 / beta2)
            )
        )
        * 0.98
    )

    # HT aerodynamic center at 25% mac
    xac_h = xm_h + 0.25 * cm_h

    # Downwash (Nelson Eq 2.23)
    deda = 2 * CLa_w / np.pi / AR_w

    # Fuselage moment slope (Raymer Eq 16.25)
    CMa_f = 0.03 * 180 / np.pi * D_f**2 * L_f / cm_w / S_w

    # Neutral point position (Raymer Eq 16.9 and Eq 16.23)
    xnp = (
        CLa_w * xac_w - CMa_f * cm_w + eta_h * S_h / S_w * CLa_h * (1 - deda) * xac_h
    ) / (CLa_w + eta_h * S_h / S_w * CLa_h * (1 - deda))

    # Static margin
    SM_fwd = (xnp - xcg_fwd_flight) / cm_w
    SM_aft = (xnp - xcg_aft_flight) / cm_w

    # VERTICAL TAIL VERIFICATION FOR OEI CONDITION

    # Compute stall factor assuming that V2 = 1.1*Vmc and V2=1.2*Vs (FAR 25.107)
    Ks = 1.2 / 1.1

    # Compute required lift for the vertical tail
    CLv = y_n / b_w * CLmaxTO / Ks**2 * T0 / W0 / n_engines / Cvt

    # Update dictionary
    airplane["xcg_fwd"] = xcg_fwd
    airplane["xcg_aft"] = xcg_aft
    airplane["xnp"] = xnp
    airplane["SM_fwd"] = SM_fwd
    airplane["SM_aft"] = SM_aft
    airplane["b_tank_b_w"] = b_tank_b_w
    airplane["CLv"] = CLv

    return None


# ----------------------------------------


def landing_gear(airplane):

    # Unpack dictionary
    x_nlg = airplane["x_nlg"]
    x_mlg = airplane["x_mlg"]
    y_mlg = airplane["y_mlg"]
    z_lg = airplane["z_lg"]
    xcg_fwd = airplane["xcg_fwd"]
    xcg_aft = airplane["xcg_aft"]
    x_tailstrike = airplane["x_tailstrike"]
    z_tailstrike = airplane["z_tailstrike"]

    # Check if there is a landing gear
    if x_nlg is not None:

        # Weight fractions on NLG for both load cases
        frac_nlg_fwd = (x_mlg - xcg_fwd) / (x_mlg - x_nlg)
        frac_nlg_aft = (x_mlg - xcg_aft) / (x_mlg - x_nlg)

        # Tipback angle (for now assume that CG is along fuselage axis)
        alpha_tipback = np.arctan((x_mlg - xcg_aft) / (-z_lg))

        # Tailstrike angle
        alpha_tailstrike = np.arctan((z_tailstrike - z_lg) / (x_tailstrike - x_mlg))

        # Overturn angle
        sgl = (xcg_fwd - x_nlg) * y_mlg / np.sqrt((x_mlg - x_nlg) ** 2 + y_mlg**2)
        phi_overturn = np.arctan(-z_lg / sgl)

    else:

        # Add dummy data
        frac_nlg_fwd = None
        frac_nlg_aft = None
        alpha_tipback = None
        alpha_tailstrike = None
        phi_overturn = None

    # Update dictionary
    airplane["frac_nlg_fwd"] = frac_nlg_fwd
    airplane["frac_nlg_aft"] = frac_nlg_aft
    airplane["alpha_tipback"] = alpha_tipback
    airplane["alpha_tailstrike"] = alpha_tailstrike
    airplane["phi_overturn"] = phi_overturn

    return None


# ========================================
# AUXILIARY FUNCTIONS


def plot3d(airplane, figname="3dview.png", az1=45, az2=-135):
    """
    az1 and az2: degrees of azimuth and elevation for the 3d plot view
    """

    from matplotlib.patches import Ellipse
    import mpl_toolkits.mplot3d.art3d as art3d

    xr_w = airplane["xr_w"]
    zr_w = airplane["zr_w"]
    b_w = airplane["b_w"]

    tct_w = airplane["tct_w"]
    tcr_w = airplane["tcr_w"]

    cr_w = airplane["cr_w"]
    xt_w = airplane["xt_w"]
    yt_w = airplane["yt_w"]
    zt_w = airplane["zt_w"]
    ct_w = airplane["ct_w"]

    xr_h = airplane["xr_h"]
    zr_h = airplane["zr_h"]

    tcr_h = airplane["tcr_h"]
    tct_h = airplane["tct_h"]

    cr_h = airplane["cr_h"]
    xt_h = airplane["xt_h"]
    yt_h = airplane["yt_h"]
    zt_h = airplane["zt_h"]
    ct_h = airplane["ct_h"]
    b_h = airplane["b_h"]

    xr_v = airplane["xr_v"]
    zr_v = airplane["zr_v"]

    tcr_v = airplane["tcr_v"]
    tct_v = airplane["tct_v"]

    cr_v = airplane["cr_v"]
    xt_v = airplane["xt_v"]
    zt_v = airplane["zt_v"]
    ct_v = airplane["ct_v"]
    b_v = airplane["b_v"]

    L_f = airplane["L_f"]
    D_f = airplane["D_f"]
    x_n = airplane["x_n"]
    y_n = airplane["y_n"]
    z_n = airplane["z_n"]
    L_n = airplane["L_n"]
    D_n = airplane["D_n"]

    if "xcg_fwd" in airplane:
        xcg_fwd = airplane["xcg_fwd"]
        xcg_aft = airplane["xcg_aft"]
    else:
        xcg_fwd = None
        xcg_aft = None

    if "xnp" in airplane:
        xnp = airplane["xnp"]
    else:
        xnp = None

    x_nlg = airplane["x_nlg"]
    y_nlg = 0
    z_nlg = airplane["z_lg"]
    x_mlg = airplane["x_mlg"]
    y_mlg = airplane["y_mlg"]
    z_mlg = airplane["z_lg"]
    x_tailstrike = airplane["x_tailstrike"]
    z_tailstrike = airplane["z_tailstrike"]

    flap_type = airplane["flap_type"]
    b_flap_b_wing = airplane["b_flap_b_wing"]
    c_flap_c_wing = airplane["c_flap_c_wing"]

    slat_type = airplane["slat_type"]
    b_slat_b_wing = airplane["b_slat_b_wing"]
    c_slat_c_wing = airplane["c_slat_c_wing"]

    b_ail_b_wing = airplane["b_ail_b_wing"]
    c_ail_c_wing = airplane["c_ail_c_wing"]

    ### PLOT

    # fig = plt.figure(fignum,figsize=(20, 10))
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    # ax.set_aspect('equal')
    ax.plot(
        [xr_w, xt_w, xt_w + ct_w, xr_w + cr_w, xt_w + ct_w, xt_w, xr_w],
        [0.0, yt_w, yt_w, 0.0, -yt_w, -yt_w, 0.0],
        [
            zr_w + cr_w * tcr_w / 2,
            zt_w + ct_w * tct_w / 2,
            zt_w + ct_w * tct_w / 2,
            zr_w + cr_w * tcr_w / 2,
            zt_w + ct_w * tct_w / 2,
            zt_w + ct_w * tct_w / 2,
            zr_w + cr_w * tcr_w / 2,
        ],
        color="blue",
    )

    ax.plot(
        [xr_h, xt_h, xt_h + ct_h, xr_h + cr_h, xt_h + ct_h, xt_h, xr_h],
        [0.0, yt_h, yt_h, 0.0, -yt_h, -yt_h, 0.0],
        [
            zr_h + cr_h * tcr_h / 2,
            zt_h + ct_h * tct_h / 2,
            zt_h + ct_h * tct_h / 2,
            zr_h + cr_h * tcr_h / 2,
            zt_h + ct_h * tct_h / 2,
            zt_h + ct_h * tct_h / 2,
            zr_h + cr_h * tcr_h / 2,
        ],
        color="green",
    )

    ax.plot(
        [xr_v, xt_v, xt_v + ct_v, xr_v + cr_v, xr_v],
        [
            tcr_v * cr_v / 2,
            tct_v * ct_v / 2,
            tct_v * ct_v / 2,
            tcr_v * cr_v / 2,
            tcr_v * cr_v / 2,
        ],
        [zr_v, zt_v, zt_v, zr_v, zr_v],
        color="orange",
    )

    ax.plot(
        [xr_v, xt_v, xt_v + ct_v, xr_v + cr_v, xr_v],
        [
            -tcr_v * cr_v / 2,
            -tct_v * ct_v / 2,
            -tct_v * ct_v / 2,
            -tcr_v * cr_v / 2,
            -tcr_v * cr_v / 2,
        ],
        [zr_v, zt_v, zt_v, zr_v, zr_v],
        color="orange",
    )

    ax.plot([0.0, L_f], [0.0, 0.0], [0.0, 0.0])
    ax.plot([x_n, x_n + L_n], [y_n, y_n], [z_n, z_n])
    ax.plot([x_n, x_n + L_n], [-y_n, -y_n], [z_n, z_n])

    # Forward CG point
    if xcg_fwd is not None:
        ax.plot([xcg_fwd], [0.0], [0.0], "ko")

    # Rear CG point
    if xcg_aft is not None:
        ax.plot([xcg_aft], [0.0], [0.0], "ko")

    # Neutral point
    if xnp is not None:
        ax.plot([xnp], [0.0], [0.0], "x")

    # Define a parametrized fuselage by setting diameter
    # values along its axis
    # xx is non-dimensionalized by fuselage length
    # dd is non-dimensionalized by fuselage diameter
    xx = [0.0, 1.24 / 41.72, 3.54 / 41.72, 7.55 / 41.72, x_tailstrike / L_f, 1.0]
    hh = [0.0, 2.27 / 4.0, 3.56 / 4.0, 1.0, 1.0, 1.07 / 4.0]
    ww = [0.0, 1.83 / 4.0, 3.49 / 4.0, 1.0, 1.0, 0.284 / 4]
    num_tot_ell = 50  # Total number of ellipses

    # Loop over every section
    for ii in range(len(xx) - 1):

        # Define number of ellipses based on the section length
        num_ell = int((xx[ii + 1] - xx[ii]) * num_tot_ell) + 1

        # Define arrays of dimensional positions, heights and widths
        # for the current section
        xdim = np.linspace(xx[ii], xx[ii + 1], num_ell) * L_f
        hdim = np.linspace(hh[ii], hh[ii + 1], num_ell) * D_f
        wdim = np.linspace(ww[ii], ww[ii + 1], num_ell) * D_f

        # Loop over every ellipse
        for xc, hc, wc in zip(xdim, hdim, wdim):
            p = Ellipse(
                (0, 0), wc, hc, angle=0, facecolor="none", edgecolor="k", lw=1.0
            )
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=xc, zdir="x")

    # ____________________________________________________________
    #                                                            \
    # MLG / NLG

    # Check if LG is activated
    d_lg = 0
    if x_nlg is not None:

        # Make landing gear dimensions based on the fuselage
        w_lg = 0.05 * D_f
        d_lg = 4 * w_lg

        mlg_len = np.linspace(y_mlg - w_lg / 2, y_mlg + w_lg / 2, 2)
        nlg_len = np.linspace(y_nlg - w_lg / 2, y_nlg + w_lg / 2, 2)

        for i in range(len(mlg_len)):
            p = Ellipse(
                (x_mlg, z_mlg),
                d_lg,
                d_lg,
                angle=0,
                facecolor="gray",
                edgecolor="k",
                lw=2,
            )
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=mlg_len[i], zdir="y")

            p = Ellipse(
                (x_mlg, z_mlg),
                d_lg,
                d_lg,
                angle=0,
                facecolor="gray",
                edgecolor="k",
                lw=2,
            )
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=-mlg_len[i], zdir="y")

            # NLG
            p = Ellipse(
                (x_nlg, z_nlg),
                d_lg,
                d_lg,
                angle=0,
                facecolor="gray",
                edgecolor="k",
                lw=1.5,
            )
            ax.add_patch(p)
            art3d.pathpatch_2d_to_3d(p, z=nlg_len[i], zdir="y")

    # Nacelle
    nc_len = np.linspace(x_n, x_n + L_n, 11)
    for i in range(len(nc_len)):
        p = Ellipse(
            (y_n, z_n), D_n, D_n, angle=0, facecolor="none", edgecolor="orange", lw=1.0
        )
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

        # Inner wall
        # p = Ellipse((y_n, z_n), D_n*0.8, D_n*0.8, angle=0,\
        # facecolor = 'none', edgecolor = 'k', lw=.1)
        # ax.add_patch(p)
        # art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

        p = Ellipse(
            (-y_n, z_n), D_n, D_n, angle=0, facecolor="none", edgecolor="orange", lw=1.0
        )
        ax.add_patch(p)
        art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

        # Inner wall
        # p = Ellipse((-y_n, z_n), D_n*0.8, D_n*0.8, angle=0, \
        # facecolor = 'none', edgecolor = 'k', lw=.1)
        # ax.add_patch(p)
        # art3d.pathpatch_2d_to_3d(p, z=nc_len[i], zdir="x")

    # Aileron
    ail_tip_margin = 0.02  # Margem entre flap e aileron em % de b_w

    # Spanwise positions (root and tip)
    yr_a = (1.0 - (ail_tip_margin + b_ail_b_wing)) * b_w / 2
    yt_a = (1.0 - (ail_tip_margin)) * b_w / 2

    cr_a = lin_interp(0, b_w / 2, cr_w, ct_w, yr_a) * c_ail_c_wing
    ct_a = lin_interp(0, b_w / 2, cr_w, ct_w, yt_a) * c_ail_c_wing

    # To find the longitudinal position of the aileron LE, we find the TE position first
    # then we subtract the chord
    xr_a = lin_interp(0, b_w / 2, xr_w + cr_w, xt_w + ct_w, yr_a) - cr_a
    xt_a = lin_interp(0, b_w / 2, xr_w + cr_w, xt_w + ct_w, yt_a) - ct_a

    zr_a = lin_interp(0, b_w / 2, zr_w, zt_w, yr_a)
    zt_a = lin_interp(0, b_w / 2, zr_w, zt_w, yt_a)

    # Airfoil thickness at aileron location
    tcr_a = lin_interp(0, b_w / 2, tcr_w, tct_w, yr_a)
    tct_a = lin_interp(0, b_w / 2, tcr_w, tct_w, yt_a)

    ax.plot(
        [xr_a, xt_a, xt_a + ct_a, xr_a + cr_a, xr_a],
        [yr_a, yt_a, yt_a, yr_a, yr_a],
        [
            zr_a + cr_a * tcr_a / 2 / c_ail_c_wing,
            zt_a + ct_a * tct_a / 2 / c_ail_c_wing,
            zt_a + ct_a * tct_a / 2 / c_ail_c_wing,
            zr_a + cr_a * tcr_a / 2 / c_ail_c_wing,
            zr_a + cr_a * tcr_a / 2 / c_ail_c_wing,
        ],
        lw=1,
        color="green",
    )

    ax.plot(
        [xr_a, xt_a, xt_a + ct_a, xr_a + cr_a, xr_a],
        [-yr_a, -yt_a, -yt_a, -yr_a, -yr_a],
        [
            zr_a + cr_a * tcr_a / 2 / c_ail_c_wing,
            zt_a + ct_a * tct_a / 2 / c_ail_c_wing,
            zt_a + ct_a * tct_a / 2 / c_ail_c_wing,
            zr_a + cr_a * tcr_a / 2 / c_ail_c_wing,
            zr_a + cr_a * tcr_a / 2 / c_ail_c_wing,
        ],
        lw=1,
        color="green",
    )

    # Slat
    if slat_type is not None:

        # slat_tip_margin = 0.02  # Margem da ponta como % da b_w
        # slat_root_margin = 0.12 # Margem da raiz como % da b_w
        # hist_c_s = 0.25        # Corda do Flap
        # hist_b_s = 1 - slat_root_margin - slat_tip_margin

        # Spanwise positions (root and tip)
        yr_s = D_f / 2
        yt_s = b_slat_b_wing * b_w / 2

        cr_s = lin_interp(0, b_w / 2, cr_w, ct_w, yr_s) * c_slat_c_wing
        ct_s = lin_interp(0, b_w / 2, cr_w, ct_w, yt_s) * c_slat_c_wing

        # Find the longitudinal position of the slat LE
        xr_s = lin_interp(0, b_w / 2, xr_w, xt_w, yr_s)
        xt_s = lin_interp(0, b_w / 2, xr_w, xt_w, yt_s)

        zr_s = lin_interp(0, b_w / 2, zr_w, zt_w, yr_s)
        zt_s = lin_interp(0, b_w / 2, zr_w, zt_w, yt_s)

        # Airfoil thickness at slat location
        tcr_s = lin_interp(0, b_w / 2, tcr_w, tct_w, yr_s)
        tct_s = lin_interp(0, b_w / 2, tcr_w, tct_w, yt_s)

        ax.plot(
            [xr_s, xt_s, xt_s + ct_s, xr_s + cr_s, xr_s],
            [yr_s, yt_s, yt_s, yr_s, yr_s],
            [
                zr_s + cr_s * tcr_s / 2 / c_slat_c_wing,
                zt_s + ct_s * tct_s / 2 / c_slat_c_wing,
                zt_s + ct_s * tct_s / 2 / c_slat_c_wing,
                zr_s + cr_s * tcr_s / 2 / c_slat_c_wing,
                zr_s + cr_s * tcr_s / 2 / c_slat_c_wing,
            ],
            lw=1,
            color="m",
        )

        ax.plot(
            [xr_s, xt_s, xt_s + ct_s, xr_s + cr_s, xr_s],
            [-yr_s, -yt_s, -yt_s, -yr_s, -yr_s],
            [
                zr_s + cr_s * tcr_s / 2 / c_slat_c_wing,
                zt_s + ct_s * tct_s / 2 / c_slat_c_wing,
                zt_s + ct_s * tct_s / 2 / c_slat_c_wing,
                zr_s + cr_s * tcr_s / 2 / c_slat_c_wing,
                zr_s + cr_s * tcr_s / 2 / c_slat_c_wing,
            ],
            lw=1,
            color="m",
        )

    # Flap outboard
    if flap_type is not None:

        # Spanwise positions (root and tip)
        yr_f = D_f / 2
        yt_f = b_flap_b_wing * b_w / 2

        cr_f = lin_interp(0, b_w / 2, cr_w, ct_w, yr_f) * c_flap_c_wing
        ct_f = lin_interp(0, b_w / 2, cr_w, ct_w, yt_f) * c_flap_c_wing

        # To find the longitudinal position of the flap LE, we find the TE position first
        # then we subtract the chord
        xr_f = lin_interp(0, b_w / 2, xr_w + cr_w, xt_w + ct_w, yr_f) - cr_f
        xt_f = lin_interp(0, b_w / 2, xr_w + cr_w, xt_w + ct_w, yt_f) - ct_f

        zr_f = lin_interp(0, b_w / 2, zr_w, zt_w, yr_f)
        zt_f = lin_interp(0, b_w / 2, zr_w, zt_w, yt_f)

        # Airfoil thickness at flap location
        tcr_f = lin_interp(0, b_w / 2, tcr_w, tct_w, yr_f)
        tct_f = lin_interp(0, b_w / 2, tcr_w, tct_w, yt_f)

        ax.plot(
            [xr_f, xt_f, xt_f + ct_f, xr_f + cr_f, xr_f],
            [yr_f, yt_f, yt_f, yr_f, yr_f],
            [
                zr_f + cr_f * tcr_f / 2 / c_flap_c_wing,
                zt_f + ct_f * tct_f / 2 / c_flap_c_wing,
                zt_f + ct_f * tct_f / 2 / c_flap_c_wing,
                zr_f + cr_f * tcr_f / 2 / c_flap_c_wing,
                zr_f + cr_f * tcr_f / 2 / c_flap_c_wing,
            ],
            lw=1,
            color="r",
        )

        ax.plot(
            [xr_f, xt_f, xt_f + ct_f, xr_f + cr_f, xr_f],
            [-yr_f, -yt_f, -yt_f, -yr_f, -yr_f],
            [
                zr_f + cr_f * tcr_f / 2 / c_flap_c_wing,
                zt_f + ct_f * tct_f / 2 / c_flap_c_wing,
                zt_f + ct_f * tct_f / 2 / c_flap_c_wing,
                zr_f + cr_f * tcr_f / 2 / c_flap_c_wing,
                zr_f + cr_f * tcr_f / 2 / c_flap_c_wing,
            ],
            lw=1,
            color="r",
        )

    # Elevator
    ele_tip_margin = 0.1  # Margem do profundor para a ponta
    ele_root_margin = 0.1  # Margem do profundor para a raiz
    hist_b_e = 1 - ele_root_margin - ele_tip_margin
    hist_c_e = 0.25

    ct_e_loc = (1 - ele_tip_margin) * (ct_h - cr_h) + cr_h
    cr_e_loc = (1 - hist_b_e - ele_tip_margin) * (ct_h - cr_h) + cr_h

    ct_e = ct_e_loc * hist_c_e
    cr_e = cr_e_loc * hist_c_e

    xr_e = (
        (1 - hist_b_e - ele_tip_margin) * (xt_h - xr_h)
        + xr_h
        + cr_e_loc * (1 - hist_c_e)
    )
    xt_e = (1 - ele_tip_margin) * (xt_h - xr_h) + xr_h + ct_e_loc * (1 - hist_c_e)

    yr_e = (1 - hist_b_e - ele_tip_margin) * b_h / 2
    yt_e = (1 - ele_tip_margin) * b_h / 2

    zr_e = (1 - hist_b_e - ele_tip_margin) * (zt_h - zr_h) + zr_h
    zt_e = (1 - ele_tip_margin) * (zt_h - zr_h) + zr_h

    ax.plot(
        [xr_e, xt_e, xt_e + ct_e, xr_e + cr_e, xr_e],
        [yr_e, yt_e, yt_e, yr_e, yr_e],
        [zr_e, zt_e, zt_e, zr_e, zr_e],
        lw=1,
        color="g",
    )

    ax.plot(
        [xr_e, xt_e, xt_e + ct_e, xr_e + cr_e, xr_e],
        [-yr_e, -yt_e, -yt_e, -yr_e, -yr_e],
        [zr_e, zt_e, zt_e, zr_e, zr_e],
        lw=1,
        color="g",
    )

    # Rudder
    ver_base_margin = 0.1  # Local da base % de b_v
    ver_tip_margin1 = 0.1  # Local da base % de b_v
    ver_tip_margin = 1 - ver_tip_margin1  # Local do topo % de b_v
    hist_c_v = 0.32

    cr_v_loc = ver_base_margin * (ct_v - cr_v) + cr_v
    ct_v_loc = ver_tip_margin * (ct_v - cr_v) + cr_v

    cr_v2 = cr_v_loc * hist_c_v
    ct_v2 = ct_v_loc * hist_c_v

    xr_v2 = ver_base_margin * (xt_v - xr_v) + xr_v + cr_v_loc * (1 - hist_c_v)
    xt_v2 = ver_tip_margin * (xt_v - xr_v) + xr_v + ct_v_loc * (1 - hist_c_v)

    zr_v2 = ver_base_margin * (zt_v - zr_v) + zr_v
    zt_v2 = ver_tip_margin * (zt_v - zr_v) + zr_v

    ax.plot(
        [xr_v2, xt_v2, xt_v2 + ct_v2, xr_v2 + cr_v2, xr_v2],
        [
            tcr_v * cr_v_loc / 2,
            tct_v * ct_v_loc / 2,
            tct_v * ct_v_loc / 2,
            tcr_v * cr_v_loc / 2,
            tcr_v * cr_v_loc / 2,
        ],
        [zr_v2, zt_v2, zt_v2, zr_v2, zr_v2],
        color="orange",
    )

    ax.plot(
        [xr_v2, xt_v2, xt_v2 + ct_v2, xr_v2 + cr_v2, xr_v2],
        [
            -tcr_v * cr_v_loc / 2,
            -tct_v * ct_v_loc / 2,
            -tct_v * ct_v_loc / 2,
            -tcr_v * cr_v_loc / 2,
            -tcr_v * cr_v_loc / 2,
        ],
        [zr_v2, zt_v2, zt_v2, zr_v2, zr_v2],
        color="orange",
    )

    # _______ONLY FRONT VIEW_______

    # Wing Lower
    # ------------------------------
    ax.plot(
        [xr_w, xt_w, xt_w + ct_w, xr_w + cr_w, xt_w + ct_w, xt_w, xr_w],
        [0.0, yt_w, yt_w, 0.0, -yt_w, -yt_w, 0.0],
        [
            zr_w - tcr_w * cr_w / 2,
            zt_w - tct_w * ct_w / 2,
            zt_w - tct_w * ct_w / 2,
            zr_w - tcr_w * cr_w / 2,
            zt_w - tct_w * ct_w / 2,
            zt_w - tct_w * ct_w / 2,
            zr_w - tcr_w * cr_w / 2,
        ],
        color="blue",
    )

    ax.plot(
        [xr_w, xr_w],
        [0.0, 0.0],
        [zr_w - tcr_w * cr_w / 2, zr_w + tcr_w * cr_w / 2],
        color="blue",
    )
    ax.plot(
        [xr_w + cr_w, xr_w + cr_w],
        [0.0, 0.0],
        [zr_w - tcr_w * cr_w / 2, zr_w + tcr_w * cr_w / 2],
        color="blue",
    )

    ax.plot(
        [xt_w, xt_w],
        [yt_w, yt_w],
        [zt_w - tct_w * ct_w / 2, zt_w + tct_w * ct_w / 2],
        color="blue",
    )
    ax.plot(
        [xt_w + ct_w, xt_w + ct_w],
        [yt_w, yt_w],
        [zt_w - tct_w * ct_w / 2, zt_w + tct_w * ct_w / 2],
        color="blue",
    )

    ax.plot(
        [xt_w, xt_w],
        [-yt_w, -yt_w],
        [zt_w - tct_w * ct_w / 2, zt_w + tct_w * ct_w / 2],
        color="blue",
    )
    ax.plot(
        [xt_w + ct_w, xt_w + ct_w],
        [-yt_w, -yt_w],
        [zt_w - tct_w * ct_w / 2, zt_w + tct_w * ct_w / 2],
        color="blue",
    )

    # ------------------------------

    # HT Lower
    # ------------------------------
    ax.plot(
        [xr_h, xt_h, xt_h + ct_h, xr_h + cr_h, xt_h + ct_h, xt_h, xr_h],
        [0.0, yt_h, yt_h, 0.0, -yt_h, -yt_h, 0.0],
        [
            zr_h - tcr_h * cr_h / 2,
            zt_h - tct_h * ct_h / 2,
            zt_h - tct_h * ct_h / 2,
            zr_h - tcr_h * cr_h / 2,
            zt_h - tct_h * ct_h / 2,
            zt_h - tct_h * ct_h / 2,
            zr_h - tcr_h * cr_h / 2,
        ],
        color="green",
    )

    ax.plot(
        [xr_h, xr_h],
        [0.0, 0.0],
        [zr_h - tcr_h * cr_h / 2, zr_h + tcr_h * cr_h / 2],
        color="green",
    )
    ax.plot(
        [xr_h + cr_h, xr_h + cr_h],
        [0.0, 0.0],
        [zr_h - tcr_h * cr_h / 2, zr_h + tcr_h * cr_h / 2],
        color="green",
    )

    ax.plot(
        [xt_h, xt_h],
        [yt_h, yt_h],
        [zt_h - tct_h * ct_h / 2, zt_h + tct_h * ct_h / 2],
        color="green",
    )
    ax.plot(
        [xt_h + ct_h, xt_h + ct_h],
        [yt_h, yt_h],
        [zt_h - tct_h * ct_h / 2, zt_h + tct_h * ct_h / 2],
        color="green",
    )

    ax.plot(
        [xt_h, xt_h],
        [-yt_h, -yt_h],
        [zt_h - tct_h * ct_h / 2, zt_h + tct_h * ct_h / 2],
        color="green",
    )
    ax.plot(
        [xt_h + ct_h, xt_h + ct_h],
        [-yt_h, -yt_h],
        [zt_h - tct_h * ct_h / 2, zt_h + tct_h * ct_h / 2],
        color="green",
    )

    # Slat Lower
    # ------------------------------
    if slat_type is not None:
        ax.plot(
            [xr_s, xt_s, xt_s + ct_s, xr_s + cr_s, xr_s],
            [yr_s, yt_s, yt_s, yr_s, yr_s],
            [
                zr_s - tcr_s * cr_s / 2 / c_slat_c_wing,
                zt_s - tct_s * ct_s / 2 / c_slat_c_wing,
                zt_s - tct_s * ct_s / 2 / c_slat_c_wing,
                zr_s - tcr_s * cr_s / 2 / c_slat_c_wing,
                zr_s - tcr_s * cr_s / 2 / c_slat_c_wing,
            ],
            lw=1,
            color="m",
        )

        ax.plot(
            [xr_s, xt_s, xt_s + ct_s, xr_s + cr_s, xr_s],
            [-yr_s, -yt_s, -yt_s, -yr_s, -yr_s],
            [
                zr_s - tcr_s * cr_s / 2 / c_slat_c_wing,
                zt_s - tct_s * ct_s / 2 / c_slat_c_wing,
                zt_s - tct_s * ct_s / 2 / c_slat_c_wing,
                zr_s - tcr_s * cr_s / 2 / c_slat_c_wing,
                zr_s - tcr_s * cr_s / 2 / c_slat_c_wing,
            ],
            lw=1,
            color="m",
        )
    # ------------------------------

    # Flap Lower
    # ------------------------------
    if flap_type is not None:
        ax.plot(
            [xr_f, xt_f, xt_f + ct_f, xr_f + cr_f, xr_f],
            [yr_f, yt_f, yt_f, yr_f, yr_f],
            [
                zr_f - tcr_f * cr_f / 2 / c_flap_c_wing,
                zt_f - tct_f * ct_f / 2 / c_flap_c_wing,
                zt_f - tct_f * ct_f / 2 / c_flap_c_wing,
                zr_f - tcr_f * cr_f / 2 / c_flap_c_wing,
                zr_f - tcr_f * cr_f / 2 / c_flap_c_wing,
            ],
            lw=1,
            color="r",
        )

        ax.plot(
            [xr_f, xt_f, xt_f + ct_f, xr_f + cr_f, xr_f],
            [-yr_f, -yt_f, -yt_f, -yr_f, -yr_f],
            [
                zr_f - tcr_f * cr_f / 2 / c_flap_c_wing,
                zt_f - tct_f * ct_f / 2 / c_flap_c_wing,
                zt_f - tct_f * ct_f / 2 / c_flap_c_wing,
                zr_f - tcr_f * cr_f / 2 / c_flap_c_wing,
                zr_f - tcr_f * cr_f / 2 / c_flap_c_wing,
            ],
            lw=1,
            color="r",
        )
    # ------------------------------

    # Aleron Lower
    # ------------------------------
    ax.plot(
        [xr_a, xt_a, xt_a + ct_a, xr_a + cr_a, xr_a],
        [yr_a, yt_a, yt_a, yr_a, yr_a],
        [
            zr_a - tcr_a * cr_a / 2 / c_ail_c_wing,
            zt_a - tct_a * ct_a / 2 / c_ail_c_wing,
            zt_a - tct_a * ct_a / 2 / c_ail_c_wing,
            zr_a - tcr_a * cr_a / 2 / c_ail_c_wing,
            zr_a - tcr_a * cr_a / 2 / c_ail_c_wing,
        ],
        lw=1,
        color="green",
    )

    ax.plot(
        [xr_a, xt_a, xt_a + ct_a, xr_a + cr_a, xr_a],
        [-yr_a, -yt_a, -yt_a, -yr_a, -yr_a],
        [
            zr_a - tcr_a * cr_a / 2 / c_ail_c_wing,
            zt_a - tct_a * ct_a / 2 / c_ail_c_wing,
            zt_a - tct_a * ct_a / 2 / c_ail_c_wing,
            zr_a - tcr_a * cr_a / 2 / c_ail_c_wing,
            zr_a - tcr_a * cr_a / 2 / c_ail_c_wing,
        ],
        lw=1,
        color="green",
    )
    # ------------------------------

    # Avoiding blanketing the rudder
    ax.plot(
        [xr_h, xr_h + b_v / np.tan(60 * np.pi / 180)],
        [0.0, 0.0],
        [zr_h, zr_h + b_v],
        "k--",
    )

    ax.plot(
        [xr_h + cr_h, xr_h + 0.6 * b_v / np.tan(30 * np.pi / 180) + cr_h],
        [0.0, 0.0],
        [zr_h, zr_h + 0.6 * b_v],
        "k--",
    )

    # Water Spray
    if x_nlg is not None:
        ax.plot(
            [x_nlg, x_nlg + 10 / np.tan(22 * np.pi / 180)],
            [0.0, 4.8],
            [0.0, 0.0],
            "k--",
        )

        # Water Spray
        ax.plot(
            [x_nlg, x_nlg + 10 / np.tan(22 * np.pi / 180)],
            [0.0, -4.8],
            [0.0, 0.0],
            "k--",
        )

    # Create cubic bounding box to simulate equal aspect ratio
    # First create o list of possible critical points along each coordinate
    X = np.array(
        [
            0,
            xr_w,
            xt_h + ct_h,
            xt_v + ct_v,
            L_f,
            xr_h + b_v / np.tan(60 * np.pi / 180),
            xr_h + 0.6 * b_v / np.tan(30 * np.pi / 180) + cr_h,
        ]
    )
    Y = np.array([-yt_w, yt_w])
    Z = np.array([-D_f / 2, zt_w, zt_h, zt_v, z_mlg - d_lg / 2, zr_h + b_v])
    max_range = np.array(
        [X.max() - X.min(), Y.max() - Y.min(), Z.max() - Z.min()]
    ).max()
    Xb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][0].flatten() + 0.5 * (
        X.max() + X.min()
    )
    Yb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][1].flatten() + 0.5 * (
        Y.max() + Y.min()
    )
    Zb = 0.5 * max_range * np.mgrid[-1:2:2, -1:2:2, -1:2:2][2].flatten() + 0.5 * (
        Z.max() + Z.min()
    )

    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], "w")

    # ax.set_box_aspect((1, 1, 1))
    ax.view_init(az1, az2)

    fig.savefig(figname, dpi=300)

    plt.show()


# ----------------------------------------


def atmosphere(z, Tba=288.15):

    """
    Funçao que retorna a Temperatura, Pressao e Densidade para uma determinada
    altitude z [m]. Essa funçao usa o modelo padrao de atmosfera para a
    temperatura no solo de Tba.
    """

    # Zbase (so para referencia)
    # 0 11019.1 20063.1 32161.9 47350.1 50396.4

    # DEFINING CONSTANTS
    # Earth radius
    r = 6356766
    # gravity
    g0 = 9.80665
    # air gas constant
    R = 287.05287
    # layer boundaries
    Ht = [0, 11000, 20000, 32000, 47000, 50000]
    # temperature slope in each layer
    A = [-6.5e-3, 0, 1e-3, 2.8e-3, 0]
    # pressure at the base of each layer
    pb = [101325, 22632, 5474.87, 868.014, 110.906]
    # temperature at the base of each layer
    Tstdb = [288.15, 216.65, 216.65, 228.65, 270.65]
    # temperature correction
    Tb = Tba - Tstdb[0]
    # air viscosity
    mi0 = 18.27e-6  # [Pa s]
    T0 = 291.15  # [K]
    C = 120  # [K]

    # geopotential altitude
    H = r * z / (r + z)

    # selecting layer
    if H < Ht[0]:
        raise ValueError("Under sealevel")
    elif H <= Ht[1]:
        i = 0
    elif H <= Ht[2]:
        i = 1
    elif H <= Ht[3]:
        i = 2
    elif H <= Ht[4]:
        i = 3
    elif H <= Ht[5]:
        i = 4
    else:
        raise ValueError("Altitude beyond model boundaries")

    # Calculating temperature
    T = Tstdb[i] + A[i] * (H - Ht[i]) + Tb

    # Calculating pressure
    if A[i] == 0:
        p = pb[i] * np.exp(-g0 * (H - Ht[i]) / R / (Tstdb[i] + Tb))
    else:
        p = pb[i] * (T / (Tstdb[i] + Tb)) ** (-g0 / A[i] / R)

    # Calculating density
    rho = p / R / T

    # Calculating viscosity with Sutherland's Formula
    mi = mi0 * (T0 + C) / (T + C) * (T / T0) ** (1.5)

    return T, p, rho, mi


# ----------------------------------------


def geo_change_sweep(x, y, sweep_x, panel_length, chord_root, chord_tip):

    """
    This function converts sweep computed at chord fraction x into
    sweep measured at chord fraction y
    (x and y should be between 0 (leading edge) and 1 (trailing edge).
    """

    sweep_y = sweep_x + np.arctan((x - y) * (chord_root - chord_tip) / panel_length)

    return sweep_y


# ----------------------------------------


def Cf_calc(Mach, altitude, length, rugosity, k_lam, Tba=288.15):
    """
    This function computes the flat plate friction coefficient
    for a given Reynolds number while taking transition into account

    k_lam: float -> Fraction of the length (from 0 to 1) where
                    transition occurs
    """

    # Dados atmosféricos
    T, p, rho, mi = atmosphere(altitude, Tba)

    # Velocidade
    v = np.sqrt(gamma_air * R_air * T) * Mach

    # Reynolds na transição
    Re_conv = rho * v * k_lam * length / mi
    Re_rug = 38.21 * (k_lam * length / rugosity) ** 1.053
    Re_trans = min(Re_conv, Re_rug)

    # Reynolds no fim
    Re_conv = rho * v * length / mi
    Re_rug = 38.21 * (length / rugosity) ** 1.053
    Re_fim = min(Re_conv, Re_rug)

    # Coeficientes de fricção
    # Laminar na transição
    Cf1 = 1.328 / np.sqrt(Re_trans)

    # Turbulento na transição
    Cf2 = 0.455 / (np.log10(Re_trans) ** 2.58 * (1 + 0.144 * Mach**2) ** 0.65)

    # Turbulento no fim
    Cf3 = 0.455 / (np.log10(Re_fim) ** 2.58 * (1 + 0.144 * Mach**2) ** 0.65)

    # Média
    Cf = (Cf1 - Cf2) * k_lam + Cf3

    return Cf


# ----------------------------------------


def FF_surface(Mach, tcr, tct, sweep, b, cr, ct, cm, x_c_max_tc=0.4):
    """
    This function computes the form factor for lifting surfaces

    INPUTS

    tcr: float -> Thickness/chord ratio at the root
    tct: float -> Thickness/chord ratio at the tip
    sweep: float -> Quarter-chord sweep angle [rad]
    b: float -> Wing span (considering both sides. Double this value for vertical tails if necessary)
    cr: float -> Root chord
    ct: float -> Tip chord
    cm: float -> Mean aerodynamic chord
    x_c_max_tc: float -> Chord fraction with maximum thickness
    """

    # Average chord fraction
    t_c = (tcr + tct) / 2

    # Sweep at maximum thickness position
    sweep_maxtc = geo_change_sweep(0.25, x_c_max_tc, sweep, b / 2, cr, ct)

    # Form factor
    FF = (
        1.34
        * Mach**0.18
        * np.cos(sweep_maxtc) ** 0.28
        * (1 + 0.6 * t_c / x_c_max_tc + 100 * (t_c) ** 4)
    )

    return FF


# ----------------------------------------


def lin_interp(x0, x1, y0, y1, x):
    """
    Linear interpolation function
    """

    y = y0 + (y1 - y0) * (x - x0) / (x1 - x0)

    return y


# ----------------------------------------


def turbojet_model(Mach, altitude, data):

    """
    Model implemented by Marcelo Yuri Sampaio de Freitas
    during his undergraduate thesis:
    "Projeto conceitual de alvo aéreo manobrável baseado na propulsão de turbojato"
    Instituto Tecnológico de Aeronáutica, 2020.
    http://www.bdita.bibl.ita.br/TGsDigitais/lista_resumo.php?num_tg=77532

    This function returns the TSFC and maximum thrust of the engine

    data should be a dictionary with the following fields:
    data['p02_p01']: compression ratio
    data['T03']: turbine inlet temperature
    data['n_comp']: compressor efficiency
    data['n_turbine']: turbine efficiency
    data['n_intake']: intake efficiency
    data['n_isen_noz']: insentropic nozzle efficiency
    data['n_mec']: mechanical efficiency
    data['n_comb']: combustion efficiency
    data['delta_pb']: pressure loss at the combustion chamber
    data['A5']: nozzle area
    """

    # Thermodynamic parameters
    Cp_air = 1005  # J/kg*K
    gamma_comb = 1.333
    Cp_comb = 1148  # J/kg*K
    LHV = 43.1e6  # J/kg fuel calorific power

    # Air properties
    [Ta, pa, rho, mi] = atmosphere(altitude, 288.15)
    sound_speed = np.sqrt(gamma_air * R_air * Ta)

    # Generic engine parameters
    p02_p01 = data["p02_p01"]
    T03 = data["T03"]
    n_comp = data["n_comp"]
    n_turbine = data["n_turbine"]
    n_intake = data["n_intake"]
    n_isen_noz = data["n_isen_noz"]
    n_mec = data["n_mec"]
    n_comb = data["n_comb"]
    delta_pb = data["delta_pb"]  # Pressure loss in combustion chamber

    # O requisito pode ser alterado para fluxo de massa
    A5 = data["A5"]
    # m_dot = 7.823

    # Intake
    T_estag = ((Mach * sound_speed) ** 2) / (2 * Cp_air)
    T01 = Ta + T_estag
    p01_pa = (1 + n_intake * T_estag / Ta) ** (gamma_air / (gamma_air - 1))
    p01 = p01_pa * pa

    # Compressor
    p02 = p02_p01 * p01
    T02 = T01 + (T01 / n_comp) * (p02_p01 ** ((gamma_air - 1) / gamma_air) - 1)

    # Combustion
    p03 = p02 * (1 - (delta_pb / 100))
    T03 = T03

    # Turbine
    deltaT = Cp_air * (T02 - T01) / Cp_comb / n_mec
    T04 = T03 - deltaT
    T04_isen = T03 - deltaT / n_turbine
    p04 = p03 * (T04_isen / T03) ** (gamma_comb / (gamma_comb - 1))

    # Exhaust nozzle
    p04_pa = p04 / pa
    p04_pc = 1 / (
        (1 - (1 / n_isen_noz) * ((gamma_comb - 1) / (gamma_comb + 1)))
        ** (gamma_comb / (gamma_comb - 1))
    )

    if p04_pa >= p04_pc:

        # Choked exhaust
        T05 = (2 / (gamma_comb + 1)) * T04
        p05 = p04 / p04_pc
        rho5 = p05 / (R_air * T05)
        C5 = np.sqrt(gamma_comb * R_air * T05)
        A5s = 1 / (rho5 * C5)
        Fs = (C5 - Mach * sound_speed) + A5s * (p05 - pa)
        m_dot = A5 / A5s
        F = Fs * m_dot  # Thrust
        f_comb = (Cp_comb * T03 - Cp_air * T02) / (LHV * n_comb - Cp_comb * T03)
        m_dot_fuel = f_comb * m_dot
        SFC = m_dot_fuel / F  # kg/s.N

    else:

        # Exhaust is not choked
        T05 = T04 - n_isen_noz * T04 * (
            1 - (1 / p04_pa) ** ((gamma_comb - 1) / gamma_comb)
        )
        C5 = np.sqrt(2 * Cp_comb * (T04 - T05))
        p05 = pa
        rho5 = p05 / (R_air * T05)
        A5s = 1 / (rho5 * C5)
        Fs = (C5 - Mach * sound_speed) + A5s * (p05 - pa)
        m_dot = A5 / A5s
        F = Fs * m_dot  # Thrust
        f_comb = (Cp_comb * T03 - Cp_air * T02) / (LHV * n_comb - Cp_comb * T03)
        m_dot_fuel = f_comb * m_dot
        SFC = m_dot_fuel / F  # kg/s.N
        # sigma = rho/1.225

    # Convert to 1/s
    C = SFC * 9.81

    return C, F


# ----------------------------------------


def standard_airplane(name="fokker100"):
    """
    The standard parameters refer to the Fokker 100, but they could be redefined for
    any new aircraft.
    """

    if name == "fokker100":

        airplane = {
            "type": "transport",  # Can be 'transport', 'fighter', or 'general'
            "S_w": 93.5,  # Wing area [m2]
            "AR_w": 8.43,  # Wing aspect ratio
            "taper_w": 0.235,  # Wing taper ratio
            "sweep_w": 17.45 * np.pi / 180,  # Wing sweep [rad]
            "dihedral_w": 5 * np.pi / 180,  # Wing dihedral [rad]
            "xr_w": 13.5,  # Longitudinal position of the wing (with respect to the fuselage nose) [m]
            "zr_w": -1.5,  # Vertical position of the wing (with respect to the fuselage nose) [m]
            "tcr_w": 0.123,  # t/c of the root section of the wing
            "tct_w": 0.096,  # t/c of the tip section of the wing
            "Cht": 1.14,  # Horizontal tail volume coefficient
            "Lc_h": 4.83,  # Non-dimensional lever of the horizontal tail (lever/wing_mac)
            "AR_h": 4.64,  # HT aspect ratio
            "taper_h": 0.39,  # HT taper ratio
            "sweep_h": 26 * np.pi / 180,  # HT sweep [rad]
            "dihedral_h": 2 * np.pi / 180,  # HT dihedral [rad]
            "zr_h": 6.359,  # Vertical position of the HT [m]
            "tcr_h": 0.1,  # t/c of the root section of the HT
            "tct_h": 0.1,  # t/c of the tip section of the HT
            "eta_h": 1.0,  # Dynamic pressure factor of the HT
            "Cvt": 0.06,  # Vertical tail volume coefficient
            "Lb_v": 0.55,  # Non-dimensional lever of the vertical tail (lever/wing_span)
            "AR_v": 1.27,  # VT aspect ratio
            "taper_v": 0.74,  # VT taper ratio
            "sweep_v": 41 * np.pi / 180,  # VT sweep [rad]
            "zr_v": 0.0,  # Vertical position of the VT [m]
            "tcr_v": 0.1,  # t/c of the root section of the VT
            "tct_v": 0.1,  # t/c of the tip section of the VT
            "L_f": 32.5,  # Fuselage length [m]
            "D_f": 3.3,  # Fuselage diameter [m]
            "x_n": 23.2,  # Longitudinal position of the nacelle frontal face [m]
            "y_n": 2.6,  # Lateral position of the nacelle centerline [m]
            "z_n": 0.0,  # Vertical position of the nacelle centerline [m]
            "L_n": 4.3,  # Nacelle length [m]
            "D_n": 1.5,  # Nacelle diameter [m]
            "n_engines": 2,  # Number of engines
            "n_engines_under_wing": 0,  # Number of engines installed under the wing
            "engine": {
                "model": "Howe turbofan",  # Check engineTSFC function for options
                "BPR": 3.04,  # Engine bypass ratio
                "Cbase": 0.57 / 3600,
            },
            "x_nlg": 3.6,  # Longitudinal position of the nose landing gear [m]
            "x_mlg": 17.8,  # Longitudinal position of the main landing gear [m]
            "y_mlg": 2.47,  # Lateral position of the main landing gear [m]
            "z_lg": -2.0,  # Vertical position of the landing gear [m]
            "x_tailstrike": 23.68,  # Longitudinal position of critical tailstrike point [m]
            "z_tailstrike": -0.84,  # Vertical position of critical tailstrike point [m]
            "c_tank_c_w": 0.4,  # Fraction of the wing chord occupied by the fuel tank
            "x_tank_c_w": 0.2,  # Fraction of the wing chord where fuel tank starts
            "clmax_w": 1.8,  # Maximum lift coefficient of wing airfoil
            "flap_type": "double slotted",  # Flap type
            "c_flap_c_wing": 0.30,  # Fraction of the wing chord occupied by flaps
            "b_flap_b_wing": 0.60,  # Fraction of the wing span occupied by flaps (including fuselage portion)
            "slat_type": None,  # Slat type
            "c_slat_c_wing": 0.00,  # Fraction of the wing chord occupied by slats
            "b_slat_b_wing": 0.00,  # Fraction of the wing span occupied by slats
            "c_ail_c_wing": 0.27,  # Fraction of the wing chord occupied by aileron
            "b_ail_b_wing": 0.34,  # Fraction of the wing span occupied by aileron
            "h_ground": 35.0
            * ft2m,  # Distance to the ground for ground effect computation [m]
            "k_exc_drag": 0.03,  # Excrescence drag factor
            "altitude_takeoff": 0.0,  # Altitude for takeoff computation [m]
            "distance_takeoff": 1800.0,  # Required takeoff distance [m]
            "altitude_landing": 0.0,  # Altitude for landing computation [m]
            "distance_landing": 1800.0,  # Required landing distance [m] (The actual Fokker100 distance is 1350 m but it is very restrictive compared to the historical regression. Therefore I kept the same TO distance since the aircraft should takeoff and land at the same runway)
            "MLW_frac": 38300 / 41500,  # Max Landing Weight / Max Takeoff Weight
            "altitude_cruise": 35000 * ft2m,  # Cruise altitude [m]
            "Mach_cruise": 0.73,  # Cruise Mach number
            "range_cruise": 1200 * nm2m,  # Cruise range [m]
            "loiter_time": 45 * 60,  # Loiter time [s]
            "altitude_altcruise": 4572,  # Alternative cruise altitude [m]
            "Mach_altcruise": 0.4,  # Alternative cruise Mach number
            "range_altcruise": 200 * nm2m,  # Alternative cruise range [m]
            "W_payload": 107 * 91 * gravity,  # Payload weight [N]
            "xcg_payload": 14.4,  # Longitudinal position of the Payload center of gravity [m]
            "W_crew": 5 * 91 * gravity,  # Crew weight [N]
            "xcg_crew": 2.5,  # Longitudinal position of the Crew center of gravity [m]
            "rho_f": 804,  # Fuel density kg/m3 (This is Jet A-1)
            #'W0_guess' : 40000*gravity # Guess for MTOW
        }

    elif name == "e145xr":

        airplane = {
            "type": "transport",  # Can be 'transport', 'fighter', or 'general'
            "S_w": 51.92,  # Wing area [m2]
            "AR_w": 7.7,  # Wing aspect ratio
            "taper_w": 0.24,  # Wing taper ratio
            "sweep_w": 22.6 * np.pi / 180,  # Wing sweep [rad]
            "dihedral_w": 5 * np.pi / 180,  # Wing dihedral [rad]
            "xr_w": 12.73,  # Longitudinal position of the wing (with respect to the fuselage nose) [m]
            "zr_w": 0.22
            - 1.14,  # Vertical position of the wing (with respect to the fuselage nose) [m]
            "tcr_w": 0.123,  # t/c of the root section of the wing
            "tct_w": 0.096,  # t/c of the tip section of the wing
            "Cht": 0.948,  # Horizontal tail volume coefficient
            "Lc_h": 4.196,  # Non-dimensional lever of the horizontal tail (lever/wing_mac)
            "AR_h": 4.63,  # HT aspect ratio
            "taper_h": 0.54,  # HT taper ratio
            "sweep_h": 22.1 * np.pi / 180,  # HT sweep [rad]
            "dihedral_h": 2 * np.pi / 180,  # HT dihedral [rad]
            "zr_h": 5.28 - 1.14,  # Vertical position of the HT [m]
            "tcr_h": 0.1,  # t/c of the root section of the HT
            "tct_h": 0.1,  # t/c of the tip section of the HT
            "eta_h": 1.0,  # Dynamic pressure factor of the HT
            "Cvt": 0.095,  # Vertical tail volume coefficient
            "Lb_v": 0.556,  # Non-dimensional lever of the vertical tail (lever/wing_span)
            "AR_v": 1.23,  # VT aspect ratio
            "taper_v": 0.745,  # VT taper ratio
            "sweep_v": 29.6 * np.pi / 180,  # VT sweep [rad]
            "zr_v": 2.06 - 1.14,  # Vertical position of the VT [m]
            "tcr_v": 0.1,  # t/c of the root section of the VT
            "tct_v": 0.1,  # t/c of the tip section of the VT
            "L_f": 28.0,  # Fuselage length [m]
            "D_f": 2.28,  # Fuselage diameter [m]
            "x_n": 20.21,  # Longitudinal position of the nacelle frontal face [m]
            "y_n": 1.95,  # Lateral position of the nacelle centerline [m]
            "z_n": 2.01 - 1.14,  # Vertical position of the nacelle centerline [m]
            "L_n": 4.30,  # Nacelle length [m]
            "D_n": 1.51,  # Nacelle diameter [m]
            "n_engines": 2,  # Number of engines
            "n_engines_under_wing": 0,  # Number of engines installed under the wing
            "engine": {
                "model": "Howe turbofan",
                "BPR": 5,  # Engine bypass ratio
                "Cbase": 0.36
                / 3600,  # Base engine TSFC [1/s] (use 'None' for Howe's values)
                # Got value from https://en.wikipedia.org/wiki/Rolls-Royce_AE_3007#cite_note-Roux-10
            },
            "x_nlg": 2.23,  # Longitudinal position of the nose landing gear [m]
            "x_mlg": 16.47,  # Longitudinal position of the main landing gear [m]
            "y_mlg": 2.07,  # Lateral position of the main landing gear [m]
            "z_lg": -1.06 - 1.14,  # Vertical position of the landing gear [m]
            "x_tailstrike": 23.68,  # Longitudinal position of critical tailstrike point [m]
            "z_tailstrike": -0.84,  # Vertical position of critical tailstrike point [m]
            "c_tank_c_w": 0.4,  # Fraction of the wing chord occupied by the fuel tank
            "x_tank_c_w": 0.2,  # Fraction of the wing chord where fuel tank starts
            "clmax_w": 1.8,  # Maximum lift coefficient of wing airfoil
            "flap_type": "double slotted",  # Flap type
            "c_flap_c_wing": 0.2,  # Fraction of the wing chord occupied by flaps
            "b_flap_b_wing": 0.6,  # Fraction of the wing span occupied by flaps (including fuselage portion)
            "slat_type": "slat",  # Slat type
            "c_slat_c_wing": 0.10,  # Fraction of the wing chord occupied by slats
            "b_slat_b_wing": 0.75,  # Fraction of the wing span occupied by slats
            "c_ail_c_wing": 0.27,  # Fraction of the wing chord occupied by ailerons
            "b_ail_b_wing": 0.34,  # Fraction of the wing span occupied by ailerons
            "h_ground": 35.0
            * ft2m,  # Distance to the ground for ground effect computation [m]
            "k_exc_drag": 0.03,  # Excrescence drag factor
            "altitude_takeoff": 0.0,  # Altitude for takeoff computation [m]
            "distance_takeoff": 1800.0,  # Required takeoff distance [m]
            "altitude_landing": 0.0,  # Altitude for landing computation [m]
            "distance_landing": 1800.0,  # Required landing distance [m]
            "MLW_frac": 0.84,  # Max Landing Weight / Max Takeoff Weight
            "altitude_cruise": 37000 * ft2m,  # Cruise altitude [m]
            "Mach_cruise": 0.8,  # Cruise Mach number
            "range_cruise": 2000 * nm2m,  # Cruise range [m]
            "loiter_time": 45 * 60,  # Loiter time [s]
            "altitude_altcruise": 4572,  # Alternative cruise altitude [m]
            "Mach_altcruise": 0.4,  # Alternative cruise Mach number
            "range_altcruise": 200 * nm2m,  # Alternative cruise range [m]
            "W_payload": 50 * 91 * gravity,  # Payload weight [N]
            "xcg_payload": 14.4,  # Longitudinal position of the Payload center of gravity [m]
            "W_crew": 3 * 91 * gravity,  # Crew weight [N]
            "xcg_crew": 2.5,  # Longitudinal position of the Crew center of gravity [m]
            "rho_f": 804,  # Fuel density kg/m3 (This is Jet A-1)
            #'W0_guess' : 24100*gravity # Guess for MTOW
        }

    elif name == "AviaoDoXerife":

        airplane = {
            "type": "transport",  # Can be 'transport', 'fighter', or 'general'
            "S_w": 78.94,  # Wing area [m2]
            "AR_w": 8.2,  # Wing aspect ratio
            "taper_w": 0.272,  # Wing taper ratio
            "sweep_w": 0.380133,  # Wing sweep [rad]
            "dihedral_w": 0.060272,  # Wing dihedral [rad]
            "xr_w": 13.5,  # Longitudinal position of the wing (with respect to the fuselage nose) [m]
            "zr_w": -1.5,  # Vertical position of the wing (with respect to the fuselage nose) [m]
            "tcr_w": 0.123,  # t/c of the root section of the wing
            "tct_w": 0.096,  # t/c of the tip section of the wing
            "Cht": 1.017,  # Horizontal tail volume coefficient
            "Lc_h": 4.027,  # Non-dimensional lever of the horizontal tail (lever/wing_mac)
            "AR_h": 4.517,  # HT aspect ratio
            "taper_h": 0.445,  # HT taper ratio
            "sweep_h": 0.477871,  # HT sweep [rad]
            "dihedral_h": 0.043517,  # HT dihedral [rad]
            "zr_h": 4.359,  # Vertical position of the HT [m]
            "tcr_h": 0.1,  # t/c of the root section of the HT
            "tct_h": 0.1,  # t/c of the tip section of the HT
            "eta_h": 1.0,  # Dynamic pressure factor of the HT
            "Cvt": 0.077,  # Vertical tail volume coefficient
            "Lb_v": 0.45,  # Non-dimensional lever of the vertical tail (lever/wing_span)
            "AR_v": 1.474,  # VT aspect ratio
            "taper_v": 0.553,  # VT taper ratio
            "sweep_v": 0.654266,  # VT sweep [rad]
            "zr_v": 0.0,  # Vertical position of the VT [m]
            "tcr_v": 0.1,  # t/c of the root section of the VT
            "tct_v": 0.1,  # t/c of the tip section of the VT
            "L_f": 25.45,  # Fuselage length [m]
            "D_f": 3.13,  # Fuselage diameter [m]
            "x_n": 23.2,  # Longitudinal position of the nacelle frontal face [m]
            "y_n": 2.6,  # Lateral position of the nacelle centerline [m]
            "z_n": 0.0,  # Vertical position of the nacelle centerline [m]
            "L_n": 1.70,  # Nacelle length [m]
            "D_n": 1.6,  # Nacelle diameter [m]
            "n_engines": 2,  # Number of engines
            "n_engines_under_wing": 0,  # Number of engines installed under the wing
            "engine": {
                "model": "Howe turbofan",  # Check engineTSFC function for options
                "BPR": 3.04,  # Engine bypass ratio
                "Cbase": 0.57 / 3600,
            },
            "x_nlg": 3.6,  # Longitudinal position of the nose landing gear [m]
            "x_mlg": 17.8,  # Longitudinal position of the main landing gear [m]
            "y_mlg": 2.47,  # Lateral position of the main landing gear [m]
            "z_lg": -2.0,  # Vertical position of the landing gear [m]
            "x_tailstrike": 23.68,  # Longitudinal position of critical tailstrike point [m]
            "z_tailstrike": -0.84,  # Vertical position of critical tailstrike point [m]
            "c_tank_c_w": 0.4,  # Fraction of the wing chord occupied by the fuel tank
            "x_tank_c_w": 0.2,  # Fraction of the wing chord where fuel tank starts
            "clmax_w": 1.5,  # Maximum lift coefficient of wing airfoil
            "flap_type": "double slotted",  # Flap type
            "c_flap_c_wing": 0.30,  # Fraction of the wing chord occupied by flaps
            "b_flap_b_wing": 0.60,  # Fraction of the wing span occupied by flaps (including fuselage portion)
            "slat_type": None,  # Slat type
            "c_slat_c_wing": 0.00,  # Fraction of the wing chord occupied by slats
            "b_slat_b_wing": 0.00,  # Fraction of the wing span occupied by slats
            "c_ail_c_wing": 0.27,  # Fraction of the wing chord occupied by aileron
            "b_ail_b_wing": 0.34,  # Fraction of the wing span occupied by aileron
            "h_ground": 35.0
            * ft2m,  # Distance to the ground for ground effect computation [m]
            "k_exc_drag": 0.03,  # Excrescence drag factor
            "altitude_takeoff": 0.0,  # Altitude for takeoff computation [m]
            "distance_takeoff": 1600,  # Required takeoff distance [m]
            "altitude_landing": 0.0,  # Altitude for landing computation [m]
            "distance_landing": 1600.0,  # Required landing distance [m] (The actual Fokker100 distance is 1350 m but it is very restrictive compared to the historical regression. Therefore I kept the same TO distance since the aircraft should takeoff and land at the same runway)
            "MLW_frac": 38300
            / 41500,  # Max Landing Weight / Max Takeoff Weight NAO SEI
            "altitude_cruise": 10000,  # Cruise altitude [m]
            "Mach_cruise": 0.79,  # Cruise Mach number
            "range_cruise": 3000 * 10**3,  # Cruise range [m]
            "loiter_time": 45 * 60,  # Loiter time [s]
            "altitude_altcruise": 4572,  # Alternative cruise altitude [m]
            "Mach_altcruise": 0.4,  # Alternative cruise Mach number
            "range_altcruise": 200 * nm2m,  # Alternative cruise range [m]
            "W_payload": 76 * 91 * gravity,  # Payload weight [N]
            "xcg_payload": 14.4,  # Longitudinal position of the Payload center of gravity [m]
            "W_crew": 4 * 91 * gravity,  # Crew weight [N]
            "xcg_crew": 2.5,  # Longitudinal position of the Crew center of gravity [m]
            "rho_f": 804,  # Fuel density kg/m3 (This is Jet A-1)
            #'W0_guess' : 40000*gravity # Guess for MTOW
        }

    elif name == "F70":
        airplane = {
            "type": "transport",  # Can be 'transport', 'fighter', or 'general'
            "S_w": 93.5,  # Wing area [m2]
            "AR_w": 8.4,  # Wing aspect ratio
            "taper_w": 0.28,  # Wing taper ratio
            "sweep_w": 19.9 * np.pi / 180,  # Wing sweep [rad]
            "dihedral_w": 3 * np.pi / 180,  # Wing dihedral [rad]
            "xr_w": 11.72,  # Longitudinal position of the wing (with respect to the fuselage nose) [m]
            "zr_w": -1.01,  # Vertical position of the wing (with respect to the fuselage nose) [m]
            "tcr_w": 0.123,  # t/c of the root section of the wing
            "tct_w": 0.096,  # t/c of the tip section of the wing
            "Cht": 0.85,  # Horizontal tail volume coefficient
            "Lc_h": 3.5,  # Non-dimensional lever of the horizontal tail (lever/wing_mac)
            "AR_h": 4.41,  # HT aspect ratio
            "taper_h": 0.42,  # HT taper ratio
            "sweep_h": 27.14 * np.pi / 180,  # HT sweep [rad]
            "dihedral_h": 4.2 * np.pi / 180,  # HT dihedral [rad]
            "zr_h": 4.79,  # Vertical position of the HT [m]
            "tcr_h": 0.1,  # t/c of the root section of the HT
            "tct_h": 0.1,  # t/c of the tip section of the HT
            "eta_h": 1.0,  # Dynamic pressure factor of the HT
            "Cvt": 0.06,  # Vertical tail volume coefficient
            "Lb_v": 0.39,  # Non-dimensional lever of the vertical tail (lever/wing_span)
            "AR_v": 1.184,  # VT aspect ratio
            "taper_v": 0.704,  # VT taper ratio
            "sweep_v": 38.57 * np.pi / 180,  # VT sweep [rad]
            "zr_v": 1.65,  # Vertical position of the VT [m]
            "tcr_v": 0.1,  # t/c of the root section of the VT
            "tct_v": 0.1,  # t/c of the tip section of the VT
            "L_f": 27.878,  # Fuselage length [m]
            "D_f": 3.3,  # Fuselage diameter [m]
            "x_n": 17.52,  # Longitudinal position of the nacelle frontal face [m]
            "y_n": 2.7,  # Lateral position of the nacelle centerline [m]
            "z_n": 0.45,  # Vertical position of the nacelle centerline [m]
            "L_n": 5.02,  # Nacelle length [m]
            "D_n": 1.64,  # Nacelle diameter [m]
            "n_engines": 2,  # Number of engines
            "n_engines_under_wing": 0,  # Number of engines installed under the wing
            "engine": {
                "model": "Howe turbofan",  # Check engineTSFC function for options
                "BPR": 3.04,  # Engine bypass ratio
                "Cbase": 0.57 / 3600,
            },
            "x_nlg": 3.64,  # Longitudinal position of the nose landing gear [m]
            "x_mlg": 15.18,  # Longitudinal position of the main landing gear [m]
            "y_mlg": 2.52,  # Lateral position of the main landing gear [m]
            "z_lg": -3.05,  # Vertical position of the landing gear [m]
            "x_tailstrike": 23.68,  # Longitudinal position of critical tailstrike point [m]
            "z_tailstrike": -0.84,  # Vertical position of critical tailstrike point [m]
            "c_tank_c_w": 0.4,  # Fraction of the wing chord occupied by the fuel tank
            "x_tank_c_w": 0.2,  # Fraction of the wing chord where fuel tank starts
            "clmax_w": 1.5,  # Maximum lift coefficient of wing airfoil
            "flap_type": "double slotted",  # Flap type
            "c_flap_c_wing": 0.30,  # Fraction of the wing chord occupied by flaps
            "b_flap_b_wing": 0.60,  # Fraction of the wing span occupied by flaps (including fuselage portion)
            "slat_type": None,  # Slat type
            "c_slat_c_wing": 0.00,  # Fraction of the wing chord occupied by slats
            "b_slat_b_wing": 0.00,  # Fraction of the wing span occupied by slats
            "c_ail_c_wing": 0.27,  # Fraction of the wing chord occupied by aileron
            "b_ail_b_wing": 0.34,  # Fraction of the wing span occupied by aileron
            "h_ground": 35.0
            * ft2m,  # Distance to the ground for ground effect computation [m]
            "k_exc_drag": 0.03,  # Excrescence drag factor
            "altitude_takeoff": 0.0,  # Altitude for takeoff computation [m]
            "distance_takeoff": 1300,  # Required takeoff distance [m]
            "altitude_landing": 0.0,  # Altitude for landing computation [m]
            "distance_landing": 1210.0,  # Required landing distance [m] (The actual Fokker100 distance is 1350 m but it is very restrictive compared to the historical regression. Therefore I kept the same TO distance since the aircraft should takeoff and land at the same runway)
            "MLW_frac": 38300
            / 41500,  # Max Landing Weight / Max Takeoff Weight NAO SEI
            "altitude_cruise": 10000,  # Cruise altitude [m]
            "Mach_cruise": 0.77,  # Cruise Mach number
            "range_cruise": 2009 * 10**3,  # Cruise range [m]
            "loiter_time": 45 * 60,  # Loiter time [s]
            "altitude_altcruise": 4572,  # Alternative cruise altitude [m]
            "Mach_altcruise": 0.4,  # Alternative cruise Mach number
            "range_altcruise": 200 * nm2m,  # Alternative cruise range [m]
            "W_payload": 10890 * gravity,  # Payload weight [N]
            "xcg_payload": 14.4,  # Longitudinal position of the Payload center of gravity [m]
            "W_crew": 4 * 91 * gravity,  # Crew weight [N]
            "xcg_crew": 2.5,  # Longitudinal position of the Crew center of gravity [m]
            "rho_f": 804,  # Fuel density kg/m3 (This is Jet A-1)
            "W0_guess": 36740 * gravity,  # Guess for MTOW
        }

        ##### DADOS QUE NÃO CONSEGUIMOS INSERIR:
        ## Empty Weight (Kg)	22784
        ## Service ceiling (m)	11000
        ## Baggage compartment volume (m³)	12,78
        ## Wing location (high,low,...)	low

        ## HT area (m²)	22,84
        ## VT area (m²)	14,41

        ## Cockpit length / Diameter	1,1545
        ## Tail cone length / Diameter	0,5333

        ## Number of galleys     2
        ## Number of lavatories    2
        ## Number of seats abreast 3
        ## Number of doors         4
        ## Number of window exits  2

        ## Total maximum thrust or power (kN)	123,2
        ## Engine TSFC or ESFC (lbm/h/lbf)	0,69

        ## Wing loading (Kg/m²)	392,9412
        ## Thrust-to-weight ratio	0,3418
        ##########################################

    elif name == "F70_XerifeEdition":
        # COMENTÁRIOS
        
        # PROCESSO DE OTIMIZAÇÃO/AJUSTE
            # Como a planilha "PrestoCabin" utiliza os valores estimados para o alcance e o MTOW
            # da aeronave, é necessário fazer um processo iterativo com o código/análises
            
        # COMENTÁRRIOS SOBRE COMO VALORES DE PARÂMETROS SÃO AJUSTADOS
            # Na planilha "PrestoCabin", os valores atuais de MTOW e de alcance são de 36000 kg e 3900 km
            # xcg_payload <- adaptação do valor do Fokker 100 para o nosso comprimento de fuselagem
            # L_f (Fuselage length [m]) <- Valor retirado da planilha PrestoCabin (ln 372)
            # D_f (Fuselage diameter [m]) <- Valor retirado da planilha PrestoCabin (ln 136)
            
        # ITENS QUE PRECISAM SER CORRIGIDOS OU APRIMORADOS
            # O comprimento da fuselagem diminuiu em comparação com o F70 original (2,53 m a menos). 
            # Por conta disso, algumas medidas e posicionamentos longitudinais precisam ser modificados
            # (posição da asa, posição das empenagens, posição do motor, posição do trem de pouso 
            # principal, comprimento da nacelle do motor)
            
        
        airplane = {
            "type": "transport",  # Can be 'transport', 'fighter', or 'general'
            "S_w": 92,  # Wing area [m2]
            "AR_w": 8.0,  # Wing aspect ratio
            "taper_w": 0.28,  # Wing taper ratio
            "sweep_w": 16.5 * np.pi / 180,  # Wing sweep [rad]
            "dihedral_w": 3.3 * np.pi / 180,  # Wing dihedral [rad]
            "xr_w": 12.82,  # Longitudinal position of the wing (with respect to the fuselage nose) [m]
            "zr_w": -1.01,  # Vertical position of the wing (with respect to the fuselage nose) [m]
            "tcr_w": 0.140,  # t/c of the root section of the wing
            "tct_w": 0.075,  # t/c of the tip section of the wing
            "Cht": 0.78,  # Horizontal tail volume coefficient
            "Lc_h": 3.7,  # Non-dimensional lever of the horizontal tail (lever/wing_mac)
            "AR_h": 4.75,  # HT aspect ratio
            "taper_h": 0.36,  # HT taper ratio
            "sweep_h": 29.14 * np.pi / 180,  # HT sweep [rad]
            "dihedral_h": 4.2 * np.pi / 180,  # HT dihedral [rad]
            "zr_h": 5.2,  # Vertical position of the HT [m]
            "tcr_h": 0.1,  # t/c of the root section of the HT
            "tct_h": 0.1,  # t/c of the tip section of the HT
            "eta_h": 1.0,  # Dynamic pressure factor of the HT
            "Cvt": 0.05,  # Vertical tail volume coefficient
            "Lb_v": 0.42,  # Non-dimensional lever of the vertical tail (lever/wing_span)
            "AR_v": 1.314,  # VT aspect ratio
            "taper_v": 0.754,  # VT taper ratio
            "sweep_v": 35.57 * np.pi / 180,  # VT sweep [rad]
            "zr_v": 3.24/2,  # Vertical position of the VT [m]
            "tcr_v": 0.1,  # t/c of the root section of the VT
            "tct_v": 0.1,  # t/c of the tip section of the VT
            "L_f": 29.27,  # Fuselage length [m] <- Otimizado de acordo com a planilha PrestoCabin (ln 372)
            "D_f": 3.24,  # Fuselage diameter [m] <- Otimizado de acordo com a planilha PrestoCabin (ln 136)
            "x_n": 22.5,  # Longitudinal position of the nacelle frontal face [m]
            "y_n": 2.80,  # Lateral position of the nacelle centerline [m]
            "z_n": 0.85,  # Vertical position of the nacelle centerline [m]
            "L_n": 3.27,  # Nacelle length [m]
            "D_n": 1.64,  # Nacelle diameter [m]
            "n_engines": 2,  # Number of engines
            "n_engines_under_wing": 0,  # Number of engines installed under the wing
            "engine": {
                "model": "Howe turbofan",  # Check engineTSFC function for options
                "BPR": 4.9,  # Engine bypass ratio
                "Cbase": 0.39/ 3600,
                "weight": 1120*gravity
            }, # Motor CF34-C5
            "x_nlg": 3.64,  # Longitudinal position of the nose landing gear [m]
            "x_mlg": 16.82,  # Longitudinal position of the main landing gear [m]
            "y_mlg": 2.52,  # Lateral position of the main landing gear [m]
            "z_lg": -3.05,  # Vertical position of the landing gear [m]
            "x_tailstrike": 22,  # Longitudinal position of critical tailstrike point [m]
            "z_tailstrike": -0.64,  # Vertical position of critical tailstrike point [m]
            "c_tank_c_w": 0.4,  # Fraction of the wing chord occupied by the fuel tank
            "x_tank_c_w": 0.2,  # Fraction of the wing chord where fuel tank starts
            "clmax_w": 1.5,  # Maximum lift coefficient of wing airfoil
            "flap_type": "double slotted",  # Flap type
            "c_flap_c_wing": 0.30,  # Fraction of the wing chord occupied by flaps
            "b_flap_b_wing": 0.60,  # Fraction of the wing span occupied by flaps (including fuselage portion)
            "slat_type": None,  # Slat type
            "c_slat_c_wing": 0.00,  # Fraction of the wing chord occupied by slats
            "b_slat_b_wing": 0.00,  # Fraction of the wing span occupied by slats
            "c_ail_c_wing": 0.27,  # Fraction of the wing chord occupied by aileron
            "b_ail_b_wing": 0.34,  # Fraction of the wing span occupied by aileron
            "h_ground": 35.0
            * ft2m,  # Distance to the ground for ground effect computation [m]
            "k_exc_drag": 0.03,  # Excrescence drag factor
            "altitude_takeoff": 0.0,  # Altitude for takeoff computation [m]
            "distance_takeoff": 1600,  # Required takeoff distance [m]
            "altitude_landing": 0.0,  # Altitude for landing computation [m]
            "distance_landing": 1600, # Required landing distance [m] (The actual Fokker100 distance is 1350 m but it is very restrictive compared to the historical regression. Therefore I kept the same TO distance since the aircraft should takeoff and land at the same runway)
            "MLW_frac": 34400/37527,  # Max Landing Weight / Max Takeoff Weight NAO SEI
            "altitude_cruise": 10000,  # Cruise altitude [m]
            "Mach_cruise": 0.76,  # Cruise Mach number
            "range_cruise": 4550 * 10**3,  # Cruise range [m]
            "loiter_time": 45 * 60,  # Loiter time [s]
            "altitude_altcruise": 4572,  # Alternative cruise altitude [m]
            "Mach_altcruise": 0.4,  # Alternative cruise Mach number
            "range_altcruise": 200 * nm2m,  # Alternative cruise range [m]
            "W_payload": 76 * 91 * gravity,  # Payload weight [N]
            "xcg_payload": 12.26,  # Longitudinal position of the Payload center of gravity [m]
            "W_crew": 4 * 91 * gravity,  # Crew weight [N]
            "xcg_crew": 2.5,  # Longitudinal position of the Crew center of gravity [m]
            "rho_f": 804,  # Fuel density kg/m3 (This is Jet A-1)
            #"W0_guess": 36740 * gravity,  # Guess for MTOW
        }

        

    return airplane
