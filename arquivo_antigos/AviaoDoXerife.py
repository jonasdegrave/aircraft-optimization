# Sample script on how to use designTool.
# Remember to save this script in the same directory as designTool.py

# IMPORTS
import designTool as dt
import numpy as np
import pprint

# Inputs of the Geometry Module

airplane = {
    # "AR_h": 4.517,
    # "AR_v": 1.474,
    # "AR_w": 8.2,
    # "Cht": 1.017,
    # "Cvt": 0.077,
    # "D_f": 3.13,
    # "D_n": 1.6,
    # "L_f": 25.45,
    # "L_n": 1.70,
    # "Lb_v": 0.45,
    # "Lc_h": 1.017,
    # "MLW_frac": 0.9228915662650602,
    # "Mach_altcruise": 0.4,
    # "Mach_cruise": 0.79,
    # "S_w": 78.94,
    # "W_crew": 4463.55, NÃO SEI O QUE É
    # "W_payload": 95519.97, NÃO SEI O QUE É
    # "altitude_altcruise": 4572, NÃO SEI O QUE É
    # "altitude_cruise": 10058.4,
    # "altitude_landing": 0.0,NÃO SEI O QUE É
    # "altitude_takeoff": 0.0,NÃO SEI O QUE É
    # "b_ail_b_wing": 0.34,NÃO SEI O QUE É
    # "b_flap_b_wing": 0.6,NÃO SEI O QUE É
    # "b_slat_b_wing": 0.0,NÃO SEI O QUE É
    # "c_ail_c_wing": 0.27,NÃO SEI O QUE É
    # "c_flap_c_wing": 0.3,NÃO SEI O QUE É
    # "c_slat_c_wing": 0.0,NÃO SEI O QUE É
    # "c_tank_c_w": 0.4,NÃO SEI O QUE É
    # "clmax_w": 1.5, #ESCOLHEMOS 'NO FLAP'
    # "dihedral_h": 0.043517,
    # "dihedral_w": 0.060272,
    # "distance_landing": 1503.0,
    # "distance_takeoff": 1330.0,
    # "engine": {"BPR": 3.04, "Cbase": 0.00015833333333333332, "model": "Howe turbofan"}, NÃO SEI O QUE É
    # "eta_h": 1.0, NÃO SEI O QUE É
    # "flap_type": "double slotted", NÃO SEI O QUE É
    # "h_ground": 10.668000000000001,NÃO SEI O QUE É
    # "k_exc_drag": 0.03,NÃO SEI O QUE É
    # "loiter_time": 2700,NÃO SEI O QUE É
    # "n_engines": 2,
    # "n_engines_under_wing": 0,
    # "range_altcruise": 370400.0,NÃO SEI O QUE É
    # "range_cruise": 2222400.0,NÃO SEI O QUE É
    # "rho_f": 804,NÃO SEI O QUE É
    # "slat_type": None,NÃO SEI O QUE É
    # "sweep_h": 0.477871,
    # "sweep_v": 0.654266,
    # "sweep_w": 0.380133,
    # "taper_h": 0.445,
    # "taper_v": 1.474,
    # "taper_w": 0.272,
    # "tcr_h": 0.1, NÃO SEI O QUE É
    # "tcr_v": 0.1, NÃO SEI O QUE É
    # "tcr_w": 0.123, NÃO SEI O QUE É
    # "tct_h": 0.1, NÃO SEI O QUE É
    # "tct_v": 0.1, NÃO SEI O QUE É
    # "tct_w": 0.096, NÃO SEI O QUE É
    # "type": "transport",NÃO SEI O QUE É
    # "x_mlg": 17.8, NÃO SEI O QUE É
    # "x_n": 23.2, NÃO SEI O QUE É
    # "x_nlg": 3.6, NÃO SEI O QUE É
    # "x_tailstrike": 23.68, NÃO SEI O QUE É
    # "x_tank_c_w": 0.2, NÃO SEI O QUE É
    # "xcg_crew": 2.5,NÃO SEI O QUE É
    # "xcg_payload": 14.4,NÃO SEI O QUE É
    # "xr_w": 13.5, NÃO SEI O QUE É
    # "y_mlg": 2.47, NÃO SEI O QUE É
    # "y_n": 2.6, NÃO SEI O QUE É
    # "z_lg": -2.0, NÃO SEI O QUE É
    # "z_n": 0.0, NÃO SEI O QUE É
    # "z_tailstrike": -0.84, NÃO SEI O QUE É
    # "zr_h": 4.359, NÃO SEI O QUE É
    # "zr_v": 0.0, NÃO SEI O QUE É
    # "zr_w": -1.5, NÃO SEI O QUE É
}

# Execute the geometry function
dt.geometry(airplane)

# Print updated dictionary
print("airplane = " + pprint.pformat(airplane))

# Generate 3D plot
dt.plot3d(airplane)
