# Sample script on how to use designTool.
# Remember to save this script in the same directory as designTool.py

# IMPORTS
import designTool as dt
import numpy as np
import pprint

# Inputs of the Geometry Module

airplane = {
    "AR_h": 4.64,
    "AR_v": 1.27,
    "AR_w": 8.43,
    "Cht": 0.94,
    "Cvt": 0.088,
    "D_f": 3.3,
    "D_n": 1.5,
    "L_f": 32.5,
    "L_n": 4.3,
    "Lb_v": 0.55,
    "Lc_h": 4.83,
    "MLW_frac": 0.9228915662650602,
    "Mach_altcruise": 0.4,
    "Mach_cruise": 0.73,
    "S_w": 93.5,
    "W_crew": 4463.55,
    "W_payload": 95519.97,
    "altitude_altcruise": 4572,
    "altitude_cruise": 10668.0,
    "altitude_landing": 0.0,
    "altitude_takeoff": 0.0,
    "b_ail_b_wing": 0.34,
    "b_flap_b_wing": 0.6,
    "b_slat_b_wing": 0.0,
    "c_ail_c_wing": 0.27,
    "c_flap_c_wing": 0.3,
    "c_slat_c_wing": 0.0,
    "c_tank_c_w": 0.4,
    "clmax_w": 1.8,
    "dihedral_h": 0.03490658503988659,
    "dihedral_w": 0.08726646259971647,
    "distance_landing": 1800.0,
    "distance_takeoff": 1800.0,
    "engine": {"BPR": 3.04, "Cbase": 0.00015833333333333332, "model": "Howe turbofan"},
    "eta_h": 1.0,
    "flap_type": "double slotted",
    "h_ground": 10.668000000000001,
    "k_exc_drag": 0.03,
    "loiter_time": 2700,
    "n_engines": 2,
    "n_engines_under_wing": 0,
    "range_altcruise": 370400.0,
    "range_cruise": 2222400.0,
    "rho_f": 804,
    "slat_type": None,
    "sweep_h": 0.4537856055185257,
    "sweep_v": 0.715584993317675,
    "sweep_w": 0.3045599544730105,
    "taper_h": 0.39,
    "taper_v": 0.74,
    "taper_w": 0.235,
    "tcr_h": 0.1,
    "tcr_v": 0.1,
    "tcr_w": 0.123,
    "tct_h": 0.1,
    "tct_v": 0.1,
    "tct_w": 0.096,
    "type": "transport",
    "x_mlg": 17.8,
    "x_n": 23.2,
    "x_nlg": 3.6,
    "x_tailstrike": 23.68,
    "x_tank_c_w": 0.2,
    "xcg_crew": 2.5,
    "xcg_payload": 14.4,
    "xr_w": 13.5,
    "y_mlg": 2.47,
    "y_n": 2.6,
    "z_lg": -2.0,
    "z_n": 0.0,
    "z_tailstrike": -0.84,
    "zr_h": 4.359,
    "zr_v": 0.0,
    "zr_w": -1.5,
}

# Execute the geometry function
dt.geometry(airplane)

# Print updated dictionary
print("airplane = " + pprint.pformat(airplane))

# Generate 3D plot
dt.plot3d(airplane)
