# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 18:30:02 2022

@author: VHSARME
"""

# ANÁLISE DE SENSIBILIDADE

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import designTool as dt
from pymoo.core.problem import Problem
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.sampling.lhs import LHS
from aux_tools_doe import corrdot

# =========================================

# SETUP


# Define design function
def des_funcs(Xinp):

    airplane = dt.standard_airplane("F70_XerifeEdition")
    Range, AR_w, S_w, sweep_w, Mach_cruise, alt_cruise, xr_w, Cht, Lc_h, Cvt, Lb_v = Xinp
    airplane["range_cruise"] = Range
    airplane["AR_w"] = AR_w
    airplane["S_w"] = S_w
    airplane["sweep_w"] = sweep_w
    airplane["Mach_cruise"] = Mach_cruise
    airplane["altitude_cruise"] = alt_cruise
    airplane["xr_w"] = xr_w
    airplane["x_mlg"] = xr_w + 4.28
    airplane["Cht"] = Cht
    airplane["Lc_h"] = Lc_h
    airplane["Cvt"] = Cvt
    airplane["Lb_v"] = Lb_v
    
    result = dt.analyze(
        airplane=airplane,
        print_log=False,  # Plot results on the terminal screen
        plot=False,  # Generate 3D plot of the aircraft
    )

    Wf_W0 = result["Wf"] / result["W0"]
    We_W0 = result["We"] / result["W0"]
    T0_W0 = result["T0"] / result["W0"]
    W0_S = result["W0"] / result["S_w"]
    SM_fwd = result["SM_fwd"]
    SM_aft = result["SM_aft"]

    # Returns
    return Wf_W0, We_W0, T0_W0, W0_S, SM_fwd, SM_aft


# Give number of input variables
n_var = 11

# Lower and upeer bounds of each input variable

Sweep_Lower = (16.5 * np.pi / 180) * 0.4
Sweep_Upper = 35 * np.pi / 180
MachCruise_Lower = 0.75
MachCruise_Upper = 0.80
AR_Lower = 6
AR_Upper = 11
SW_Lower = 60
SW_Upper = 130
Range_Lower = 4500000
Range_Upper = 5500000
AltCruise_Lower = 10000
AltCruise_Upper = 12000
xr_w_Lower = 10
xr_w_Upper = 15
Cht_Lower = 0.75
Cht_Upper = 1.1
Lch_Lower = 3.7 * 0.9
Lch_Upper = 3.7 * 1.1
Cvt_Lower = 0.04
Cvt_Upper = 0.12
Lbv_Lower = 0.42 * 0.9
Lbv_Upper = 0.42 * 1.1

lb = [Range_Lower, AR_Lower, SW_Lower, Sweep_Lower, MachCruise_Lower, AltCruise_Lower, xr_w_Lower,
      Cht_Lower, Lch_Lower, Cvt_Lower, Lbv_Lower]
ub = [Range_Upper, AR_Upper, SW_Upper, Sweep_Upper, MachCruise_Upper, AltCruise_Upper, xr_w_Upper,
      Cht_Upper, Lch_Upper, Cvt_Upper, Lbv_Upper]

# Desired number of samples
n_samples = 100

# Sampling type
# sampler = FloatRandomSampling()
sampler = LHS()

# Plot type (0-simple, 1-complete)
plot_type = 1
# =========================================

# EXECUTION

# Set random seed to make results repeatable
np.random.seed(123)

# Initialize problem with lower and upper
problem = Problem(n_var=n_var, xl=lb, xu=ub)

# Draw samples
X = sampler(problem, n_samples).get("X")

# Samples are originally between 0 and 1,
# so we need to scale them to the desired interval
# for ii in range(n_inputs):
#    X[:,ii] = lb[ii] + (ub[ii] - lb[ii])*X[:,ii]

# Execute all cases and store outputs
y1_samples = np.zeros(n_samples)
y2_samples = np.zeros(n_samples)
y3_samples = np.zeros(n_samples)
y4_samples = np.zeros(n_samples)
y5_samples = np.zeros(n_samples)
y6_samples = np.zeros(n_samples)

for ii in range(n_samples):

    # Evaluate sample
    (y1, y2, y3, y4, y5, y6) = des_funcs(X[ii, :])

    # Store the relevant information
    y1_samples[ii] = y1
    y2_samples[ii] = y2
    y3_samples[ii] = y3
    y4_samples[ii] = y4
    y5_samples[ii] = y5
    y6_samples[ii] = y6

# Create a pandas dataframe with all the information
df_w = pd.DataFrame(
    {   "AR_w": X[:, 1],
        "S_w": X[:, 2],
        "sweep_w": X[:, 3],
        "xr_w": X[:,6],
        "Wf_W0": y1_samples,
        "We_W0": y2_samples,
        "T0 / W0": y3_samples,
        "W0 / Sw": y4_samples,
        "SM_fwd": y5_samples,
        "SM_aft": y6_samples,
    }
)

# Plot the correlation matrix
sns.set(style="white", font_scale=0.8)

if plot_type == 0:

    # Simple plot
    fig = sns.pairplot(df_w, corner=True)

elif plot_type == 1:

    # Complete plot
    # based on: https://stackoverflow.com/questions/48139899/correlation-matrix-plot-with-coefficients-on-one-side-scatterplots-on-another
    fig = sns.PairGrid(df_w, diag_sharey=False)
    fig.map_lower(sns.regplot, lowess=True, line_kws={"color": "black"})
    fig.map_diag(sns.histplot)
    fig.map_upper(corrdot)

# Plot window
plt.tight_layout()
plt.show()

# Create a pandas dataframe with all the information
df_op = pd.DataFrame(
    {   "alt_cruise": X[:, 5],
        "mach": X[:, 4],
        "Wf_W0": y1_samples,
        "We_W0": y2_samples,
        "T0 / W0": y3_samples,
        "W0 / Sw": y4_samples,
        "SM_fwd": y5_samples,
        "SM_aft": y6_samples,
    }
)

# Plot the correlation matrix
sns.set(style="white", font_scale=0.8)

if plot_type == 0:

    # Simple plot
    fig = sns.pairplot(df_op, corner=True)

elif plot_type == 1:

    # Complete plot
    # based on: https://stackoverflow.com/questions/48139899/correlation-matrix-plot-with-coefficients-on-one-side-scatterplots-on-another
    fig = sns.PairGrid(df_op, diag_sharey=False)
    fig.map_lower(sns.regplot, lowess=True, line_kws={"color": "black"})
    fig.map_diag(sns.histplot)
    fig.map_upper(corrdot)

# Plot window
plt.tight_layout()
plt.show()