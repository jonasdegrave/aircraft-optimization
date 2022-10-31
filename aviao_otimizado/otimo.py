# -*- coding: utf-8 -*-
"""
Created on Wed Oct 26 18:30:02 2022

@author: VHSARME
"""

# ANÁLISE DE SENSIBILIDADE

#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import designTool as dt
from pymoo.core.problem import Problem
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.sampling.lhs import LHS
from aux_tools_doe import corrdot

#=========================================

# SETUP


# Define design function
def des_funcs(Xinp):
    
    airplane = dt.standard_airplane("F70_XerifeEdition")
    Range = Xinp
    airplane['range_cruise'] = Range
    #airplane['AR_w'] = AR_w
    #airplane['S_w'] = S_w
    result = dt.analyze(
        airplane=airplane,
        print_log=False,  # Plot results on the terminal screen
        plot=False,  # Generate 3D plot of the aircraft
    )

    W0 = result['W0']/9.81
    W_fuel = result['Wf']/9.81
    #W_empty = result['We']/9.81
    deltaS_wlan = result['deltaS_wlan']
    T0_W0 = result['T0']/result['W0']
    W0_S = result['W0']/result['S_w']
    
    

    # Returns
    return W0, W_fuel, deltaS_wlan, T0_W0, W0_S

# Give number of input variables
n_var = 1

# Lower and upeer bounds of each input variable
lb = [3000]
ub = [6000]

# Desired number of samples
n_samples = 100

# Sampling type
#sampler = FloatRandomSampling()
sampler = LHS()

# Plot type (0-simple, 1-complete)
plot_type = 1
#=========================================

# EXECUTION

# Set random seed to make results repeatable
np.random.seed(123)

# Initialize problem with lower and upper
problem = Problem(n_var=n_var, xl=lb, xu=ub)

# Draw samples
X = sampler(problem, n_samples).get("X")

# Samples are originally between 0 and 1,
# so we need to scale them to the desired interval
#for ii in range(n_inputs):
#    X[:,ii] = lb[ii] + (ub[ii] - lb[ii])*X[:,ii]

# Execute all cases and store outputs
y1_samples = np.zeros(n_samples)
y2_samples = np.zeros(n_samples)
y3_samples = np.zeros(n_samples)
y4_samples = np.zeros(n_samples)
y5_samples = np.zeros(n_samples)
for ii in range(n_samples):

    # Evaluate sample
    (y1,y2,y3,y4,y5) = des_funcs(X[ii,:])

    # Store the relevant information
    y1_samples[ii] = y1
    y2_samples[ii] = y2
    y3_samples[ii] = y3
    y4_samples[ii] = y4
    y5_samples[ii] = y5

# Create a pandas dataframe with all the information
df = pd.DataFrame({'Range' : X[:,0],
                   #'AR_w' : X[:,0],
                   #'S_w' : X[:,1],
                   'W0' : y1_samples,
                   'W_fuel' : y2_samples,
                   'deltaS_wlan' : y3_samples,
                   'T0 / W0' : y4_samples,
                   'W0 / Sw' : y5_samples})

# Plot the correlation matrix
sns.set(style='white', font_scale=1.1)

if plot_type == 0:

    # Simple plot
    fig = sns.pairplot(df,corner=True)

elif plot_type == 1:

    # Complete plot
    # based on: https://stackoverflow.com/questions/48139899/correlation-matrix-plot-with-coefficients-on-one-side-scatterplots-on-another
    fig = sns.PairGrid(df, diag_sharey=False)
    fig.map_lower(sns.regplot, lowess=True, line_kws={'color': 'black'})
    fig.map_diag(sns.histplot)
    fig.map_upper(corrdot)

# Plot window
plt.tight_layout()
plt.show()