'''
INSTITUTO TECNOLÓGICO DE AERONÁUTICA
DIVISÃO DE ENGENHARIA AERONÁUTICA

This script generates correlation plots

REQUIRED PACKAGES (if not using Anaconda)
pip3 install pandas
pip3 install seaborn
pip3 install statsmodels
pip3 install pymoo==0.6.0

Cap. Ney Sêcco 2022
'''

#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pymoo.core.problem import Problem
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.sampling.lhs import LHS
from aux_tools_doe import corrdot

#=========================================

# SETUP

# Define design function
def des_funcs(Xinp):

    x1, x2, x3 = Xinp

    y1 = x1**2 + (x2-2)**2 - 0.0001*x3

    y2 = x1 - x2

    # Returns
    return y1, y2

# Give number of input variables
n_var = 3

# Lower and upeer bounds of each input variable
lb = [-5.0, -5.0, -5.0]
ub = [ 5.0,  5.0,  5.0]

# Desired number of samples
n_samples = 200

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
for ii in range(n_samples):

    # Evaluate sample
    (y1,y2) = des_funcs(X[ii,:])

    # Store the relevant information
    y1_samples[ii] = y1
    y2_samples[ii] = y2

# Create a pandas dataframe with all the information
df = pd.DataFrame({'x1' : X[:,0],
                   'x2' : X[:,1],
                   'x3' : X[:,2],
                   'y1' : y1_samples,
                   'y2' : y2_samples})

# Plot the correlation matrix
sns.set(style='white', font_scale=1.4)

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

fig.savefig('doe.png')
