'''
INSTITUTO TECNOLÓGICO DE AERONÁUTICA
DIVISÃO DE ENGENHARIA AERONÁUTICA

This script introduces how to use pymoo to
generate random or LHS samples for DOE

REQUIRED PACKAGES (if not using Anaconda)
pip3 install pymoo==0.6.0

Cap. Ney Sêcco 2022
'''


#IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from pymoo.core.problem import Problem
from pymoo.operators.sampling.rnd import FloatRandomSampling
from pymoo.operators.sampling.lhs import LHS

#=========================================

# EXECUTION

# Parameters
n_var = 2
n_samples = 20
lb = [0]*n_var # Lower bounds
ub = [1]*n_var # Upper bounds

# Set random seed to make results repeatable
np.random.seed(123)

# Initialize problem with lower and upper
problem = Problem(n_var=n_var, xl=lb, xu=ub)

# Initialize samplers
rnd_sampler = FloatRandomSampling()
lhs_sampler = LHS()

# Draw samples
X_rnd = rnd_sampler(problem, n_samples).get("X")
X_lhs = lhs_sampler(problem, n_samples).get("X")

# Plot samples
fig = plt.figure()

plt.subplot(121)
plt.plot(X_rnd[:,0], X_rnd[:,1], 'o')
plt.title('Random', fontsize=12)
plt.xlabel('X1', fontsize=12)
plt.ylabel('X2', fontsize=12)

plt.subplot(122)
plt.plot(X_lhs[:,0], X_lhs[:,1], 'o')
plt.title('LHS', fontsize=12)
plt.xlabel('X1', fontsize=12)
plt.ylabel('X2', fontsize=12)

plt.tight_layout()
plt.show()