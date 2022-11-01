"""
PROGRAMA DE ESPECIALIZAÇÃO EM ENGENHARIA AERONÁUTICA
OTIMIZAÇÃO MULTIDISCIPLINAR

Código com aplicações de ferramentas do Python para otimização
restringida

Cap. Ney Sêcco 2020
"""

# IMPORTS
import numpy as np
from scipy.optimize import minimize
import auxmod_mdo as am
import matplotlib.pyplot as plt

# EXECUTION

# Define history list
xlist = []
flist = []
g1list = []
g2list = []

# Define objective function


def objfun(x):
    f = 1 - np.exp(-0.1 * (x[0] ** 2 + x[1] ** 2))
    xlist.append(x)
    flist.append(f)
    return f


def objfungrad(x):
    gradf = np.array(
        [
            np.exp(-0.1 * (x[0] ** 2 + x[1] ** 2)) * 0.2 * x[0],
            np.exp(-0.1 * (x[0] ** 2 + x[1] ** 2)) * 0.2 * x[1],
        ]
    )
    return gradf


# Define constraints


def confun(x):
    g1 = x[0] + x[1] - 2
    g2 = x[1] - 1

    g1list.append(g1)
    g2list.append(g2)
    return g1, g2


def confungrad(x):
    gradg1 = np.array([1, 1])
    gradg2 = np.array([0, 1])
    return gradg1, gradg2


# Create list of constraints

con1 = {"type": "ineq", "fun": confun, "jac": confungrad}

cons = [con1]

# Define starting point
x0 = np.array([3.0, 2.0])

# Define bounds
bounds = [[-4.0, 4.0], [-3.0, 3.0]]  # x[0] bounds  # x[1] bounds

# Run optimizer
result = minimize(
    objfun, x0, jac=objfungrad, constraints=cons, bounds=bounds, method="slsqp"
)

# Print results
print(result)
xopt = result.x

# Plot optimization history
fig = plt.figure()
plt.subplot(311)
plt.plot(xlist, "o-")
plt.ylabel("x", fontsize=20)
plt.subplot(312)
plt.plot(flist, "o-")
plt.ylabel("f", fontsize=20)
plt.subplot(313)
plt.plot(g1list, "o-")
plt.plot(g2list, "o-")
plt.plot([0, len(g1list) - 1], [0, 0], "gray", linewidth=0.5)
plt.ylabel("g", fontsize=20)
plt.xlabel("evaluations", fontsize=20)
plt.tight_layout()

plt.show()
