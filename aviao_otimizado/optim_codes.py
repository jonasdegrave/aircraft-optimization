"""
PROGRAMA DE ESPECIALIZAÇÃO EM ENGENHARIA AERONÁUTICA
OTIMIZAÇÃO MULTIDISCIPLINAR

Código com aplicações de ferramentas do Python para otimização
não restringida

Cap. Ney Sêcco 2020
"""

# IMPORTS
import numpy as np
from scipy.optimize import minimize, differential_evolution

# EXECUTION

# ==============================================

# Define number of input variables
nvar = 2

# Define objective function
def objfun(x):

    fun = np.sum(100 * (x[1:] - x[:-1] ** 2) ** 2 + (1 - x[:-1]) ** 2)

    return fun


def objfun_grad(x):

    grad = np.zeros(len(x))
    grad[0] = -400 * (x[1] - x[0] ** 2) * x[0] - 2 * (1 - x[0])
    grad[1:-1] = (
        -400 * (x[2:] - x[1:-1] ** 2) * x[1:-1]
        - 2 * (1 - x[1:-1])
        + 200 * (x[1:-1] - x[:-2] ** 2)
    )
    grad[-1] = 200 * (x[-1] - x[-2] ** 2)

    return grad


# ==============================================

### DIFFERENTIAL EVOLUTION

# Define bounds to determine initial population
bounds = [[-2, +2]] * nvar

# Solve the optimization problem
result = differential_evolution(
    objfun, bounds, polish=False, seed=1, atol=1e-10, maxiter=5000
)

# Print the results
print("")
print("DE")
print(result)

### NELDER-MEAD

# Define initial guess
x0 = np.zeros(nvar)

# Set optimization options
options = {"maxiter": 200000, "fatol": 1e-12, "adaptive": True}

# Run the optimization algorithm
result = minimize(objfun, x0, method="Nelder-Mead", options=options)

# Print results
print("")
print("NELDER-MEAD")
print(result)

### CG - Without gradient

# Set optimization options
options = {"maxiter": 200000}

# Run the optimization algorithm
result = minimize(objfun, x0, method="CG", tol=1e-6, options=options)

# Print results
print("")
print("CG without grad")
print(result)

### CG - With gradient

# Set optimization options
options = {"maxiter": 200000}

# Run the optimization algorithm
result = minimize(objfun, x0, jac=objfun_grad, method="CG", tol=1e-6, options=options)

# Print results
print("")
print("CG with grad")
print(result)

### BFGS - With gradient

# Set optimization options
options = {"maxiter": 200000}

# Run the optimization algorithm
result = minimize(objfun, x0, jac=objfun_grad, method="BFGS", tol=1e-6, options=options)

# Print results
print("")
print("BFGS with grad")
print(result)
