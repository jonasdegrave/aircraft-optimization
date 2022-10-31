import numpy as np
from pymoo.core.problem import ElementwiseProblem

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

from pymoo.termination import get_termination

from pymoo.optimize import minimize

import designTool as dt

class MyProblem(ElementwiseProblem):

    def __init__(self):
        
        self.airplane = dt.standard_airplane("F70_XerifeEdition")
        
        Sweep_Lower = ( 16.5 * np.pi / 180 ) * 0.4
        Sweep_Upper = ( 16.5 * np.pi / 180 ) * 1.6
        Dihedral_Lower = ( 3.3 * np.pi / 180 ) * 0.4
        Dihedral_Upper = ( 3.3 * np.pi / 180 ) * 1.6
        Range_Lower = ( 4700 * 10**3 ) * 0.4
        Range_Upper = ( 4700 * 10**3 ) * 1.6
        
        super().__init__(n_var=3,
                         n_obj=1,
                         n_ieq_constr=12,
                         xl=np.array([ Sweep_Lower, Dihedral_Lower, Range_Lower ]),
                         xu=np.array([ Sweep_Upper, Dihedral_Upper, Range_Upper ]))

    def _evaluate(self, x, out, *args, **kwargs):
        
        self.airplane["range_cruise"] = x[0]
        self.airplane["dihedral_w"] = x[1]
        self.airplane["sweep_w"] = x[2]
        
        airplane = dt.analyze(
            airplane=self.airplane,
            print_log=False,  # Plot results on the terminal screen
            plot=False,  # Generate 3D plot of the aircraft
        )
        
        
        f1 = airplane["range_cruise"]
        
        g1 =    airplane["W0"]  - 39000*9.81
        g2 = -  airplane["deltaS_wlan"]
        g3 =    airplane["SM_fwd"] - 0.3
        g4 = - (airplane["SM_aft"] - 0.05)
        g5 =    airplane["frac_nlg_fwd"] - 0.18
        g6 = - (airplane["frac_nlg_aft"] - 0.05)
        g7 = - (airplane["alpha_tipback"] * (180.0 / np.pi) - 15)
        g8 = - (airplane["alpha_tailstrike"]  * (180.0 / np.pi) - 10)
        g9 =   (airplane["phi_overturn"] * (180.0 / np.pi) - 63)
        g10 =   airplane["b_tank_b_w"] - 0.95
        g11 =   airplane["b_w"] - 36
        g12 =   airplane["CLv"] - 0.75
        
        out["F"] = [f1]
        out["G"] = [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12]


problem = MyProblem()

algorithm = NSGA2(
    pop_size=40,
    n_offsprings=10,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=True
)

termination = get_termination("n_gen", 40)

res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               save_history=True,
               verbose=True)

X = res.X
F = res.F