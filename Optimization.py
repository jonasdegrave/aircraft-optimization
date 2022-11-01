###############################################################################
########################## IMPORTANDO AS BIBLIOTECAS ##########################
###############################################################################

# Numpy
import numpy as np

# ElementwiseProblem usado para definir o nosso problema
from pymoo.core.problem import ElementwiseProblem

# Bibliotecas para o algoritmo genético utilizado na otimização
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

# Biblioteca para o critério de parada
from pymoo.termination import get_termination
from pymoo.termination.default import DefaultSingleObjectiveTermination

# Biblioteca para o minimizador
from pymoo.optimize import minimize

# DesignTools
import designTool as dt

# Matplotlib
import matplotlib.pyplot as plt

###############################################################################
############################## FUNÇÕES PRÓPRIAS ###############################
###############################################################################

# Função para normalizar os parâmetros e alimentar ao otimizador.
def Normalize(x, xmin, xmax):

    y = (x - xmin) / (xmax - xmin)

    return y


# Função para voltar os parâmetros ao domínio original
# e alimentar ao dicionário.
def Denormalize(y, xmin, xmax):

    x = y * (xmax - xmin) + xmin

    return x


###############################################################################
################################# OTIMIZADOR ##################################
###############################################################################

########------------------ DEFININDO O NOSSO PROBLEMA ------------------#######


class MyProblem(ElementwiseProblem):
    def __init__(self):

        self.airplane = dt.standard_airplane("F70_XerifeEdition")

        # 11 variáveis de entrada
        # 1 função objetivo
        # 16 restrições

        # Variáveis de entrada estão normalizadas entre 0 e 1.

        super().__init__(
            n_var=11,
            n_obj=1,
            n_ieq_constr=16,
            xl=np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            xu=np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        )

    def _evaluate(self, x, out, *args, **kwargs):

        ####### VARIÁVEIS

        # Enflechamento
        Sweep_Lower = 6.6 * np.pi / 180
        Sweep_Upper = 35 * np.pi / 180

        # Mach de cruzeiro
        MachCruise_Lower = 0.75
        MachCruise_Upper = 0.80

        # Alongamento
        AR_Lower = 6
        AR_Upper = 11

        # Área de asa
        SW_Lower = 60
        SW_Upper = 130

        # Alcance em metros
        Range_Lower = 4500000
        Range_Upper = 5500000

        # Altitude de cruzeiro em metros
        AltCruise_Lower = 11000
        AltCruise_Upper = 14000

        # Posição longitudinal da asa
        xr_w_Lower = 10
        xr_w_Upper = 15

        # Coeficiente de volume da empenagem horizontal
        Cht_Lower = 0.75
        Cht_Upper = 1.1

        # Braço de alavanca adimensional da empenagem horizontal
        Lch_Lower = 3.7 * 0.9
        Lch_Upper = 3.7 * 1.1

        # Coeficiente de volume da empenagem vertical
        Cvt_Lower = 0.04
        Cvt_Upper = 0.12

        # Braço de alavanca adimensional da empenagem vertical
        Lbv_Lower = 0.42 * 0.9
        Lbv_Upper = 0.42 * 1.1

        ####### Retornando as variáveis para o domínio original de ordem de grandeza.

        x[0] = Denormalize(x[0], Sweep_Lower, Sweep_Upper)
        x[1] = Denormalize(x[1], MachCruise_Lower, MachCruise_Upper)
        x[2] = Denormalize(x[2], AR_Lower, AR_Upper)
        x[3] = Denormalize(x[3], SW_Lower, SW_Upper)
        x[4] = Denormalize(x[4], Range_Lower, Range_Upper)
        x[5] = Denormalize(x[5], AltCruise_Lower, AltCruise_Upper)
        x[6] = Denormalize(x[6], xr_w_Lower, xr_w_Upper)
        x[7] = Denormalize(x[7], Cht_Lower, Cht_Upper)
        x[8] = Denormalize(x[8], Lch_Lower, Lch_Upper)
        x[9] = Denormalize(x[9], Cvt_Lower, Cvt_Upper)
        x[10] = Denormalize(x[10], Lbv_Lower, Lbv_Upper)

        ####### Atualizando o dicionário com o valor das variáveis na geração atual.

        self.airplane["sweep_w"] = x[0]
        self.airplane["Mach_cruise"] = x[1]
        self.airplane["AR_w"] = x[2]
        self.airplane["S_w"] = x[3]
        self.airplane["range_cruise"] = x[4]
        self.airplane["altitude_cruise"] = x[5]
        self.airplane["xr_w"] = x[6]
        self.airplane["x_mlg"] = x[6] + 4.28
        self.airplane["Cht"] = x[7]
        self.airplane["Lc_h"] = x[8]
        self.airplane["Cvt"] = x[9]
        self.airplane["Lb_v"] = x[10]

        ####### Realizando o 'analyze' com os valores da geração atual.

        airplane = dt.analyze(
            airplane=self.airplane,
            print_log=False,  # Plot results on the terminal screen
            plot=False,  # Generate 3D plot of the aircraft
        )

        ####### Definindo a função objetivo a ser minimizada.

        f1 = -airplane["range_cruise"]

        ####### Definindo as restrições a serem respeitadas pelo otimizador.

        g1 = (airplane["W0"] / (39000 * 9.81)) - 1
        g2 = -airplane["deltaS_wlan"]
        g3 = airplane["SM_fwd"] / 0.3 - 1
        g4 = -(airplane["SM_aft"] / 0.05 - 1)
        g5 = airplane["frac_nlg_fwd"] / 0.18 - 1
        g6 = -(airplane["frac_nlg_aft"] / 0.05 - 1)
        g7 = -((airplane["alpha_tipback"] * (180.0 / np.pi) / 15) - 1)
        g8 = -((airplane["alpha_tailstrike"] * (180.0 / np.pi) / 10) - 1)
        g9 = (airplane["phi_overturn"] * (180.0 / np.pi) / 63) - 1
        g10 = airplane["b_tank_b_w"] / 0.95 - 1
        g11 = airplane["b_w"] / 36 - 1
        g12 = airplane["CLv"] / 0.75 - 1
        g13 = airplane["T0"] / 130000 - 1
        g14 = (airplane["xr_v"] + airplane["cr_v"]) / (airplane["L_f"] - 0.5) - 1
        g15 = abs(airplane["xr_h"] - airplane["xt_v"]) / 0.1 - 1
        g16 = abs(airplane["cr_h"] - airplane["ct_v"]) / 0.5 - 1

        ####### Saídas do otimizador.

        out["F"] = [f1]
        out["G"] = [
            g1,
            g2,
            g3,
            g4,
            g5,
            g6,
            g7,
            g8,
            g9,
            g10,
            g11,
            g12,
            g13,
            g14,
            g15,
            g16,
        ]


####### O problema está agora definido.
problem = MyProblem()

####### Definindo o algoritmo genético utilizado e os seus parâmetros:
# Tamanho da população inicial: 200;
# Quantidade de filhos por geração: 80;
# Possibilidade de crossover com 90% de chance de ocorrência;
# Possibilidade de mutação;
# Eliminação de dados duplicados.

algorithm = NSGA2(
    pop_size=200,
    n_offsprings=80,
    sampling=FloatRandomSampling(),
    crossover=SBX(prob=0.9, eta=15),
    mutation=PM(eta=20),
    eliminate_duplicates=True,
)

####### Definindo o critério de parada.
#  Critério: parada na geração 200.

termination = DefaultSingleObjectiveTermination(
    xtol=1e-5, ftol=50, period=100, n_max_gen=2000, n_max_evals=100000
)
####### Compilando o otimizador com a definição do problema, algoritmo e
#######critério de parada.

res = minimize(problem, algorithm, termination, seed=3, save_history=True, verbose=True)

X = res.X


F = res.F

n_evals = np.array([e.evaluator.n_eval for e in res.history])
opt = np.array([e.opt[0].F for e in res.history])

plt.title("Convergence")
plt.plot(n_evals, -opt, "--")
plt.show()


Sweep_Lower = 6.6 * np.pi / 180
Sweep_Upper = 35 * np.pi / 180

MachCruise_Lower = 0.75
MachCruise_Upper = 0.80

AR_Lower = 6
AR_Upper = 11

SW_Lower = 60
SW_Upper = 130

Range_Lower = 4500000
Range_Upper = 5500000

AltCruise_Lower = 11000
AltCruise_Upper = 14000

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

Sweep = (Denormalize(X[0], Sweep_Lower, Sweep_Upper)) * (180 / np.pi)
MachCruise = Denormalize(X[1], MachCruise_Lower, MachCruise_Upper)
AR_w = Denormalize(X[2], AR_Lower, AR_Upper)
S_w = Denormalize(X[3], SW_Lower, SW_Upper)
Range = Denormalize(X[4], Range_Lower, Range_Upper)
AltCruise = Denormalize(X[5], AltCruise_Lower, AltCruise_Upper)
xr_w = Denormalize(X[6], xr_w_Lower, xr_w_Upper)
Cht = Denormalize(X[7], Cht_Lower, Cht_Upper)
Lch = Denormalize(X[8], Lch_Lower, Lch_Upper)
Cvt = Denormalize(X[9], Cvt_Lower, Cvt_Upper)
Lbv = Denormalize(X[10], Lbv_Lower, Lbv_Upper)
T0 = 130000 + res.G[12] * 130000
