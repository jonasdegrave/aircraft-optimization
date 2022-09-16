import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

# CONSTANTS
ft2m = 0.3048
kt2ms = 0.514444
lb2N = 4.44822
nm2m = 1852.0
gravity = 9.81
gamma_air = 1.4
R_air = 287

# Load the standard airplane
airplane = dt.standard_airplane("F70_XerifeEdition")

# dt.geometry(airplane)
# Modify one parameter (if necessary)
# airplane["range_cruise"] =  2000 * 10**3

# Execute the geometry function from the designTools module (dt)
airplane = dt.analyze(
    airplane=airplane,
    print_log=True,  # Plot results on the terminal screen
    plot=True,  # Generate 3D plot of the aircraft
)
print("")
print("AVISOS:")
print("")

if airplane["W0"] / gravity > 39000:
    print("MTOW EXCEDIDO")

if airplane["deltaS_wlan"] < 0:
    print("Área de asa muito pequena")

if airplane["SM_fwd"] > 0.3:
    print("Margem estática máxima (CG dianteiro) muito grande")

if airplane["SM_aft"] < 0.05:
    print("Margem estática mínima (CG traseiro) muito baixa")

if airplane["CLv"] > 0.75:
    print("Empenagem vertical muito pequena - condição de motor inoperante")

if airplane["frac_nlg_fwd"] > 0.18:
    print("Trem de pouso de nariz muito carregado - condição de CG dianteiro")

if airplane["frac_nlg_aft"] < 0.05:
    print("Trem de pouso de nariz pouco carregado - condição de CG traseiro")

if airplane["alpha_tipback"] * 180.0 / np.pi < 15:
    print("Perigo de tombar para trás na decolagem - tipback angle pequeno")

if airplane["alpha_tailstrike"] * 180.0 / np.pi < 10:
    print("Perigo de tailstrike")

if airplane["phi_overturn"] * 180.0 / np.pi > 63:
    print("Trem de pouso em configuração instável lateralmente")

if airplane["b_tank_b_w"] > 0.95:
    print("Volume de tanque de combustível insuficiente")
