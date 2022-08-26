# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 17:24:33 2022

@author: Microsoft Windows
"""

import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

###---------------------- Exercício 1 ----------------------###
# Load the standard airplane
airplane = dt.standard_airplane("fokker100")

# Modify one parameter
mach = np.arange(0.6, 0.9, 0.01)
sweep_w = [10, 15, 20, 25, 30, 35, 40]

cede = np.zeros([len(mach), len(sweep_w)])
CLmaximo = np.zeros([len(mach), len(sweep_w)])

for i in range(0, len(sweep_w)):

    for m in range(0, len(mach)):

        airplane["sweep_w"] = sweep_w[i] * np.pi / 180.0
        dt.geometry(airplane)

        CD, CLmax, dragDict = dt.aerodynamics(
            airplane,
            Mach=mach[m],
            altitude=11000,
            CL=0.5,
            W0_guess=40000 * 9.81,
            n_engines_failed=0,
            highlift_config="clean",
            lg_down=0,
            h_ground=0,
        )
        cede[m, i] = CD
        CLmaximo[m, i] = CLmax
# Call the plotting function to make sure the aircraft is correct
# dt.plot3d(airplane)

# Exercicio 1
plt.figure(1)
plt.title("Mach x $C_{D}$", fontsize="large")
plt.xlabel("Mach")
plt.ylabel("$C_{D}$")
for i in range(0, len(sweep_w)):

    plt.plot(mach, cede[:, i], label="$\Lambda$ = " + str(sweep_w[i]) + "$^\circ$")

plt.legend()

###---------------------- END Exercício 1 ----------------------###
###---------------------- Exercício 2 ----------------------###

CLmaximo = CLmaximo[0, :]

plt.figure(2)
plt.title("$\Lambda_{W}$ x $C_{L_{max}}$", fontsize="large")
plt.xlabel("$\Lambda_{W}$")
plt.ylabel("$C_{L_{max}}$")
plt.plot(sweep_w, CLmaximo)

###---------------------- END Exercício 2 ----------------------###
###---------------------- Exercício 3 ----------------------###

# TEÓRICO

###---------------------- END Exercício 3 ----------------------###
###---------------------- Exercício 4 ----------------------###

# Cruise: CLmax = linspace(-0.5, 1.5454459846664212)
# Takeoff: CLmax = linspace(-0.5, 2.1726437504368183)
# Landing: Clmax = linspace(-0.5, 2.5907755942837496)

airplane = dt.standard_airplane("fokker100")
dt.geometry(airplane)

CL_List = [
    np.linspace(-0.5, 1.5454459846664212),
    np.linspace(-0.5, 2.1726437504368183),
    np.linspace(-0.5, 2.5907755942837496),
]
CD_List = [
    np.zeros(len(CL_List[0])),
    np.zeros(len(CL_List[1])),
    np.zeros(len(CL_List[2])),
]

Mach_List = [0.75, 0.2, 0.2]
Altitude_List = [11000, 0, 0]
lgdown_List = [0, 1, 1]
hground_List = [0, 10, 10]
highliftconfig_List = ["clean", "takeoff", "landing"]

Label_List = ["Class 2: Cruise", "Class 2: Takeoff", "Class 2: Landing"]
plt.figure(3)
plt.title(
    "$C_{D}$ x $C_{L}$",
    fontsize="large",
)
plt.xlabel("$C_{D}$")
plt.ylabel("$C_{L}$")

for j in range(len(CL_List)):

    count = 0
    for i in CL_List[j]:

        CD_List[j][count], CLmax, dragDict = dt.aerodynamics(
            airplane,
            Mach=Mach_List[j],
            altitude=Altitude_List[j],
            CL=i,
            W0_guess=40000 * 9.81,
            n_engines_failed=0,
            highlift_config=highliftconfig_List[j],
            lg_down=lgdown_List[j],
            h_ground=hground_List[j],
        )
        count = count + 1

    plt.plot(CD_List[j], CL_List[j], label=Label_List[j])

plt.legend()

###---------------------- END Exercício 4 ----------------------###
###---------------------- Exercício 5 ----------------------###

Swet_over_Sref_historical = 6.26
c_fe = 0.0026

############################ Cruise ############################

CD0_Cruise = Swet_over_Sref_historical * c_fe
AR = airplane["AR_w"]
e_Cruise = 0.85

K_Cruise = 1 / (AR * e_Cruise * np.pi)

CD_Cruise = CD0_Cruise + K_Cruise * CL_List[0] ** 2

############################ END Cruise ############################
############################ Takeoff ############################

CD0_Takeoff = CD0_Cruise + 0.015 + 0.02
e_Takeoff = 0.775
K_Takeoff = 1 / (AR * e_Takeoff * np.pi)

CD_Takeoff = CD0_Takeoff + K_Takeoff * CL_List[1] ** 2

############################ END Cruise ############################
############################ Landing ############################

CD0_Landing = CD0_Cruise + 0.065 + 0.02
e_Landing = 0.725
K_Landing = 1 / (AR * e_Landing * np.pi)

CD_Landing = CD0_Landing + K_Landing * CL_List[2] ** 2

############################ END Landing ############################
############################ Plotting ############################

plt.plot(CD_Cruise, CL_List[0], "--", label="Class 1: Cruise")
plt.plot(CD_Takeoff, CL_List[1], "--", label="Class 1: Takeoff")
plt.plot(CD_Landing, CL_List[2], "--", label="Class 1: Landing")

plt.legend()

#plt.show()


############################ End Plotting ############################
###---------------------- END Exercício 5----------------------###

###---------------------- Exercício 6 ----------------------###

# Preencher tabela com os resultados do 5 e do 4 (PRA UM AVIÃO PARECIDO COM FOKKER 100)

###---------------------- Exercício 7 ----------------------###

#print("RODOU ATÉ AQUI???")
airplane = dt.standard_airplane("AviaoDoXerife")
dt.geometry(airplane)

CL_List = [
    np.linspace(-0.5, 1.5454459846664212),
    np.linspace(-0.5, 2.1726437504368183),
    np.linspace(-0.5, 2.5907755942837496),
]
CD_List = [
    np.zeros(len(CL_List[0])),
    np.zeros(len(CL_List[1])),
    np.zeros(len(CL_List[2])),
]

Mach_List = [0.75, 0.2, 0.2]
Altitude_List = [11000, 0, 0]
lgdown_List = [0, 1, 1]
hground_List = [0, 10, 10]
highliftconfig_List = ["clean", "takeoff", "landing"]

Label_List = ["Cruise", "Takeoff", "Landing"]
plt.figure(4)
plt.title(
    "$C_{D}$ x $C_{L}$",
    fontsize="large",
)
plt.xlabel("$C_{D}$")
plt.ylabel("$C_{L}$")

for j in range(len(CL_List)):

    count = 0
    for i in CL_List[j]:

        CD_List[j][count], CLmax, dragDict = dt.aerodynamics(
            airplane,
            Mach=Mach_List[j],
            altitude=Altitude_List[j],
            CL=i,
            W0_guess=40000 * 9.81,
            n_engines_failed=0,
            highlift_config=highliftconfig_List[j],
            lg_down=lgdown_List[j],
            h_ground=hground_List[j],
        )
        count = count + 1

    plt.plot(CD_List[j], CL_List[j], label=Label_List[j])

plt.legend()
plt.show()

# Preencher tabela com os resultados do 5 e do 4 (PRO NOSSO AVIÃO)
