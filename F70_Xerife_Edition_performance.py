
# Análise de desempenho + Área de asa mínima do F70 Xerife Edition

import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

# ANÁLISE AERODINÂMICA

airplane_xerife_aero = dt.standard_airplane("F70_XerifeEdition_Optimized")
dt.geometry(airplane_xerife_aero)

CL_List = [
    np.linspace(-0.5, 1.5532879704864728),
    np.linspace(-0.5, 2.1595375198410105),
    np.linspace(-0.5, 2.57513395170655),
]
CD_List = [
    np.zeros(len(CL_List[0])),
    np.zeros(len(CL_List[1])),
    np.zeros(len(CL_List[2])),
]

CLmax_List = [0, 0, 0]

Mach_List = [0.76, 0.2, 0.2]
Altitude_List = [12623.454814339644, 0, 0]
lgdown_List = [0, 1, 1]
hground_List = [0, 10, 10]
highliftconfig_List = ["clean", "takeoff", "landing"]


Label_List = ["Cruise", "Takeoff", "Landing"]
plt.figure(1)
plt.title(
    "$C_{D}$ x $C_{L}$",
    fontsize="large",
)
plt.xlabel("$C_{D}$")
plt.ylabel("$C_{L}$")

for j in range(len(CL_List)):

    count = 0
    for i in CL_List[j]:

        CD_List[j][count], CLmax_List[j], dragDict = dt.aerodynamics(
            airplane_xerife_aero,
            Mach=Mach_List[j],
            altitude=Altitude_List[j],
            CL=i,
            W0_guess=38000 * 9.81,
            n_engines_failed=0,
            highlift_config=highliftconfig_List[j],
            lg_down=lgdown_List[j],
            h_ground=hground_List[j],
        )
        count = count + 1

    plt.plot(CD_List[j], CL_List[j], label=Label_List[j])

plt.legend()
plt.show()

## ANÁLISE DE DESEMPENHO

airplane_new = dt.standard_airplane("F70_XerifeEdition_Optimized")
airplane_old = dt.standard_airplane("F70_XerifeEdition")
dt.geometry(airplane_new)


SW = np.arange(75, 256, 1)
m = np.arange(8)

W0_vec = np.zeros(len(SW))
T0_vec = np.zeros((len(SW), len(m)))
deltaSWlan = np.zeros(len(SW))

count = 0
countfig = 1

plt.figure(countfig)
plt.title("T0/W0 x W0/Sw")
plt.xlabel("W0/Sw")
plt.ylabel("T0/W0")


for i in SW:

    airplane_new["S_w"] = i

    # Execute the geometry function from the designTools module (dt)
    dt.geometry(airplane_new)

    # Call the plotting function to make sure the aircraft is correct
    # dt.plot3d(airplane_new)

    # Execute the thrust on module
    W0_guess = 39000 * 9.81
    T0_guess = 0.3 * W0_guess
    dt.thrust_matching(W0_guess, T0_guess, airplane_new)

    W0_vec[count] = airplane_new["W0"]
    T0_vec[count, :] = airplane_new["T0vec"][:]
    deltaSWlan[count] = airplane_new["deltaS_wlan"]

    count = count + 1

AREA = np.interp(0, [deltaSWlan[6], deltaSWlan[7]], [SW[6], SW[7]])
PESO = np.interp(AREA, [SW[6], SW[7]], [W0_vec[6], W0_vec[7]])
print("Sw_min = ", AREA)


mylabels = [
    "Takeoff",
    "Cruise",
    "Takeoff Climb",
    "Transition Climb",
    "Second Climb",
    "En Route Climb",
    "Balked L. Climb",
    "Balked Climb OEI",
    "Landing",
    "Our project optimized",
    "Our project baseline"
]

plt.figure(2)
plt.title("T0/W0 x W0/Sw")
plt.xlabel("W0/Sw")
plt.ylabel("T0/W0")
for j in m:

    plt.plot(W0_vec[:] / SW, T0_vec[:, j] / W0_vec[:], label=mylabels[j])

airplane_xerife = dt.standard_airplane("F70_XerifeEdition_Optimized")
dt.analyze(airplane_xerife)
dt.analyze(airplane_old)
W0_over_Sw_opt = airplane_xerife["W0"] / airplane_xerife["S_w"]
T0_over_W0_opt = airplane_xerife["T0"] / airplane_xerife["W0"]
W0_over_Sw_base = airplane_old["W0"] / airplane_old["S_w"]
T0_over_W0_base = airplane_old["T0"] / airplane_old["W0"]
plt.axvline(PESO / AREA, label=mylabels[8])
plt.scatter(W0_over_Sw_opt, T0_over_W0_opt, label=mylabels[9])
plt.scatter(W0_over_Sw_base, T0_over_W0_base, label=mylabels[10])
plt.legend(loc = 'upper left')
plt.show()

# CÁLCULO DISTÂNCIA DE DECOLAGEM

TOP = (W0_over_Sw_opt/(CLmax_List[1]*T0_over_W0_opt))
S_TOFL = 0.2387*TOP

# CÁLCULO DISTÂNCIA DE POUSO

v_stall = np.sqrt(2*airplane_xerife['W0']*airplane_xerife['MLW_frac']/(1.225*CLmax_List[2]*airplane_xerife['S_w']))
v_a = 1.3*v_stall
S_FL = (v_a/1.701)**2
