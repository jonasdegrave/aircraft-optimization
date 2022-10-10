# Análise aerodinâmica do F70 Xerife Edition

# -*- coding: utf-8 -*-
import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

airplane = dt.standard_airplane("F70_XerifeEdition")
dt.geometry(airplane)

CL_List = [
    np.linspace(-0.5, 1.5532879704864728),
    np.linspace(-0.5, 2.170063659288774),
    np.linspace(-0.5, 2.5812474518236415),
]
CD_List = [
    np.zeros(len(CL_List[0])),
    np.zeros(len(CL_List[1])),
    np.zeros(len(CL_List[2])),
]

CLmax_List = [0, 0, 0]

Mach_List = [0.76, 0.2, 0.2]
Altitude_List = [11800, 0, 0]
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
            airplane,
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