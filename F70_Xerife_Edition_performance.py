# Análise de desempenho + Área de asa mínima do F70 Xerife Edition

import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

airplane_new = dt.standard_airplane("F70_XerifeEdition")
main_offset = 0
nlg_offset = 0.4
airplane_new["xr_w"] =  airplane_new["xr_w"] - main_offset
airplane_new["xcg_payload"] = 12.06
airplane_new["x_mlg"] = airplane_new["x_mlg"] - main_offset + 0.3
airplane_new["x_n"] = airplane_new["x_n"] - main_offset + 0.1
airplane_new["x_nlg"] = airplane_new["x_nlg"] + nlg_offset
## Q6-Q1

# Modify one parameter (if necessary)
SW = np.arange(85, 256, 1)
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
    W0_guess = 30000 * 9.81
    T0_guess = 0.3 * W0_guess
    dt.thrust_matching(W0_guess, T0_guess, airplane_new)

    W0_vec[count] = airplane_new["W0"]
    T0_vec[count, :] = airplane_new["T0vec"][:]
    deltaSWlan[count] = airplane_new["deltaS_wlan"]

    count = count + 1


AREA = np.interp(0, [deltaSWlan[27], deltaSWlan[28]], [SW[27], SW[28]])
PESO = np.interp(AREA, [SW[27], SW[28]], [W0_vec[27], W0_vec[28]])
print("Sw_min = ", AREA)
print("W0/Sw = ", PESO/AREA)


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
]

for j in m:

    plt.plot(W0_vec[:] / SW, T0_vec[:, j] / W0_vec[:], label=mylabels[j])

print("Q6-Q1")
plt.axvline(PESO / AREA, label=mylabels[8])
plt.legend()


# LIMITANTES: CRUISE, TAKEOFF AND LANDING
## Q6-Q3

countfig = countfig + 1

plt.figure(countfig)
plt.title("SW x W0")
plt.xlabel("SW")
plt.ylabel("W0")

print("Q6-Q3")
plt.plot(SW, W0_vec)

plt.show()