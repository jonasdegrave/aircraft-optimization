# Análise de desempenho + Área de asa mínima do F70 Xerife Edition

import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

airplane_new = dt.standard_airplane("F70_XerifeEdition")
#airplane_new["altitude_cruise"]= 11000  # Cruise altitude [m]
dt.geometry(airplane_new)

## Q6-Q1

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

AREA = np.interp(0, [deltaSWlan[5], deltaSWlan[6]], [SW[5], SW[6]])
PESO = np.interp(AREA, [SW[5], SW[6]], [W0_vec[5], W0_vec[6]])
print("Sw_min = ", AREA)
#print("W0/Sw = ", PESO/AREA)


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
    "Our project"
]

for j in m:

    plt.plot(W0_vec[:] / SW, T0_vec[:, j] / W0_vec[:], label=mylabels[j])

airplane_aux = dt.standard_airplane("F70_XerifeEdition")
dt.analyze(airplane_aux)
W0_over_Sw = airplane_aux["W0"]/airplane_aux["S_w"]
T0_over_W0 = airplane_aux["T0"]/airplane_aux["W0"]
plt.axvline(PESO / AREA, label=mylabels[8])
plt.scatter(W0_over_Sw, T0_over_W0, label=mylabels[9])
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