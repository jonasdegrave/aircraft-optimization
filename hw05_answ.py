import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

# Load the standard airplane
airplane = dt.standard_airplane("fokker100")

# ######################### Q1

# Modify one parameter (if necessary)
SW = np.arange(80, 120, 1)
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

    airplane["S_w"] = i

    # Execute the geometry function from the designTools module (dt)
    dt.geometry(airplane)

    # Call the plotting function to make sure the aircraft is correct
    # dt.plot3d(airplane)

    # Execute the thrust on module
    W0_guess = 40000 * 9.81
    T0_guess = 0.3 * W0_guess
    dt.thrust_matching(W0_guess, T0_guess, airplane)

    W0_vec[count] = airplane["W0"]
    T0_vec[count, :] = airplane["T0vec"][:]
    deltaSWlan[count] = airplane["deltaS_wlan"]

    count = count + 1


AREA = np.interp(0, [deltaSWlan[1], deltaSWlan[2]], [SW[1], SW[2]])
PESO = np.interp(AREA, [SW[1], SW[2]], [W0_vec[1], W0_vec[2]])

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


plt.axvline(PESO / AREA, label=mylabels[8])

# OS LIMITANTES SÃO: SECOND SEGMENT CLIMB PARA BAIXA CARGA ALAR, TAKEOFF PARA CARGAS ALARES CRESCENTES ATÉ QUE O LIMITANTE PASSA A SER LANDING.

# ######################### Q2

MTOW = (W0_vec[14] + W0_vec[13]) / 2

SCATTER_SW = 93.5
SCATTER_T0 = 120500
plt.scatter(MTOW / SCATTER_SW, SCATTER_T0 / MTOW)
plt.legend()

# ######################### Q3

countfig = countfig + 1

plt.figure(countfig)
plt.title("SW x W0")
plt.xlabel("SW")
plt.ylabel("W0")

plt.plot(SW, W0_vec)

############################ Q6

airplane_new = dt.standard_airplane("AviaoDoXerife")

## Q6-Q1

# Modify one parameter (if necessary)
SW = np.arange(80, 120, 1)
m = np.arange(8)

W0_vec = np.zeros(len(SW))
T0_vec = np.zeros((len(SW), len(m)))
deltaSWlan = np.zeros(len(SW))

count = 0
countfig += 1

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
    W0_guess = 40000 * 9.81
    T0_guess = 0.3 * W0_guess
    dt.thrust_matching(W0_guess, T0_guess, airplane_new)

    W0_vec[count] = airplane_new["W0"]
    T0_vec[count, :] = airplane_new["T0vec"][:]
    deltaSWlan[count] = airplane_new["deltaS_wlan"]

    count = count + 1


AREA = np.interp(0, [deltaSWlan[1], deltaSWlan[2]], [SW[1], SW[2]])
PESO = np.interp(AREA, [SW[1], SW[2]], [W0_vec[1], W0_vec[2]])

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

## Q6-Q3

countfig = countfig + 1

plt.figure(countfig)
plt.title("SW x W0")
plt.xlabel("SW")
plt.ylabel("W0")

print("Q6-Q3")
plt.plot(SW, W0_vec)

plt.show()
