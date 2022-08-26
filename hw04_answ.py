import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

# Load the standard airplane
airplane = dt.standard_airplane('fokker100')

# Modify one parameter (if necessary)
#airplane['AR_w'] = 8.43

# Execute the geometry function from the designTools module (dt)
dt.geometry(airplane)

# Call the plotting function to make sure the aircraft is correct
#dt.plot3d(airplane)

# Execute the weight estimation module
W0_guess = 40000*9.81
T0_guess = 0.3*W0_guess
W0, We, Wf, Mf_cruise, xcg_e = dt.weight(W0_guess, T0_guess, airplane)


print('Weights. in Newtons:')

print('W0:',W0)
print('We:',We)
print('Wf:',Wf)
print('Wp:',airplane['W_payload'])
print('Wc:',airplane['W_crew'])

print('breakdown:')
for key in ['W_w','W_h','W_v','W_f','W_nlg','W_mlg','W_eng','W_allelse']:
    print(key+': ',airplane[key])

countfig = 1

# Q1
pesos = np.array([airplane['W_payload'], airplane['W_crew'], Wf, airplane['W_w'], airplane['W_h'], airplane['W_v'], airplane['W_f'], airplane['W_mlg'], airplane['W_eng'], airplane['W_allelse'], airplane['W_nlg']]/W0)
mycolors = ["#A9A9A9", "#8B008B", "#483D8B", "#8B0000", "#556B2F", "#FF8C00", "#008000", "#FFD700", "#FF69B4", "#FFA07A", "#A52A2A"]
mylabels = ["Payload", "Crew", "Fuel", "Wing", "HS", "VS", "Fuselage", "MLG", "Engine", "Others", "NLG"]

plt.figure(countfig)
plt.pie(pesos, labels = mylabels, colors = mycolors, autopct='%1.1f%%')
plt.axis('equal')
plt.legend()
countfig = countfig+1

# Q2

AR = np.arange(6,12.5,0.5)
count = 0
W0_q2 = np.zeros(len(AR))

print(AR[10])

for i in AR:
    airplane['AR_w'] = i
    dt.geometry(airplane)
    W0_aux, We, Wf, Mf_cruise, xcg_e = dt.weight(W0_guess, T0_guess, airplane)
    W0_q2[count] = W0_aux
    count = count + 1

plt.figure(countfig)
plt.plot(W0_q2,AR)
countfig = countfig+1


plt.show()

