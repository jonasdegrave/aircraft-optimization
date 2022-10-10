# Análise de peso do F70 Xerife Edition

# -*- coding: utf-8 -*-
import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

airplane = dt.standard_airplane("F70_XerifeEdition")
dt.geometry(airplane)

# =============================================================================
# airplane["range_cruise"] =  4500 * 10**3
# airplane["xr_w"] = 12.82
# airplane["Lc_h"] = 3.7
# airplane["Lb_v"] = 0.42
# airplane["xcg_payload"] = 12.26
# airplane["x_mlg"] = airplane["xr_w"] + 4
# =============================================================================

W0_guess = 5e3 * airplane["S_w"]
T0_guess = 0.3*W0_guess
W0, We, Wf, Mf_cruise, xcg_e = dt.weight(W0_guess, T0_guess, airplane)


print('Weights in Kg (F70_XerifeEdition):')

print('W0:',W0/9.81)
print('We:',We/9.81)
print('Wf:',Wf/9.81)
print('Wp:',airplane['W_payload']/9.81)
print('Wc:',airplane['W_crew']/9.81)

print('breakdown:')
for key in ['W_w','W_h','W_v','W_f','W_nlg','W_mlg','W_eng','W_allelse']:
    print(key+': ',airplane[key]/9.81)
    
    pesos = np.array([airplane['W_payload'], airplane['W_crew'], Wf, airplane['W_w'], airplane['W_h'], airplane['W_v'], airplane['W_f'], airplane['W_mlg'], airplane['W_eng'], airplane['W_allelse'], airplane['W_nlg']]/W0)
    mycolors = ["#A9A9A9", "#8B008B", "#483D8B", "#8B0000", "#556B2F", "#FF8C00", "#008000", "#FFD700", "#FF69B4", "#FFA07A", "#A52A2A"]
    mylabels = ["Payload", "Crew", "Fuel", "Wing", "HS", "VS", "Fuselage", "MLG", "Engine", "Others", "NLG"]

plt.figure()
plt.pie(pesos, labels = mylabels, colors = mycolors, autopct='%1.1f%%')
plt.axis('equal')
plt.legend()