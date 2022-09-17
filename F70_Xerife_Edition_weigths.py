# Análise de peso do F70 Xerife Edition

# -*- coding: utf-8 -*-
import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

airplane = dt.standard_airplane("F70_XerifeEdition")
dt.geometry(airplane)

W0_guess = 38000*9.81
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