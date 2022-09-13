import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

# Load the standard airplane
airplane = dt.standard_airplane('fokker100')

# Modify one parameter (if necessary)
airplane['S_w'] = 93.5

# Execute the geometry function from the designTools module (dt)
airplane = dt.analyze(airplane = None,
                      print_log = True, # Plot results on the terminal screen
                      plot = True, # Generate 3D plot of the aircraft
                      )