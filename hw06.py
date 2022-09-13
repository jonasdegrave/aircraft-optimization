import designTool as dt
import matplotlib.pyplot as plt
import numpy as np
import pprint

# Load the standard airplane
airplane = dt.standard_airplane("AviaoDoXerife")

# Modify one parameter (if necessary)
# airplane['S_w'] = 93.5

# Execute the geometry function from the designTools module (dt)
airplane = dt.analyze(
    airplane=airplane,
    print_log=True,  # Plot results on the terminal screen
    plot=True,  # Generate 3D plot of the aircraft
)
