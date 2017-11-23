from pt_conditions import *
from fugacity_calculator import *
from flow_stress_calculator import *
from fugacity_grid import *

import pdb

#User inputs
#Depth in km

depth = range(10,21,1) #10

#Density in gm/cc
density = 2.7
#Geothermal gradient in C/km
geothermal_gradient = 30

#Grain size range in microns
grain_size = range(10,21,1)
#Shear zone width in meters
width = 20 



pt = PTCalculator(depth, density, geothermal_gradient)

p,t  = pt.pt_calculator()

print(p)
print(t)

p = pt.pt_calculator_pressure_value(p,10) #Reassign pressure value to a single depth, rather than a range. (here we use 10km)

f = FlowStressCalculator(t, p)
f.calculate_fugacity() #only works for a single pressure value
f.calculate_differential_stress(grain_size) #As currently written grain sizes must be in a python list
sr= f.calculate_strain_rate()

sl = f.calculate_slip_rate(width)

f.plot_strain_slip_rates()

t2 = np.array(t)#Pressure and temperature must be in numpy arrays for grid plots
p2 = np.array(p)#only works over a range of pressures

pdb.set_trace()

fg = FugacityGrid(t2,p2)
fg.fugacity_grid_plot()




