from pt_conditions import *
from fugacity_calculator import *
from flow_stress import *
from fugacity_grid import *


#User inputs
#Depth in km
depth = range(10,21,1)
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
print p
print t

f = FlowStressCalculator(t, p)
f.calculate_fugacity()
f.calculate_differential_stress(grain_size) #As currently written grain sizes must be in a python list
sr= f.calculate_strain_rate()

sl = f.calculate_slip_rate(width)

f.plot_strain_slip_rates()

t2 = np.array(t)
p2 = np.array(p)
fg = FugacityGrid(t2,p2)
fg.fugacity_grid_plot()