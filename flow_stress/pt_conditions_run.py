from pt_conditions import *

depth = [20,14]
density = 2.7
geothermal_gradient = 30



pt = PTCalculator(depth, density, geothermal_gradient)

p,t  = pt.pt_calculator()
print p,t
