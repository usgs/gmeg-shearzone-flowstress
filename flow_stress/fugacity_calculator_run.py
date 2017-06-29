#Run fugacity calculator class
from flow_stress.fugacity_calculator import *

t = 723
p = 400000000


f = FugacityCalculator(t,p)
f.calculate_coefficient_table()
f.fugacity_optimizer()