from fugacity import FugacityCalculator 
from strain_and_slip import calculate_differential_stress, calculate_strain_rate, calculate_slip_rate

import pdb

temp = [400]
pressure = [400]
grain_size = 9.2
width = 30

f = FugacityCalculator(temp, pressure)
fc = f.calculate_fugacity()

s = calculate_differential_stress(grain_size)
sr = calculate_strain_rate(s, fc, flow_law='H01')
#pdb.set_trace()
v = calculate_slip_rate(sr, width)


print(fc)
print(s)
print(sr)
print(v)