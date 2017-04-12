from fugacity import FugacityCalculator 
from strain_and_slip import calculate_differential_stress, calculate_strain_rate, calculate_slip_rate

#import pdb

temp = [450]
pressure = [400]
grain_size = [10]
width = 10

f = FugacityCalculator(temp, pressure)
fc = f.calculate_fugacity()

s = calculate_differential_stress(grain_size)
sr = calculate_strain_rate(s, fc, flow_law='H01')
#pdb.set_trace()
v = calculate_slip_rate(sr, width)
#v = calculate_slip_rate_constant_strain_rate(5.01815303e-12, width)


#print(fc)
print(s)
print('strain rate')
print(sr)
print('slip rate')
print(v)