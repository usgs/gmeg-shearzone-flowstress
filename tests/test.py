from fugacity import FugacityCalculator 
from strain_and_slip import calculate_differential_stress, calculate_strain_rate, calculate_slip_rate

#import pdb

temp = [450]
pressure = [400]
grain_size = [7]
width = [5]

f = FugacityCalculator(temp, pressure)
fc = f.calculate_fugacity()

s = calculate_differential_stress(grain_size)
sr = calculate_strain_rate(s, fc, flow_law='H01')
#pdb.set_trace()
v = calculate_slip_rate(sr, width)
#v = calculate_slip_rate_constant_strain_rate(5.01815303e-12, width)

print('P,T,F, Length: ' + str(len(fc)))
print(fc)
print('differential stress, Length: ' + str(len(s)))
print(s)
print('strain rate, Length: ' + str(len(sr)))
print(sr)
print('slip rate, Length: ' + str(len(v)))
print(v)

