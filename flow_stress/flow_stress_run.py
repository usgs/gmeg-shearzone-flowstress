
from flow_stress.flow_stress import *#FlowStressCalculator, plot_strain_slip_rates


temp = range(300, 600)
pressure = [400]
grain_size = range(3, 35, 3)#range(5,10)#[10,12,13,16,27,29]##[5,6,7,8,9,10,11,12,13,14,15,30] #List of grainsizes 
width = [30]


f = FlowStressCalculator(temp,pressure)
fugacity = f.calculate_fugacity()

differential_stress = f.calculate_differential_stress(grain_size)#, paleopiezometer = 'T77')

strain_rate= f.calculate_strain_rate()

slip_rate = f.calculate_slip_rate(width)


#convert K back to C for plotting
temperature_C = K2C(f.temperature)


plot_strain_slip_rates(f.temperature, strain_rate, slip_rate, grain_size)