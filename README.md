# Flowstress

Developing a flow stress, strain rate and slip rate calculator and simulator for shear zones. 

Functions include calculators for strain rate that incorporate the temperature - pressure dependence of water fugacity at lithostatic pore fluid pressures for as many as nine seperate flow laws and four seperate paleopiezometers. The slip rate calculator takes variables of shear zone width and strain rate. 

PS.js courtesy of Tony Withers (http://publish.uwo.ca/~awither5/fugacity/index.htm) and is the basis for the fugacity depency in this code. 


Code also includes basic plotting functions, some examples here: 

[Strain and Slip Rates](figs/Strain_Slip_Rates.pdf "Plot of strain rate and slip rate over a range of temperatures, calculations account for variation in fugacity with varying temperature and a pressure of 400 MPa. Slip rates correspond to a shear zone 30 meters wide.")


[Paleopiezometers](figs/Paleopiezometers.pdf "Plot showing the four paleopiezometers included with the code.")


### Examples
```
temperature = range(300, 600) #C
pressure = [400] #MPa
grain_size = range(5,25) #microns
width = [30] #m
f = FlowStressCalculator(temperature, pressure)
f.calculate_fugacity()
f.calculate_differential_stress(grain_size)
f.calculate_strain_rate()
f.calculate_slip_rate(width)

fig = plot_strain_slip_rates(f.temperature, strain_rate, slip_rate)
```




