# Flowstress

Developing a flow stress, strain rate and slip rate calculator and simulator for shear zones. 

Functions include calculators for strain rate that incorporate the temperature - pressure dependence of water fugacity at lithostatic pore fluid pressures for as many as nine seperate flow laws and four seperate paleopiezometers. The slip rate calculator takes variables of shear zone width and strain rate. 

PS.js courtesy of Tony Withers (http://publish.uwo.ca/~awither5/fugacity/index.htm) and is the basis for the fugacity depency in this code. 


Code also includes basic plotting functions, some examples here: 





![screen shot 2017-06-26 at 9 20 35 am](https://user-images.githubusercontent.com/18178879/27549580-3e47df88-5a51-11e7-89a7-a1103a3b4af3.png)


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




>>>>>>> 7d55d8d664a0ebaffda712f188478473395a7b09
