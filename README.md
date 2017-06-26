# Flowstress

Developing a flow stress, strain rate and slip rate calculator and simulator for shear zones. 

Functions include calculators for strain rate that incorporate the temperature - pressure dependence of water fugacity at lithostatic pore fluid pressures for as many as nine seperate flow laws and four seperate paleopiezometers. The slip rate calculator takes variables of shear zone width and strain rate. 

PS.js courtesy of Tony Withers (http://publish.uwo.ca/~awither5/fugacity/index.htm) and is the basis for the fugacity depency in this code. 


Code also includes basic plotting functions:

Plot of strain and slip rates for a range of temperatures at 400 MPa using the Hirth et al. (2001) flow law. 

![screen shot 2017-06-26 at 9 18 22 am](https://user-images.githubusercontent.com/18178879/27549623-5e2c3db2-5a51-11e7-84d4-2c4c6003f5d6.png)

Code currently includes four seperate paleopiezometers, default is Stipp and Tullis (2003) with the Holyoke and Kronenberg (2010) correction.

![screen shot 2017-06-26 at 9 20 35 am](https://user-images.githubusercontent.com/18178879/27549580-3e47df88-5a51-11e7-89a7-a1103a3b4af3.png)


### Examples
```
temperature = range(300, 600) #C
pressure = [400] #MPa
grain_size = range(5,25) #microns
width = [30] #m
f = FlowStressCalculator(temperature, pressure) #Currently works for a range of temperature at a single pressure value.
f.calculate_fugacity() #Cacluates the temperature-pressure dependence of water fugacity over this temperature range
f.calculate_differential_stress(grain_size) #default is Holyoke and Kronenberg (2010)
f.calculate_strain_rate() #Default flow law is Hirth et al. (2001)
f.calculate_slip_rate(width) #add with width of a shear zone

fig = plot_strain_slip_rates(f.temperature, strain_rate, slip_rate)
```


