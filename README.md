# Flowstress

Developing a flow stress, strain rate and slip rate calculator and simulator for shear zones. 

Functions include calculators for strain rate that incorporate the temperature - pressure dependence of water fugacity at lithostatic pore fluid pressures for as many as nine seperate flow laws and four seperate paleopiezometers. The slip rate calculator takes variables of shear zone width and strain rate. 

PS.js courtesy of Tony Withers (http://publish.uwo.ca/~awither5/fugacity/index.htm) and is the basis for the fugacity depency in this code. 


Code also includes basic plotting functions, some examples here: 

<img width="1104" alt="ssr" src="https://user-images.githubusercontent.com/18178879/27503496-6b11ee40-5831-11e7-90f9-4a561189c6a5.png">
Plot of strain rate and slip rate over a range of temperatures, calculations account for variation in fugacity with varying temperature and a pressure of 400 MPa. Slip rates correspond to a shear zone 30 meters wide.


<img width="1091" alt="slr" src="https://user-images.githubusercontent.com/18178879/27503579-05d0ee72-5832-11e7-9f91-0917a21a3734.png">
Plot showing the four paleopiezometers included with the code.


# Examples
temperature = range(300, 600) #C
pressure = [400] #MPa
grain_size = range(5,25) #microns
f = FlowStressCalculator(temperature, pressure)
f.calculate_fugacity()
f.calculate_differential_stress(grain_size)
f.calculate_strain_rate()
f.calculate_slip_rate()

fig = plot_strain_slip_rates(f.temperature, strain_rate, slip_rate)





