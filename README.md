# Flowstress

Developing a flow stress, strain rate and slip rate calculator and simulator for shear zones. 

Functions include calculators for strain rate that incorporate the temperature - pressure dependence of water fugacity at lithostatic pore fluid pressures for as many as nine seperate flow laws and four seperate paleopiezometers. The slip rate calculator takes variables of shear zone width and strain rate. PS.js courtesy of Tony Withers (http://publish.uwo.ca/~awither5/fugacity/index.htm) and is the basis for the fugacity depency in this code. 

Code also includes basic plotting functions:

Plot of strain and slip rates for a range of temperatures at 400 MPa using the Hirth et al. (2001) flow law. 

![screen shot 2017-06-26 at 9 34 12 am](https://user-images.githubusercontent.com/18178879/27549955-b550f94c-5a52-11e7-900c-9b20ff36f156.png)

Code currently includes four seperate paleopiezometers, default is Stipp and Tullis (2003) with the Holyoke and Kronenberg (2010) correction.

![screen shot 2017-06-26 at 9 20 35 am](https://user-images.githubusercontent.com/18178879/27549580-3e47df88-5a51-11e7-89a7-a1103a3b4af3.png)


Fugacity as a function of pressure and temperature:
![screen shot 2017-06-29 at 6 09 17 pm](https://user-images.githubusercontent.com/18178879/27716861-1e7478ea-5cf6-11e7-9ab5-bdaef92f89bf.png)


### Examples
#### Plotting strain and slip rate as a function of temperature and grain size
```
#User inputs 
temperature = range(300, 600) #C
pressure = [400] #MPa
grain_size = range(5,25) #microns
width = [30] #m


f = FlowStressCalculator(temperature, pressure) #Currently works for a range of temperature at a single pressure value.

#Cacluates the temperature-pressure dependence of water fugacity over this temperature range.
f.calculate_fugacity() 
f.calculate_differential_stress(grain_size) #default is Holyoke and Kronenberg (2010)
f.calculate_strain_rate() #Default flow law is Hirth et al. (2001)
f.calculate_slip_rate(width) #add with width of a shear zone

plot_strain_slip_rates(f.temperature, strain_rate, slip_rate)

```

#### Plotting fugacity as a function of pressure and temperature
```
#User inputs must be numpy arrays for fugacity grid plotter
t = np.arange(300, 600, 50)
p = np.arange(200, 700, 50)

fg = FugacityGrid(t, p)
fg.fugacity_grid_plot()
```

## Run Testsuite
```pip install nosetests```

```nosetests```
