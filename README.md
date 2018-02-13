# Flowstress Calculator

Developing a flow stress, strain rate and slip rate calculator and simulator for shear zones. 

Functions include calculators for strain rate that incorporate the temperature - pressure dependence of water fugacity at lithostatic pore fluid pressures for as many as nine seperate flow laws and four seperate paleopiezometers. The slip rate calculator takes variables of shear zone width and strain rate. PS.js courtesy of Tony Withers (http://publish.uwo.ca/~awither5/fugacity/index.htm) and is the basis for the fugacity depency in this code. 

Code also includes basic plotting functions:

Plot of strain and slip rates for a range of temperatures at 400 MPa using the Hirth et al. (2001) flow law. 

![screen shot 2017-06-26 at 9 34 12 am](https://user-images.githubusercontent.com/18178879/27549955-b550f94c-5a52-11e7-900c-9b20ff36f156.png)

Code currently includes four seperate paleopiezometers, default is Stipp and Tullis (2003).

![screen shot 2017-06-26 at 9 20 35 am](https://user-images.githubusercontent.com/18178879/27549580-3e47df88-5a51-11e7-89a7-a1103a3b4af3.png)


Fugacity as a function of pressure and temperature:

![screen shot 2017-06-29 at 6 09 17 pm](https://user-images.githubusercontent.com/18178879/27716861-1e7478ea-5cf6-11e7-9ab5-bdaef92f89bf.png)


# Examples
## Import and define pressure, temperature conditions
```
from flow_stress.flow_stress_calculator import FlowStressCalculator
from flow_stress.pt_conditions import *
from flow_stress.fugacity_grid import *

# Decide what pressure temperature ranges you want to use, 
# takes depth in km, density in g/cm3 and geothermal gradient in C/km. 
# Outputs are: pressure, temperature
print(PTCalculator(10, 2.7, 30).pt_calculator()) 
print(PTCalculator(18, 2.7, 30).pt_calculator())
```

## Plot the fugacity grid

```
# Use these ranges to define your pressure temperature ranges
pressure = range(260,480,10)
temperature = range(300,540,10)

FugacityGrid(temperature,pressure).fugacity_grid_plot() #Creates fugacity plot
```

## Calculate flow stress and plot

```
# Define grain size (microns) and shear zone width (meters) and a single pressure to calculate strain and slip rates
grain_size = range(5,26,2)
width = 20 #Width in meters
pressure = [400] # Plot only works at a single pressure over a range of temperatures

f = FlowStressCalculator(temperature, pressure)
fugacity = f.calculate_fugacity()
differential_stress = f.calculate_differential_stress(grain_size)  # Default is Holyoke and Kronenberg (2010)
strain_rate = f.calculate_strain_rate() # Default flow law is Hirth et al. (2001)
slip_rate = f.calculate_slip_rate(width) 
f.plot_strain_slip_rates()
```

## To explore options

```
#To see all piezometers and flow laws
from flow_stress.fugacity_calculator import *
print(FLOW_LAWS)
print(PIEZOMETERS)
```

## Run Testsuite
```pip install nosetests```

```nosetests```
