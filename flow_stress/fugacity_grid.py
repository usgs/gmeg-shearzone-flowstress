
#Functions to plot a grid of fugacity values for a range of temperatures and pressures.

import numpy as np
import matplotlib.pyplot as plt
from scipy.constants.constants import C2K, K2C
from scipy import optimize as opt # for optimization
from fugacity_calculator import *
from flow_stress_calculator import *



class FugacityGrid(FlowStressCalculator):

    def __init__(self, temperature, pressure):
        
        FlowStressCalculator.__init__(self, temperature, pressure)
        self.T, self.P = np.meshgrid(temperature, pressure)


    def fugacity_grid_optimizer(self, temperature, pressure):

        self.temperature = C2K(np.array(temperature))
        self.pressure = np.array(pressure)*1.0E6
        

        calculate_coefficient_table(self.temperature)

        def volume_limits(v):
            return eos(self.temperature, v)- self.pressure

        volume = opt.brentq(volume_limits, 5, 30) #Volume in cc/mo

        self.fugacity = PSfug(self.pressure, self.temperature, volume)#Calculate fugacity 
            
        return self.fugacity



    def fugacity_grid_plot(self):

        self.tC = K2C(self.temperature)
        self.pM = np.array(self.pressure)/1.0E6

        fo_v = np.vectorize(self.fugacity_grid_optimizer)
        result_array = fo_v(self.T, self.P)

        self.extent = [self.tC.min(), self.tC.max(), self.pM.max(), self.pM.min()]
        
        im = plt.imshow(result_array, cmap=plt.cm.jet, aspect="auto",extent=self.extent)#interpolation="none",
        plt.xlabel('Temperature (C)')
        plt.ylabel('Pressure (MPa)')
        plt.colorbar(im, label='Fugacity (MPa)')

        plt.show()






