import numpy as np
import matplotlib.pyplot as plt
from scipy.constants.constants import C2K, K2C
from scipy import optimize as opt # for optimization
from fugacity_calculator import *


#t = np.array([723.15, 773.15, 823.15, 873.15])#450, 500, 550 C
# t = np.arange(600, 800, 5)
# tC = K2C(t)
# #p = np.array([400000000, 450000000, 500000000, 550000000])
# p = np.arange(200000000, 700000000, 5000000)
# T, P = np.meshgrid(t,p)

class FugacityGrid():

    def __init__(self, temperature, pressure):
        self.tC = temperature
        self.pM = pressure
        self.temperature = C2K(np.array(temperature))
        self.pressure = np.array(pressure)*1.0E6
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

        #T, P = np.meshgrid(self.temperature, self.pressure)

        fo_v = np.vectorize(self.fugacity_grid_optimizer)
        result_array = fo_v(self.T, self.P)

        self.extent = [self.tC.min(), self.tC.max(), self.pM.max(), self.pM.min()]
        
        im = plt.imshow(result_array, cmap=plt.cm.jet, aspect="auto",extent=self.extent)#interpolation="none",
        plt.xlabel('Temperature (C)')
        plt.ylabel('Pressure (MPa)')
        plt.colorbar(im, label='Fugacity (MPa)')

        plt.show()






