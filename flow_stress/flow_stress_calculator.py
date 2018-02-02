import numpy as np
from scipy.constants import convert_temperature
import math
from scipy import optimize as opt  # for optimization
import matplotlib.pyplot as plt
from matplotlib import cm
from fugacity_calculator import *

# CALCULATE FLOW STRESS
# Takes imputs of pressure and temperature converts them from MPa and C to
# Pa and K, calculates fugacity, differential stress, flow stress and slip rate
# and plots the results. 
# Pa and K

class FlowStressCalculator():

    def __init__(self, temperature_values, pressure_values):

        self.temperature = convert_temperature(temperature_values, 'c','k')#C2K(np.array(temperature_values))

        self.pressure = np.array([pressure_values])*1.0E6
        self.grain_size = []
        self.fugacity = []
        self.differential_stress = []
        self.strain_rate = []
        self.slip_rate = []
        self.width = []

    def calculate_fugacity(self):

        # for t, p in zip(self.temperature, self.pressure):
        for t in self.temperature:
            for p in self.pressure:
                calculate_coefficient_table(t)
                fug = fugacity_optimizer(t, p)
                self.fugacity.append(fug)
        return self.fugacity

    def calculate_differential_stress(self, grain_size, paleopiezometer='ST03'):

        self.grain_size = grain_size
        for grain in self.grain_size:
            part = math.exp((math.log(grain)-math.log(PIEZOMETERS[paleopiezometer][
                            'Constant_B']))/PIEZOMETERS[paleopiezometer]['Exponent'])
            self.differential_stress.append(part)
        return self.differential_stress

    def calculate_strain_rate(self, flow_law='H01'):

        for stress in self.differential_stress:
            for t, f in zip(self.temperature, self.fugacity):
                sr_i = (FLOW_LAWS[flow_law]['A']*np.power(stress, FLOW_LAWS[flow_law]['n'])
                        * np.power(f, 1)*np.exp(-FLOW_LAWS[flow_law]['Q']/(8.3144598*t)))
                self.strain_rate.append(sr_i)

        return self.strain_rate

    def calculate_slip_rate(self, width):  # width in m, output of mm/yr
        
        self.width = width
        width = [width]
        for w in width:
            for strain in self.strain_rate:
                vel = w*1000*31536000*strain
                self.slip_rate.append(vel)

        return self.slip_rate

    def plot_strain_slip_rates(self):

        temperature_C = convert_temperature(self.temperature,'k','c')

        # Group the data for plotting
        sr_grouped = [self.strain_rate[
            x:x+len(self.temperature)] for x in range(0, len(self.strain_rate), len(self.temperature))] #xrange
        sl_grouped = [self.slip_rate[
            x:x+len(self.temperature)] for x in range(0, len(self.slip_rate), len(self.temperature))] #xrange

        fig = plt.figure()
        sub1 = fig.add_subplot(2, 1, 1)
        sub1.set_yscale('log')
        sub1.text(0.01, 0.85, '[A] Strain Rates',
                  transform=sub1.transAxes, fontsize=10)
        sub1.set_ylabel('Strain Rate (1/s)')
        sub1.legend(loc='upper left')

        sub2 = fig.add_subplot(2, 1, 2, sharex=sub1)
        sub2.set_yscale('log')
        sub2.text(0.01, 0.85, '[B] Slip Rates',
                  transform=sub2.transAxes, fontsize=10)
        sub2.text(0.75, 0.15, 'Width: %s m' % self.width, transform=sub2.transAxes, fontsize=10)
        sub2.set_xlabel('Temperature ($\degree$C)')
        sub2.set_ylabel("Velocity (mm/yr)")

        plt.setp(sub1.get_xticklabels(), visible=False)
        plt.subplots_adjust(hspace=0.2)

        colors = iter(cm.jet(np.linspace(0, 1, len(self.grain_size))))
        colors2 = iter(cm.jet(np.linspace(0, 1, len(self.grain_size))))
        for i, grain in enumerate(self.grain_size):
            sub1.plot(temperature_C, sr_grouped[i], color=next(
                colors), label=(str(grain) + ' $\mu$m'))
            sub1.legend(loc='center left', bbox_to_anchor=(1, 0),
                        title="Grain Size", fontsize='x-small')
            sub2.plot(temperature_C, sl_grouped[i], color=next(colors2))

        sub1.set_xlim(left=min(temperature_C), right=max(temperature_C))
        sub2.set_xlim(left=min(temperature_C), right=max(temperature_C))
    
        plt.show()
        return fig



# def chunker(temperature, data):
#     chunk = [data[x:x+len(temperature)] for x in xrange(0, len(data), len(temperature))]

#     return chunk

# Export PDF file
# def export_pdf(fig, title='title'):
#     pdf = PdfPages(title+'.pdf')
#     pdf.savefig(fig)  # , bbox_inches='tight')
#     pdf.close()
