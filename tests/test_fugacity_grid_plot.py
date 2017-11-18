import unittest
import sys

from flow_stress.fugacity_grid import *

depth  = range(10,21,1)
density = 2.7
geothermal_gradient = 30

# temperature = [450, 403]
# pressure = [400]




class test_Fugacity_Grid(unittest.TestCase):

	def test_if_pressure_temperature_values_are_created_correctly(self):
		self.pt = PTCalculator(depth, density, geothermal_gradient)
		self.p, self.t  = self.pt.pt_calculator()

		self.assertEqual(self.p[0], 264.6000000000001)
		self.assertEqual(self.t[0], 300)


	def test_if_fugacity_grid_initiates(self):

        self.T, self.P = np.meshgrid(temperature, pressure)

        self.assertEqual(self.P[0], 264.60000000000008)
        self.assertEqual(self.T[0], 300)

