import unittest
import sys

from flow_stress.pt_conditions import *

depth = 10
density = 2.7
geothermal_gradient = 30

class Test_pt_conditions(unittest.TestCase):

	def test_pt_conditions_initiates(self):
		self.pt = PTCalculator(depth, density, geothermal_gradient)

		self.assertEqual(self.pt.gravity, 9.8)
		self.assertEqual(self.pt.density, 2700)
		self.assertEqual(self.pt.depth, 10)
		self.assertEqual(self.pt.geothermal_gradient, 30)

	def test_if_pt_calculator_works(self):
		self.pt = PTCalculator(depth, density, geothermal_gradient)
		self.p, self.t = self.pt.pt_calculator()

		self.assertEqual(self.p[0], 264.6000000000001)
		self.assertEqual(self.t[0], 300)

	def test_if_pt_calculator_pressure_value_works(self):
		self.single_depth = 10
		self.pt = PTCalculator(depth, density, geothermal_gradient)
		self.p, self.t = self.pt.pt_calculator()
		self.single_p = self.pt.pt_calculator_pressure_value(self.p, self.single_depth)

		self.assertEqual(self.single_p, 264.6000000000001)

		#if single_pressure_value == False:
	#		print('You will have multiple pressure values')
