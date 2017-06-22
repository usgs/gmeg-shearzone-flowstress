# import numpy as np
# from scipy.constants.constants import C2K
from fugacity import *
import unittest 
import pdb

temperature = [450, 500]
pressure = [400]
grain_size =[10]
width = [30]




class Test_flow_stress_calculator(unittest.TestCase):


	def test_if_flow_stress_calculator_initializes(self):
		self.fs_calculation = FlowStressCalculator(temperature, pressure)

		self.assertEqual(self.fs_calculation.temperature[1], 773.15)#723.15
		self.assertEqual(self.fs_calculation.pressure, 400000000)
		self.assertEqual(self.fs_calculation.grain_size, [])
		self.assertEqual(self.fs_calculation.fugacity,[])
		self.assertEqual(self.fs_calculation.differential_stress,[])
		self.assertEqual(self.fs_calculation.strain_rate,[])

	def test_if_flow_stress_calculator_works(self):
		self.fs_calculation = FlowStressCalculator(temperature, pressure)

		self.fs_calculation.calculate_fugacity()
		self.fs_calculation.calculate_differential_stress(grain_size) 
		self.sr = self.fs_calculation.calculate_strain_rate()
		self.slip = self.fs_calculation.calculate_slip_rate(width)
		

		self.assertEqual(self.fs_calculation.fugacity,[116.71651725233453, 149.72095732050246])
		self.assertEqual(self.fs_calculation.differential_stress, [78.75972939002484])
		self.assertEqual(self.sr, [5.0181530333149198e-12, 2.7498290139749202e-11])
		self.assertEqual(self.slip, [4.7475742217585797, 26.015582335413924])
		# pdb.set_trace()

		

		