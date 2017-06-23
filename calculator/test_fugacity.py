# import numpy as np
# from scipy.constants.constants import C2K
from fugacity import *
import unittest 
import pdb

temperature = [450, 403]
pressure = [400]
grain_size =[3, 35]
width = [30]




class Test_flow_stress_calculator(unittest.TestCase):


	def test_if_flow_stress_calculator_initializes(self):
		self.fs_calculation = FlowStressCalculator(temperature, pressure)

		self.assertEqual(self.fs_calculation.temperature[1], 676.14999999999998)#723.15
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
		

		self.assertEqual(self.fs_calculation.fugacity,[116.71651725233453, 87.86334635835435])
		self.assertEqual(self.fs_calculation.differential_stress, [204.77990936946233, 29.14098182806724])
		self.assertEqual(self.sr[0], [2.2933821401993624e-10])
		self.assertEqual(self.slip, [4.7475742217585797, 26.015582335413924])
		#pdb.set_trace()
		self.assertEqual(self.strain_rate, 0)

		# pdb.set_trace()

		

		