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

	def setUp(self):
		self.fs_calculation = FlowStressCalculator(temperature, pressure, grain_size)

	def test_if_flow_stress_calculator_initializes(self):
		self.assertEqual(self.fs_calculation.temperature[1], 773.15)#723.15
		self.assertEqual(self.fs_calculation.pressure, 400000000)
		self.assertEqual(self.fs_calculation.grain_size, [10])
		self.assertEqual(self.fs_calculation.fugacity,[])
		self.assertEqual(self.fs_calculation.differential_stress,[])
		self.assertEqual(self.fs_calculation.strain_rate,[])

	def test_if_fugacity_calculator_works(self):
		self.fs_calculation.calculate_fugacity()
		self.assertEqual(self.fs_calculation.fugacity,[116.71651725233453, 149.72095732050246])

	def test_if_differential_stress_is_calculated(self):
		self.ds = self.fs_calculation.calculate_differential_stress() 
		self.assertEqual(self.ds, [78.75972939002484])
		
	def test_if_strain_rate_is_calculated(self):

		self.sr = self.fs_calculation.calculate_strain_rate(flow_law='H01')
	
		self.assertEqual(self.sr, [5.0181530333149198e-12])

	def test_if_slip_rate_is_calculated(self):

		self.slip = self.fs_calculation.calculate_slip_rate([5.01815303331492e-12], width)

		self.assertEqual(self.slip, [4.7475742217585797])