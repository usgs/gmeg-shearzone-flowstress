from fugacity import *
import unittest 
import pdb

class Test_fugacity_calculator(unittest.TestCase):

	def setUp(self):
		self.fugacity_calculation = FugacityCalculator(temperature, pressure)

	def test_if_fugacity_calculator_initializes(self):
		self.assertEqual(self.fugacity_calculation[0][0], 723.15)

