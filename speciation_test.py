import unittest
import os
from speciation import (
	read_gc_speciation_file,
	read_gc_inner_standard_file,
	calculate_speciation_in_moles_per_total_mass
	)
from calibration import (
	load_fitted_calibration_factor_functions,
	)

class test_speciation(unittest.TestCase):

	def setUp(self):

		self.calibration_factor_function_dict = load_fitted_calibration_factor_functions()

	def test_read_gc_speciation_file(self):

		gc_measurement_path = os.path.join('data', 'measurement')
		gc_speciation_file = os.path.join(
									gc_measurement_path, 
									'test_condition', 'gc_speciation',
									'sample0.txt'
									)

		gc_speciation_data_dict = read_gc_speciation_file(gc_speciation_file)

		self.assertEqual(len(gc_speciation_data_dict), 12)
		self.assertEqual(gc_speciation_data_dict['toluene'], 29494571)
		self.assertEqual(gc_speciation_data_dict['chlorothiophene'], 101958153)
		self.assertEqual(gc_speciation_data_dict['undecane'], 277875246)
		self.assertEqual(gc_speciation_data_dict['PDD'], 536041757)

	def test_read_gc_inner_standard_file(self):

		gc_measurement_path = os.path.join('data', 'measurement')
		gc_inner_standard_file = os.path.join(
										gc_measurement_path, 
										'test_condition', 'gc_inner_standard',
										'sample0.csv'
										)

		gc_inner_standard_data_dict = read_gc_inner_standard_file(gc_inner_standard_file)

		self.assertEqual(len(gc_inner_standard_data_dict), 3)
		self.assertAlmostEqual(gc_inner_standard_data_dict['total_liquid_mass(g)'], 0.2, 1)
		self.assertEqual(gc_inner_standard_data_dict['inner_standard'], 'chlorothiophene')
		self.assertAlmostEqual(gc_inner_standard_data_dict['inner_standard_mass(g)'], 0.02, 2)

	def test_calculate_speciation_in_moles_per_total_mass(self):

		gc_measurement_path = os.path.join('data', 'measurement')
		gc_speciation_file = os.path.join(
									gc_measurement_path, 
									'test_condition', 'gc_speciation',
									'sample0.txt'
									)

		gc_speciation_data_dict = read_gc_speciation_file(gc_speciation_file)

		gc_inner_standard_file = os.path.join(
										gc_measurement_path, 
										'test_condition', 'gc_inner_standard',
										'sample0.csv'
										)

		gc_inner_standard_data_dict = read_gc_inner_standard_file(
										gc_inner_standard_file
										)

		speciation_dict_in_moles_per_total_mass = calculate_speciation_in_moles_per_total_mass(
													gc_speciation_data_dict,
													gc_inner_standard_data_dict,
													self.calibration_factor_function_dict
													)
		self.assertIn('PDD', speciation_dict_in_moles_per_total_mass)
		self.assertAlmostEqual(speciation_dict_in_moles_per_total_mass['PDD']/(1e-2), .15, 1)
