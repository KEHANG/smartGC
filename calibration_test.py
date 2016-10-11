from calibration import (
					read_calibration_data, 
					read_calibration_file,
					read_calibration_species_constants,
					get_calibration_factor_function_dict_by_linear_regression
					)
import unittest
import os

class test_calibration(unittest.TestCase):

	def test_read_calibration_data(self):

		calibration_species_list = ['PDD', 'undecane', 'toluene', 'chlorothiophene']

		calibration_data_dict = read_calibration_data(calibration_species_list)

		self.assertEqual(len(calibration_data_dict), 4)
		self.assertEqual(len(calibration_data_dict['PDD']), 2)
		self.assertEqual(len(calibration_data_dict['toluene'][0]), 4)
		self.assertEqual(calibration_data_dict['chlorothiophene'][1][1], 261284787)

	def test_read_calibration_file(self):

		calibration_file = os.path.join('data', 'calibration', 'PDD.txt')
		injection_volumes, peak_areas = read_calibration_file(calibration_file)

		self.assertEqual(len(injection_volumes), 4)
		self.assertEqual(len(peak_areas), 4)
		self.assertEqual(peak_areas[0], 131261886)

	def test_read_calibration_species_constants(self):

		calibration_species_constants_dict = read_calibration_species_constants()

		self.assertEqual(len(calibration_species_constants_dict['PDD']), 2)
		self.assertEqual(calibration_species_constants_dict['PDD']['density'], .856)
		self.assertEqual(calibration_species_constants_dict['PDD']['MW'], 246.43)


	def test_get_calibration_factor_function_dict_by_linear_regression(self):

		calibration_species_list = ['PDD']

		calibration_data_dict = read_calibration_data(calibration_species_list)

		calibration_factor_function_dict = get_calibration_factor_function_dict_by_linear_regression(
											calibration_data_dict)

		self.assertEqual(len(calibration_factor_function_dict), 1)
		self.assertAlmostEqual(calibration_factor_function_dict['PDD'][0]/1e15, 0.32, 1)
