from calibration import (
		read_calibration_data, 
		prepare_calibration_factor_functions
		)

from speciation import prepare_speciation_in_moles_per_total_mass

import unittest

class test_smartGC(unittest.TestCase):

	def test_smartGC_workflow(self):

		# A want to use GC data to get calibration factor functions 
		# for each species he measured
		calibration_factor_function_dict = prepare_calibration_factor_functions()

		self.assertIn('PDD', calibration_factor_function_dict)
		self.assertEqual(len(calibration_factor_function_dict['PDD']), 1)

		# next step: A wants to get mole-mass concentrations for the major species
		# (undecane, toluene, PDD, chlorothiophene) by feeding experimental 
		# GC speciation and GC inner standard files 
		# for sample0 under 450C pyrolysis condition

		condition = 'test_condition'
		sample = 'sample0'

		speciation_dict_in_moles_per_total_mass = prepare_speciation_in_moles_per_total_mass(
													condition,
													sample,
													calibration_factor_function_dict
													)

		self.assertIn('PDD', speciation_dict_in_moles_per_total_mass)
		self.assertAlmostEqual(speciation_dict_in_moles_per_total_mass['PDD']/(1e-2), .15, 1)
