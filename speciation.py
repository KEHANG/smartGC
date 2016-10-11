import os
import json
from calibration import (
	load_fitted_calibration_factor_functions,
	read_calibration_species_constants
	)

def prepare_speciation_in_moles_per_total_mass(
										condition,
										sample,
										calibration_factor_function_dict=None
										):
	gc_measurement_path = os.path.join('data', 'measurement')
	gc_speciation_file = os.path.join(
								gc_measurement_path, 
								condition, 'gc_speciation',
								'{0}.txt'.format(sample)
								)
	gc_inner_standard_file = os.path.join(
								gc_measurement_path,
								condition, 'gc_inner_standard',
								'{0}.csv'.format(sample)
								)

	gc_speciation_data_dict = read_gc_speciation_file(gc_speciation_file)

	gc_inner_standard_data_dict = read_gc_inner_standard_file(gc_inner_standard_file)

	if calibration_factor_function_dict is None:
		calibration_factor_function_dict = load_fitted_calibration_factor_functions()

	speciation_dict_in_moles_per_total_mass = calculate_speciation_in_moles_per_total_mass(
												gc_speciation_data_dict,
												gc_inner_standard_data_dict,
												calibration_factor_function_dict
												)
	# save speciation data in mole per total liquid mass (moles/g)
	save_results_path = os.path.join(gc_measurement_path,
									condition, 
									'speciation_results'
									)
	if not os.path.exists(save_results_path):
		os.mkdir(save_results_path)

	speciation_in_moles_per_total_mass_save_file = os.path.join(
												save_results_path,
												'{0}.json'.format(sample)
												)


	with open(speciation_in_moles_per_total_mass_save_file, 'w') as write_out:
		json.dump(speciation_dict_in_moles_per_total_mass, write_out, indent=2)

	return speciation_dict_in_moles_per_total_mass

def read_gc_speciation_file(gc_speciation_file):

	gc_speciation_data_dict = {}
	with open(gc_speciation_file, 'r') as read_in:
		start_read = False
		
		for line in read_in:
			if start_read and line.strip() != '':
				fields = line.split('\t')
				peak_area = float(fields[4].strip())
				species = fields[7].strip()
				if species not in gc_speciation_data_dict:
					gc_speciation_data_dict[species] = peak_area
				else:
					gc_speciation_data_dict[species] += peak_area

			if line.startswith('Peak'):
				start_read = True

	return gc_speciation_data_dict

def read_gc_inner_standard_file(gc_inner_standard_file):

	gc_inner_standard_data_dict = {}
	with open(gc_inner_standard_file, 'r') as read_in:
		start_read = False
		
		for line in read_in:
			if start_read and line.strip() != '':
				fields = line.split(',')
				total_liquid_mass = float(fields[0].strip())
				inner_standard = fields[1].strip()
				inner_standard_mass = float(fields[2].strip())

				gc_inner_standard_data_dict['total_liquid_mass(g)'] = total_liquid_mass
				gc_inner_standard_data_dict['inner_standard'] = inner_standard
				gc_inner_standard_data_dict['inner_standard_mass(g)'] = inner_standard_mass

			if line.startswith('total_liquid_mass'):
				start_read = True

	return gc_inner_standard_data_dict

def calculate_speciation_in_moles_per_total_mass(
										gc_speciation_data_dict,
										gc_inner_standard_data_dict,
										calibration_factor_function_dict
										):

	speciation_dict_in_moles_per_total_mass = {}

	# get inner standard
	inner_standard = gc_inner_standard_data_dict['inner_standard']
	inner_standard_calibration_factor = calibration_factor_function_dict[inner_standard][0]
	inner_standard_peak_area = gc_speciation_data_dict[inner_standard]
	inner_standard_mass = gc_inner_standard_data_dict['inner_standard_mass(g)']
	inner_standard_MW = read_calibration_species_constants()[inner_standard]['MW']
	total_liquid_mass = gc_inner_standard_data_dict['total_liquid_mass(g)']


	for species in gc_speciation_data_dict:

		if species != inner_standard and (species in calibration_factor_function_dict):

			species_peak_area = gc_speciation_data_dict[species]
			species_calibration_factor = calibration_factor_function_dict[species][0]
			species_moles_per_total_mass = (species_peak_area/species_calibration_factor)\
									/(inner_standard_peak_area/inner_standard_calibration_factor)\
									*inner_standard_mass/inner_standard_MW/total_liquid_mass

			speciation_dict_in_moles_per_total_mass[species] = species_moles_per_total_mass

	return speciation_dict_in_moles_per_total_mass

if __name__ == '__main__':
	
	condition = 'pyrolysis_450C'
	samples = [
			'new_sample1_bf_instd',
			'new_sample1_aft_instd', 
			'new_sample1_aft_instd_repeat',
			'new_sample2_bf_instd',
			'new_sample2_aft_instd_repeat', 
			'new_sample3_bf_instd',
			'new_sample3_aft_instd', 
			'new_sample3_aft_instd_repeat'
			]
	for sample in samples:
		prepare_speciation_in_moles_per_total_mass(
										condition,
										sample
										)

