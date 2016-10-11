import os
import statsmodels.api as sm
import json

def read_calibration_file(calibration_file):

	injection_volumes = []
	peak_areas = []
	with open(calibration_file, 'r') as read_in:
		start_read = False
		
		for line in read_in:
			if start_read and line.strip() != '':
				fields = line.split('\t')
				injection_volume = float(fields[0].strip())
				peak_area = float(fields[4].strip())
				injection_volumes.append(injection_volume)
				peak_areas.append(peak_area)

			if line.startswith('volume'):
				start_read = True

	return injection_volumes, peak_areas

def read_calibration_species_constants():

	calibration_species_constants_file = os.path.join(
												'data', 
												'calibration', 
												'calibration_species_constants.csv'
												)

	calibration_species_constants_dict = {}
	with open(calibration_species_constants_file, 'r') as read_in:
		start_read = False
		
		for line in read_in:
			if start_read and line.strip() != '':
				fields = line.split(',')
				species = fields[0].strip()
				density = float(fields[1].strip())
				MW = float(fields[2].strip())

				calibration_species_constants_dict[species] = {}
				calibration_species_constants_dict[species]['density'] = density
				calibration_species_constants_dict[species]['MW'] = MW

			if line.startswith('species'):
				start_read = True

	return calibration_species_constants_dict

def read_calibration_data(calibration_species_list=None, calibration_path=None):
	"""
	This method reads GC calibration files and output 
	a dictionary named `calibration_data_dict` with 
	key: `species` and value: tuple `(injection_moles, peak_areas)`
	"""

	if calibration_path is None:
		calibration_path = os.path.join('data', 'calibration')

	if calibration_species_list is None:
		calibration_species_list = []

		for f in os.listdir(calibration_path):
			if os.path.isfile(os.path.join(calibration_path, f)) and f.endswith('.txt'):
				calibration_species_list.append(f.split('.txt')[0])

	# load calibration species constants such as density, MW
	calibration_species_constants_dict = read_calibration_species_constants()

	calibration_data_dict = {}
	for species in calibration_species_list:
		calibration_file = os.path.join(calibration_path, species+'.txt')

		injection_volumes, peak_areas = read_calibration_file(calibration_file)

		# convert volume/uL to moles
		density = calibration_species_constants_dict[species]['density'] # unit: g/cm3
		MW = calibration_species_constants_dict[species]['MW'] # unit: g/mol
		injection_moles = [volume*1e-3*density/MW for volume in injection_volumes]

		calibration_data_dict[species] = (injection_moles, peak_areas)

	return calibration_data_dict

def get_calibration_factor_function_dict_by_linear_regression(calibration_data_dict, zero_intercept=True):
	"""
	This method fits the calibration experimental data into the form 
	`peak_area = a*injection_moles` or `peak_area = a*injection_moles + b`
	with the choice of `zero_intercept`, which is True by default.
	
	Here I chose to only use the first 3 data points is due to the fact that
	these 3 points covers the concentration range of species in experiment. But 
	user can play with it with his own situation.

	"""

	calibration_factor_function_dict = {}
	for species in calibration_data_dict:
		injection_moles, peak_areas = calibration_data_dict[species]
		injection_moles.insert(0, 0.0)
		peak_areas.insert(0, 0.0)

		y_to_fit = peak_areas[:]
		x_to_fit = injection_moles[:]

		if not zero_intercept:
			x_to_fit = sm.add_constant(x_to_fit)
		
		model = sm.OLS(y_to_fit, x_to_fit)
		calibration_factor_function_dict[species] = list(model.fit().params)

	return calibration_factor_function_dict


def prepare_calibration_factor_functions(calibration_species_list=None, calibration_path=None):

	# read calibration data
	calibration_data_dict = read_calibration_data(calibration_species_list, calibration_path)

	# linear regression
	calibration_factor_function_dict = get_calibration_factor_function_dict_by_linear_regression(
											calibration_data_dict)	

	# save calibration factor functions
	calibration_factor_function_save_file = os.path.join(
												'data', 
												'calibration', 
												'fitted_calibration_factor_functions.json')

	with open(calibration_factor_function_save_file, 'w') as write_out:
		json.dump(calibration_factor_function_dict, write_out, indent=2)

	# return calibration factor function dict
	return calibration_factor_function_dict

def load_fitted_calibration_factor_functions():

	calibration_factor_function_dict = {}

	fitted_calibration_factor_functions_path = os.path.join(
												'data', 
												'calibration', 
												'fitted_calibration_factor_functions.json')
	with open(fitted_calibration_factor_functions_path, 'r') as read_in:
		calibration_factor_function_dict = json.load(read_in)

	return calibration_factor_function_dict

if __name__ == '__main__':
	prepare_calibration_factor_functions()

