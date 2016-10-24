import os
from calibration import read_calibration_species_constants
import json

def save_exptl_data_to_chemkin_simulation_format(condition_before, condition_after):
	"""
	This method tries to gether speciation results before and after experiment
	to generate files for downstream chemkin simulation:
	
	input: files containing mol/g data for all samples under both 
	condition_before and condition_after.
	
	output: exptl_data files containing mol/g data before and after with time. Each file is 
	for one single sample. Format would be 
		{ "Time": [],
		"PDD": [],
		"toluene": [],
		...
		}

	"""

	# get all the files needed
	speciation_results_condition_before_path = os.path.join('data', 
														'measurement', 
														condition_before,
														'speciation_results'
														)

	speciation_results_condition_after_path = os.path.join('data', 
														'measurement', 
														condition_after,
														'speciation_results'
														)
	speciation_before_dict = {} # key: sample, value: [speciation1, speciation2]
	for f in os.listdir(speciation_results_condition_before_path):
		if os.path.isfile(os.path.join(speciation_results_condition_before_path, f)) \
			and f.endswith('.json'):
			sample = f.split('_bf')[0]
			with open(os.path.join(speciation_results_condition_before_path, f), 'r') as read_in:
				speciation = json.load(read_in)
				if sample not in speciation_before_dict:
					speciation_before_dict[sample] = [speciation]
				else:
					speciation_before_dict[sample].append(speciation)

	speciation_after_dict = {} # key: sample, value: [speciation1, speciation2]
	for f in os.listdir(speciation_results_condition_after_path):
		if os.path.isfile(os.path.join(speciation_results_condition_after_path, f)) \
			and f.endswith('.json') and ('_aft' in f):
			sample = f.split('_aft')[0]
			with open(os.path.join(speciation_results_condition_after_path, f), 'r') as read_in:
				speciation = json.load(read_in)
				if sample not in speciation_after_dict:
					speciation_after_dict[sample] = [speciation]
				else:
					speciation_after_dict[sample].append(speciation)

	# get the detailed condition specification of the experiment
	end_time, temperature, pressure = read_condition_details(condition_after)

	# construct exptl_data dictionary
	exptl_data_dict = {} # key: sample, value: { "Time": [..,..], "PDD": [..,..],...}
	for sample, speciations_after in speciation_after_dict.iteritems():
		exptl_data = {}
		for speciation in speciations_after:
			for species, moles_per_total_mass in speciation.iteritems():
				if species not in exptl_data:
					exptl_data[species] = [moles_per_total_mass]
				else:
					exptl_data[species].append(moles_per_total_mass)

		for species in exptl_data:
			exptl_data[species] = sum(exptl_data[species])/len(exptl_data[species])

		exptl_data_dict[sample] = exptl_data

	for sample, exptl_data in exptl_data_dict.iteritems():
		speciations_before = speciation_before_dict[sample]

		for species in exptl_data:
			moles_per_total_mass_list = []
			for speciation in speciations_before:
				try:
					moles_per_total_mass_list.append(speciation[species])
				except:
					moles_per_total_mass_list.append(0)

			exptl_data[species] = [sum(moles_per_total_mass_list)/len(moles_per_total_mass_list),
									exptl_data[species]
									]
		# normalize initial mol/g
		exptl_data = normalize_initial_moles_per_total_mass(exptl_data)


		#
		exptl_data['Time'] = [0, end_time] # unit: hour

	# save into json file, one per each sample
	save_dir = os.path.join('data', 'measurement', condition_after, 'exptl_data_for_simulation')
	if not os.path.exists(save_dir):
		os.mkdir(save_dir)
	for sample, exptl_data in exptl_data_dict.iteritems():
		save_filename = '{0}_{1}.json'.format(sample, condition_after.split('_')[1])
		with open(os.path.join(save_dir, save_filename), 'w') as write_out:
			json.dump(exptl_data, write_out, indent=2)

def normalize_initial_moles_per_total_mass(exptl_data):

	mass = 0
	calibration_species_constants_dict = read_calibration_species_constants()
	for species, moles_per_total_masses in exptl_data.iteritems():
		initial_moles_per_total_mass = moles_per_total_masses[0]
		if initial_moles_per_total_mass > 0:
			species_MW = calibration_species_constants_dict[species]['MW']
			mass += species_MW * initial_moles_per_total_mass

	for species, moles_per_total_masses in exptl_data.iteritems():
		moles_per_total_masses[0] = moles_per_total_masses[0]/mass

	return exptl_data

def read_condition_details(condition):

	# get condition.csv file
	condition_file = os.path.join('data', 
								'measurement', 
								condition,
								'condition.csv'
								)

	# parse condition.csv file
	with open(condition_file, 'r') as read_in:
		start_read = False
		
		for line in read_in:
			if start_read and line.strip() != '':
				fields = line.split(',')
				end_time = float(fields[0].strip())
				temperature = float(fields[1].strip())
				pressure = float(fields[2].strip())
				break

			if line.startswith('Time'):
				start_read = True

	return end_time, temperature, pressure

if __name__ == '__main__':

	save_exptl_data_to_chemkin_simulation_format('before_pyrolysis', 'pyrolysis_450C')