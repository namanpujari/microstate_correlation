'''
Validator script made for the puropse of checking whether the connection between
microstate_analysis.py/correl_stat.py and head3.lst/fort.38 (expected output) is 
legitimate.
'''

from __future__ import print_function	
# General imports...
from correl_stats import MicrostateCorrelation
from microstate_analysis import MicrostateAnalysis
from utils import function_timer
import numpy as np
import os
import sys

'''
1.  Verification function

	2 dictionaries.
	First should be update to head3.lst update in microstate_analysis
	so that the dictionary ALSO contains occupancies. head3.lst already has
	the conformer ids.

	update my own classs to take in fort.38, converts into dictionary where keys
	are conformer ids and the values are occupancies. This will be checked against
	the dictionary created in the previous step
	
	can return list of ids that do not match up, list of ids that match up (in terms of matching occupancies) 
	and list of missing ids
	
2.  Expand trajectory based on state count next to each record.
	
	E.g. if we have [1, 32 ... 43, 21, 31] then 2 as state count then the trajectory
	would be expanded into 
	
	[[1, 32 ... 43, 21, 31]
	 [1, 32 ... 32, 21, 31]]

step2 out check if different chains have different ids
'''

def convert_to_occupancies(count_dictionary. return_case=False):
	''' Follows the instruction shown in step 2 of the next_steps.txt file.
	Runs an instance of MicrstateCorrelation, calculated occupancies, then compares
	with the one displayed in fort.38 (actual occupancies)
	'''
	# count_dictionary is with residue
	# as key and list of conf ID and its instances as value. THESE VALUES WILL NOT BE ENTIRELY CORRECT
	# AS WE ARE NOT TAKING INTO ACCOUNT THE EXPANDED TRAJECTORY (WHERE WE CONSIDER STATE COUNT FOR 
	# EACH RECORD). Let's do it without changing it first
	self.occupancy_dict = {}
	for key in count_dictionary.keys():
		residue_conf_count_list = count_dictionary[key]
		occupancy_dict[key] = np.array([])
		for row in range(residue_conf_count_list.shape[0]):
			# We will use total_records to calculate conformer occupancy for now.
			# When the problem mentioned earlier is solved, we can transition to 
			# state counts
			conf_id = residue_conf_count_list[row][0]
			occupancy = residue_conf_count_list[row][1] / MicrostateCorrelation.total_records # total_record to change to state_counts
			############################################## ^ IT WOULD BE BETTER TO MAKE THIS  ####
			##############################################     FUNCTIONS INTO METHODS IN THE  ####
			##############################################       MICROSTATECORRELATION CLASS  ####
			occupancy_array = np.array([conf_id, occupancy])
			np.append(occupancy_dict[key], occupancy_array, axis = 0)
	# By now we have a new dictionary, with each key having a value in the form a list showing conformer id
	# and its corresponding occupancy.
	if(return_case = True):
		return occupancy_dict


def compare_with_fort38(fort38_location, theoretical_occ, experimental_counts):
	conformers_extra = []
	conformers_missing = []
	conformers_incorrect = []
	experiental_occ = convert_to_occupancies(experimental_counts)
	with open(fort38_location, 'r') as fort38:
		interesection = np.interesect1d(np.array(theoretical_occ.keys()), np.array(experimental_occ.keys()))
		union = np.union1d(np.array(theoretical_occ.keys()), np.array(experimental_occ.keys()))
		print(interesection)
		print(union)
		
		return conformers_incorrect, conformers_missing

if __name__ == '__main__':
	correlation_instance = MicrostateCorrelation('ms.dat', 'head3.lst')
	# Extracting all the relevant properties/values/methods from the class, to apply 
	# in this context
	conf_id_data = correlation_instance.conformer_data
	correlation_instance.conformer_count()
	conf_counts = correlation_instance.conf_count_dict
	print(type(conf_counts))
 

	# Now we can execute the functions above.
	compare_with_fort38('fort.38', conf_id_data, conf_counts)




