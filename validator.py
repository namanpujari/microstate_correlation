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

def do_comparision_step2(fort38_location):
	''' Follows the instruction shown in step 2 of the next_steps.txt file.
	Runs an instance of MicrstateCorrelation, calculated occupancies, then compares
	with the one displayed in fort.38 (actual occupancies)
	'''
	correlation_instance = MicrostateCorrelation('ms.dat', 'head3.lst')
	conf_counts = MicrostateCorrelation.conformer_count() # Returns dictionary with residue
	# as key and list of conf ID and its instances as value. THESE VALUES WILL NOT BE ENTIRELY CORRECT
	# AS WE ARE NOT TAKING INTO ACCOUNT THE EXPANDED TRAJECTORY (WHERE WE CONSIDER STATE COUNT FOR 
	# EACH RECORD). Let's do it without changing it first
	occupancy_dict = {}
	for key in conf_counts.keys():
		residue_conf_count_list = conf_counts[key]
		occupancy_dict[key] = np.array([])
		for row in residue_conf_count_list.shape[0]:
			# We will use total_records to calculte conformer occupancy for now.
			# When the problem mentioned earlier is solved, we can transition to 
			# state counts
			conf_id = residue_conf_count_list[row][0]
			occupancy = residue_conf_count_list[row][1] / MicrostateCorrelation.total_records
			occupancy_array = np.array([conf_id, occupancy])
			np.append(occupancy_dict[key], occupancy_array, axis = 0)
	# By now we have a new dictionary, with each key having a value in the form a list showing conformer id
	# and its corresponding occupancy.
	with open(fort38_location, 'r') as fort38:
		occupancy_dict_from_fort38 = {}
		for line in fort38_location.readlines()[1:]:
			# To be continued...
	return conformers_incorrect, conformers_missing



