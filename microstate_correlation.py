"""
Created on Fri Dec 29 16:06:39 2017

@author: Naman

Remarks: important class properties

Properties: 

residue_list (a python array containing all of the residues
in the ms.dat file that is analyzed)
trajectory (a total_records * n_res size ndarray
that stores the conformer of a residues by column, that are in a particular
microstate. This is the most important feature of ms.dat and the parser)
state_count (each record can have particular conformer configuration, and how many
times that microstate was considered. The latter part is stored in state_count, and
its size, by basic logic, is total_records)
energies (contains the energy at the given record, is a np.array of size total_records)
"""

from __future__ import print_function
from __future__ import division	
import os
import sys
import numpy as np
from scipy.stats import itemfreq
from microstate_analysis import MicrostateAnalysis
from utils import function_timer


class MicrostateCorrelation(object):
	def __init__(self, ms_dat_file, head3_lst_file, fort38_location):
		# Basic inheritance
		self.MA = MicrostateAnalysis(ms_dat_file, head3_lst_file, head3required = True)
		self.MA.generate_byte_indices()
		self.MA.parse_records()
		self.fort38_location = fort38_location
		self.fort38_by_id = {}
		self.fort38_by_name = {}
		#self.total_records = self.MA.total_records
		#self.state_counts = self.MA.state_counts
		#self.trajectory = self.MA.trajectory
		#self.res_list = self.MA.residue_list

	@function_timer
	def conformer_count(self, trajectory_as_file=False, return_case=False):
		if(trajectory_as_file == True):
			# Look for a trajectory data file, since we do not
			# need to run the parser again, in this case. 
			# read, convert to np.array
			# self.trajectory = newly created np.array
			sys.exit('Not yet implemented')
		# We can use an in-built, scipy.stats method known as itemfreq
		self.conf_count_dict = {} # Initiate python dictionary for counting the 
		# conformer counts and matching it with respective residues. This is 
		# a useful data structure in this case.
		self.PARSED_COUNT_BY_ID = {}
		for residue in range(0, self.MA.trajectory.shape[1], 1):
			res_conf_count = itemfreq(self.MA.trajectory[:, residue].ravel())
			self.conf_count_dict[self.MA.residue_list[residue]] = res_conf_count
			for row in range(res_conf_count.shape[0]):
				# We will use total_records to calculate conformer occupancy for now.
				# When the problem mentioned earlier is solved, we can transition to 
				# state counts
				conf_id = res_conf_count[row][0]
				self.PARSED_COUNT_BY_ID[conf_id] = res_conf_count[row][1]
		if(return_case == True): # Add functionality so it is useable outside
		# the scope of the class 
			return self.conf_count_dict

	@function_timer
	def get_occupancies(self, count_arr=None, return_case=False):
		''' 
		Follows the instruction shown in step 2 of the next_steps.txt file.
		Runs an instance of MicrstateCorrelation, calculated occupancies, then compares
		with the one displayed in fort.38 (actual occupancies)
		'''
		self.extract_fort38() # so that we can use the id->name and name->id functions
		if(count_arr != None):
			self.conf_count_dict = count_arr
		# count_dictionary is with residue
		# as key and list of conf ID and its instances as value. THESE VALUES WILL NOT BE ENTIRELY CORRECT
		# AS WE ARE NOT TAKING INTO ACCOUNT THE EXPANDED TRAJECTORY (WHERE WE CONSIDER STATE COUNT FOR 
		# EACH RECORD). Let's do it without changing it first
		self.occupancy_dict = {}
		self.PARSED_OCCUPANCY_BY_ID = {}
		for key in self.conf_count_dict.keys():
			residue_conf_count_list = self.conf_count_dict[key]
			self.occupancy_dict[key] = np.array([])
			for row in range(residue_conf_count_list.shape[0]):
				# We will use total_records to calculate conformer occupancy for now.
				# When the problem mentioned earlier is solved, we can transition to 
				# state counts
				conf_id = residue_conf_count_list[row][0]
				occupancy = residue_conf_count_list[row][1] / float(self.MA.total_records) # total_record to change to state_counts
				occupancy_array = np.array([conf_id, occupancy])
				np.append(self.occupancy_dict[key], occupancy_array, axis = 0)
				self.PARSED_OCCUPANCY_BY_ID[conf_id] = occupancy
		if(return_case == True):
			return self.occupancy_dict, self.PARSED_OCCUPANCY_BY_ID # the first is a dictionary keyed by residues, with each value containing
			# a n X 2 list of its conformer ids and corresponding energies. the second is a much bigger dictionary keyed by conformer id,
			# with each value being its energy.

	# -----------------|
	# UTILITY FUNCTIONS|
	# -----------------|

	def extract_fort38(self):
		with open(self.fort38_location, 'r') as fort38:
			lines = fort38.readlines()[1:] # starting from the second line
			for index, line in enumerate(lines):
				data_id = np.array([str(line.split()[0]), float(line.split()[1])])
				data_name = np.array([int(index+1), float(line.split()[1])])
				self.fort38_by_id[int(index+1)] = data_id
				self.fort38_by_name[str(line.split()[0])] = data_name
				# two dictionaries have been populated. id to name is keyed by id and 
				# has values in the form of [conformer name, energy] and name to id is keyed
				# by name, and has values in the form of [id number, energy]

	def retrieve_from_name(self, conf_name):
		return self.fort38_by_name[conf_name][0]

	def retrieve_from_id(self, conf_id):
		return self.fort38_by_id[conf_id][0]

	# ---------------------------|
	# PARSER VALIDATION FUNCTIONS|
	# ---------------------------|

	def fort38_compare(self):
		# Let's compare parsed occupancies with the fort.38 occupancies. 
		# There is no use in relying on head3lst at this point, as, according
		# to the file, there is no energy in any conformer.
		# conformers_extra = th - th INTERSECT exp
		# conformers_missing = exp - th INTERSECT exp
		# conformers_incorrect = th INTERSECT exp WRONG VALUES
		theoretical_conformers = self.fort38_by_id.keys() # Each value is 2 member array. The second
		# value is the occupancy, unlike the case for PARSED_OCCUPANCY_BY_ID, which simply contains
		# floating point values for each of its keys.
		experimental_conformers = self.PARSED_OCCUPANCY_BY_ID.keys()
		intersection = np.intersect1d(theoretical_conformers, experimental_conformers)
		conformers_missing = np.setdiff1d(theoretical_conformers, intersection)
		conformers_extra = np.setdiff1d(experimental_conformers, intersection)

		RESULTS = { "extra": conformers_extra, 
					"missing": [], # Reason for keeping this array empty will be explained later
					"incorrect": [],
					"correct": [] }

		# Now we wish to validate the conformers that are understood by the parser to at least exist
		# We can approach this problem simply by using a for loop...
		for candidate_conformer in intersection:
			# Let's not be too strict with our judgement. Since data like this is bound to have
			# some sort of error, we will allow a 5% threshold for validation.
			th_value = self.fort38_by_id[candidate_conformer][1]
			exp_value = self.PARSED_OCCUPANCY_BY_ID[candidate_conformer]
			if((np.abs(float(th_value) - float(exp_value)) / float(th_value)) <= 0.05):
				RESULTS["correct"] = np.append(RESULTS["correct"], candidate_conformer) #RESULTS["correct"].append(candidate_conformer)
			else: RESULTS["incorrect"] = np.append(RESULTS["incorrect"], candidate_conformer) #RESULTS["incorrect"].append(candidate_conformer)

		# The algorithm for scipy.itemfreq does not return any sign of conformers that have never
		# occured, which might have made up the majority of conformers_missing. Let's check...
		for candidate_conformer in conformers_missing:
			if(float(self.fort38_by_id[candidate_conformer][1]) == 0.0): 
				RESULTS["correct"] = np.append(RESULTS["correct"], candidate_conformer)
			else: RESULTS["missing"] = np.append(RESULTS["missing"], candidate_conformer) #RESULTS["missing"].append(candidate_conformer)

		return RESULTS

		
if __name__ == '__main__':
	# Execute the functions here
	correl = MicrostateCorrelation('ms.dat', 'head3.lst', 'fort.38')
	print('Counting...')
	correl.conformer_count()
	correl.get_occupancies()

	# # TEST
	# res = correl.fort38_compare()
	# print("There are "  + str(len(res["correct"])) + " correct conformers")
	# print("There are "  + str(len(res["incorrect"])) + " incorrect conformers")
	# print("There are "  + str(len(res["extra"])) + " extra conformers")
	# print("There are "  + str(len(res["missing"])) + " missing conformers")

	# print("Approximately " + str(int(np.ceil(10 * (len(res["correct"]) / 
	# 		(len(res["incorrect"]) + len(res["correct"])))))) + "0% correct conformers")   
	
	# TEST
	# common, extra, missing = correl.fort38_compare()
	# print('THE COMMON CONFORMERS, THAT ARE NOT IN FORT 38, ARE:')
	# print(common)
	# print('THE AMOUNT OF COMMON CONFORMERS IS ' + str(len(common)))
	# print('THE EXTRA CONFORMERS, THAT ARE NOT IN FORT 38, ARE:')
	# print(extra)
	# print('THE AMOUNT OF EXTRA CONFORMERS IS ' + str(len(extra)))
	# print('THE MISSING CONFORMERS, THAT ARE NOT IN PARSER CONFORMERS, ARE:')
	# print(missing)
	# print('THE AMOUNT OF MISSING CONFORMERS IS ' + str(len(missing)))

	# TEST
	# a = correl.retrieve_from_name("NTR01A0014_001") # retrieve from name returns an array in the form of [id, occupancy]
	# # this occupancy is obtained from fort.38
	# b = correl.retrieve_from_id(1) # retrieve from id returns a similar array, but it simply uses the conformer_data dictionary
	# # that you created in microstate_analysis.py. I looked up how the dictionary was made; you use head3lst to obtain the occupancies
	# # for each conformer name.

	# print("This is the id for conformer name NTR01A0014_001: \n\n" + str(a))
	# print("\n This is the name for conformer id 1: \n\n" + str(b))

	# # TEST
	# for i in np.sort(correl.PARSED_OCCUPANCY_BY_ID.keys()):
	# 	 print(str(i) + "    " + str(correl.PARSED_OCCUPANCY_BY_ID[i]))
	# 	 print(len(correl.PARSED_OCCUPANCY_BY_ID.keys()))
	# 	 break

	# # TEST
	# for i in np.sort(correl.PARSED_COUNT_BY_ID.keys()):
	# 	 print(str(i) + "    " + str(correl.PARSED_COUNT_BY_ID[i]))
	# 	 print(len(correl.PARSED_COUNT_BY_ID.keys()))
	# 	 break

	# # TEST 
	# print(correl.MA.trajectory[-3:])
	# some_random_numbers = [0, 1, 2, 32, 56, 76, 101, 132, 211, 242, 323, 390]
	
	# for i in some_random_numbers:
	# 	print(correl.MA.residue_list[i])
	# 	print(correl.conf_count_dict[correl.MA.residue_list[i]])

