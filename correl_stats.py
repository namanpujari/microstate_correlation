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
import os
import sys
import numpy as np
from scipy.stats import itemfreq
from microstate_analysis import MicrostateAnalysis
from utils import function_timer

class MicrostateCorrelation(object):
	def __init__(self, ms_dat_file, head3_lst_file):
		# Basic inheritance
		self.MA = MicrostateAnalysis(ms_dat_file, head3_lst_file, head3required = True)
		self.MA.generate_byte_indices()
		self.MA.parse_records()
		self.total_records = self.MA.total_records
		self.state_counts = self.MA.state_counts
		self.trajectory = self.MA.trajectory
		self.res_list = self.MA.residue_list
		self.conformer_data = self.MA.conformer_data
		# Handed over the trajectory to THIS particular class
		
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
		for residue in range(0, self.trajectory.shape[1], 1):
			res_conf_count = itemfreq(self.trajectory[:, residue].ravel())
			self.conf_count_dict[self.res_list[residue]] = res_conf_count
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
		if(count_arr != None):
			self.conf_count_dict = count_arr
		# count_dictionary is with residue
		# as key and list of conf ID and its instances as value. THESE VALUES WILL NOT BE ENTIRELY CORRECT
		# AS WE ARE NOT TAKING INTO ACCOUNT THE EXPANDED TRAJECTORY (WHERE WE CONSIDER STATE COUNT FOR 
		# EACH RECORD). Let's do it without changing it first
		self.occupancy_dict = {}
		self.id_to_occ_dict = {}
		for key in self.conf_count_dict.keys():
			residue_conf_count_list = self.conf_count_dict[key]
			self.occupancy_dict[key] = np.array([])
			for row in range(residue_conf_count_list.shape[0]):
				# We will use total_records to calculate conformer occupancy for now.
				# When the problem mentioned earlier is solved, we can transition to 
				# state counts
				conf_id = residue_conf_count_list[row][0]
				occupancy = residue_conf_count_list[row][1] / float(self.total_records) # total_record to change to state_counts
				occupancy_array = np.array([conf_id, occupancy])
				np.append(self.occupancy_dict[key], occupancy_array, axis = 0)
				self.id_to_occ_dict[conf_id] = occupancy
		# By now we have a new dictionary, with each key having a value in the form a list showing conformer id
		# and its corresponding occupancy.
		if(return_case == True):
			return self.occupancy_dict
		return self.id_to_occ_dict

if __name__ == '__main__':
	# Execute the functions here
	correl = MicrostateCorrelation('ms.dat', 'head3.lst')
	print('Counting...')
	correl.conformer_count()
	correl.get_occupancies()

	for i in np.sort(correl.id_to_occ_dict.keys()):
		print(str(i) + "    " + str(correl.id_to_occ_dict[i]))


	# TEST 
	# print(correl.trajectory[-3:])
	# some_random_numbers = [1, 2, 32, 56, 76, 101, 132, 211, 242, 323, 390]
	
	# for i in some_random_numbers:
	# 	print(correl.conf_count_dict.keys()[i])
	# 	print(correl.conf_count_dict[correl.res_list[i]])

