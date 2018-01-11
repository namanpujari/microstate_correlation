"""
Created on Fri Dec 29 16:06:39 2017

@author: Naman

Remarks: important class properties

Properties: 

residue_list (a python array containing all of the residues
in the ms.dat file that is analyzed)
trajectory (a n_records * n_res size ndarray
that stores the conformer of a residues by column, that are in a particular
microstate. This is the most important feature of ms.dat and the parser)
state_count (each record can have particular conformer configuration, and how many
times that microstate was considered. The latter part is stored in state_count, and
its size, by basic logic, is n_records)
energies (contains the energy at the given record, is a np.array of size n_records)
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
		self.MA = MicrostateAnalysis(ms_dat_file, head3_lst_file)
		self.MA.generate_byte_indices()
		self.MA.parse_records()
		self.trajectory = self.MA.trajectory
		self.res_list = self.MA.residue_list
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
		conf_count_dict = {} # Initiate python dictionary for counting the 
		# conformer counts and matching it with respective residues. This is 
		# a useful data structure in this case.
		for residue in range(0, self.trajectory.shape[1], 1):
			res_conf_count = itemfreq(self.trajectory[:, residue].ravel())
			conf_count_dict[self.res_list[residue]] = res_conf_count
		self.conf_counts = conf_count_dict # Transfer into a class property
		if(return_case == True): # Add functionality so it is useable outside
		# the scope of the class
			return conf_count_dict
		
if __name__ == '__main__':
	# Execute the functions here
	correl = MicrostateCorrelation('ms.dat', 'dummy_head3.txt')
	print('Counting...')
	correl.conformer_count()
	
	# TEST 
	some_random_numbers = [1, 2, 32, 56, 101, 211, 323]
	for i in some_random_numbers:
		print(correl.conf_counts.keys()[i])
		print(correl.conf_counts[correl.res_list[i]])
		print()
