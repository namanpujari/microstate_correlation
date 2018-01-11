"""
Moduule containing useful classes for the analysis of eqiulibrium distribution of 
conformations and protomers, stored as microstate data files. Comments on code segments
made by Naman Pujari, in the interest of better understanding the inherited code. Any 
comment that seems to narrate what a module does, or what a specific 
piece of code does (apart from documentation) is most probably a part of this effort.
THE CODE ITSELF IS BY KAMRAN HAIDER 
"""

from __future__ import print_function # Can be used to call print as a function, as seen
# in the following code: print("Hello World"). OUTPUT: Hello World. The python2.7 version
# of the very same is print "Hello World".
from __future__ import division # Is used to call the "future" implementation of division and 
# is TRUE division. Old '/' indicated floor(a/b), where the '/' in a/b is TRUE division.
from builtins import zip # Can allow us to create iterables out of tuples, or non-singular objects
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div # Floor division
import os
import struct # A useful python module for retreiving information in a binary file. It
# can also handle the task of doing the opposite: packing user-readable data into binary
# buffers as per the buffer protocol.

import numpy as np # Basic mathematical module used by data scientist far and wide

import networkx as nx # According to its github.io site, "NetworkX is a Python package
# for the creation, manipulation, and study of the structure, dynamics, and functions 
# of complex networks."

from utils import function_timer, print_progress_bar, is_grotthuss # utils is a python
# script also written by Kamran and includes certain utilities functions, that he calls
# in line 29

class MicrostateAnalysis(object):
    """A class representing microstate data and various analysis methods"""
    # The class requires 2 arguments, as seen from the __init__ method. Since none of the
    # arguments are optional, a head3lst file is MANDATORY along with the ms.dat file
    def __init__(self, ms_data_file, head3lst_file, head3required=False):
        self.ms_data_file = ms_data_file # A variable assigned to the ms.dat file with which
        # we are working. 
        self.byte_indices = None 
        self.total_microstates = 0 # Variable name is self-explanatory
        self.total_records = 0 # Variable name is self-explanatory.

        res_list = [] # Python array dedicated to the list of all residues read up to
        # this point.
        ########################### OPEN FILE ################################
        with open(self.ms_data_file, "rb") as md: # rb mode specializing in importing non rich-text
            # files (binary files, etc.)
            bytes_n_res = md.read(4) # Holds the first 4 characters, or in other words, the 
            # first four BYTES into a variable named bytes_n_res
            n_res = struct.unpack('i', bytes_n_res)[0] # Struct unpack converts byte-strings into 
            # readable data. The first variable 'i' indicates that the datatype we are seeking is in 
            # fact an integer. The second indicates the byte-sequence we want to obtain the data from.
            # This returns a tuple of (possibly) multiple values, with the first value always being 
            # the integer we desire. Hence we add a [0] to the end of the statement, giving us only 
            # the integer.
            for i in range(n_res):
                resname = str(md.read(8), 'utf-8') # Converts to string, an 8-byte sequence of characters
                res_list.append(resname) # Simply appends the residue name to the res_list
            self.n_res = n_res # Makes n_res a class property
        ########################## CLOSE FILE ################################
        self.residue_list = res_list # Make the res_list a class property with the name residue_list
        self.residue_hb_matrix = np.zeros((n_res, n_res), dtype=float) # Creates a square matrix based 
        # off of the amount of residues read by the file interpreter. Note that the size of the 
        # residue_list array SHOULD be the same as the value of n_res. Initiates the mentioned array
        # with 0.0 values (as they are floats).

        #with open(ms_gold_file, "r") as ms_gold:
        #    self.n_res = len([res.strip() for res in ms_gold.readlines() if len(res) != 0])
        
        if head3required == True:
            conf_data = {} # Initialized dictionary, conf_data. (Data for the conformers)
            with open(head3lst_file, "r") as h3: # Opened in "r" mode so that the readlines() and split()
            # functionalities are available for use. 
                for line in h3.readlines()[1:]: # Gather the contents of all lines into an array, starting
                # from the second line, to the end of the file. Iterate along the lines in that array.
                    data = line.split() # For each line, split each unique, non-white space entry, as a
                    # a member of a dedicated array for the current line. 
                    conf_id = int(data[0]) # The first unique value in the line is the conformer name/id.
                    # It is converted to an int and is passed to the conf_id variable.
                    conf_data[conf_id] = [data[1], data[3]]
                    # The data for that particular conformer is passed as an array containing the second 
                    # unique value and the 4th unique value in the line dedicated to that conformer id.
                    # This is stored in a dictionary fashion.
            self.conformer_data = conf_data 
            # The conformer_data is made a class property.
        elif head3required == False: 
            print('head3.lst not supplied and not required')
        
        # NOTE: It appears that Kamran does not really use the head3.lst file for any remaining
        # functions in this class. For testing purposes, we may simply remove the requirement of 
        # supplying a head3.lst file. 
        

    def generate_byte_indices(self, sample_frequency=100, filename=None):
        """
        Parameters
        ----------
        n_res : TYPE
            Description
        sample_frequency : int, optional
            Description
        filename : None, optional
            Description
        
        Returns
        -------
        rec_indices : list
            A list of the starting bytes for each record in the microstate data file
        """
        if filename is None:
            # The logic above is introduced so that this particular method can be available outside
            # the context of the class. That is, we will be able to use any other ms.dat file without
            # first initiated an instance of the class MicrostateAnalysis dedicated to it. 
            filename = self.ms_data_file
        

        start_byte = 4 + (8 * self.n_res) # The first 4 bytes (int32 datatype) and then we have a slew
        # of 8-byte length residue names. These are overlooked since they have been taken care of in
        # the __init__ function. 
        bytes_per_record = (self.n_res * 2) + 20
        # Here record is a "record of a microstate", which includes all the conformer ids and after
        # this sequence of conform ids, some values (including energy)
        file_size = os.path.getsize(filename) # getsize() is an inbuilt python functions that 
        # returns, in bytes, the size of the file we are dealing with. 

        # n_records = (file_size - start_byte) / bytes_per_record
        
        # Sample frequency is regarding the amount of checks of the records. This is done
        # knowing the fact that there is an equilibriation step before the mc sampling (that produces
        # the microstates)
        rec_indices = list(range(start_byte, file_size, sample_frequency * bytes_per_record))
        # Using an array, we smartly obtain the starting indices of each microstate record. The 
        # sample_frequency variable jumps across the file whilst always landing on the start of a 
        # record (insured by a step value of sample_frequency * bytes_per_record)
        self.total_records = len(rec_indices)
        # len() functions to get the amount of records obtained with the given sample_frequency
        # If you think about it, this is a WAY smarter method than just computing
        # the difference between file_size and start_byte, then dividing by the length of one record. 
        # Using the array method we have both the indices and the amount of records in the file. 
        # self.total_records is simply the number of records parsed with the preset 
        # sampling_frequency value.
        self.byte_indices = rec_indices
        # The record_indices are made into a class property with the name "byte_indices"
        # Hence self.byte_indices == rec_indices is True at all times after the usage of this
        # method.
        
        # RECORDS != MICROSTATES. A record may indicate a conf id configuration, and then how many
        # times that microstate has existed. So in some ways, RECORDS = some value * microstates
        
    @function_timer
    def parse_records(self):
        """
        
        """
        # This is where the real parsing happens.
        trajectory = np.zeros([self.total_records, self.n_res], dtype=int)
        state_counts = np.zeros([self.total_records], dtype=int)
        energies = np.zeros([self.total_records], dtype=float)
        #print trajectory
        progress_counter = 0 # Ties into the progress bar that Kamran created in utils.py
        print_progress_bar(progress_counter, self.total_records) # The function called outside the 
        # for loop to initiate with an empty bar. 
        with open(self.ms_data_file, "rb") as ms:
            for index, record in enumerate(self.byte_indices): # enumerate returns an iterable that
                # provides coders with a counter along with the values that they wish to iterate
                # across. Hence enumerate(self.byte_indices) returns for each entry in the list, 
                # its index as one variable and its value as another. 
                ms.seek(record) # seek() basically goes to some value in the file. Here we go
                # the location specified by the current record value (which is in bytes, and is 
                # the beginning of a microstate). 
                bytes_conf_ids = ms.read(2 * self.n_res) # Starting from the record index, the 
                # conformer id is exactly 2 * self.n_res bytes long. This is largely due to the fact
                # that each microstate contains a conformer of each residue. 
                bytes_energies_1 = ms.read(8) # Read in the energy corresponding to that particular 
                # microstate. 
                ms.seek(ms.tell() + 8) # Skip ahead 8 bytes, because we do not care about the 
                # information located at that position.
                energy = struct.unpack("d", bytes_energies_1)[0] # This is where we will convert the
                # binary data to a double datatype. Note, a double is a float with much higher 
                # precision and range. Energy binary is converted to a decimal. 
                bytes_state_count = ms.read(4) # The remaining 4 bytes are stored into some value
                # unknown to me right now, but possibly relating to a positive number telling us 
                # the amount of times a microstate has occcured. 
                trajectory[index, :] = np.asarray(struct.unpack(str(self.n_res) + "H", bytes_conf_ids))
                # The bytes containing conformer ids are unpacked into small unsigned short datatype
                # segments. These are then stored in the row denoted by index, which is defined by the
                # for-loop. For this to work it is required that the product of np.asarray(...) has
                # the same amount of columns as the trajector, which is the number of residues. Hence
                # my hypothesis at this time is that the RHS contains conformer ids for each residue.

                #print(struct.unpack(str(self.n_res) + "H", bytes_conf_ids)[-2:])

                state_count = struct.unpack("i", bytes_state_count)[0] # Converts the binary data
                # of the state count to an integer.
                self.total_microstates += state_count # Class property total_microstate is increased
                # by the state count of the particular microstate the for-loop is currently on. 
                state_counts[index] += state_count # The value for microstate occurence in the current
                # record is recorded in the index we are on, for the list state_counts. 
                energies[index] += energy # Similiarly, the energy corresponding to that microstate
                # is also recorded.
                progress_counter += 1 # Update that a step has been completed. 
                print_progress_bar(progress_counter, self.total_records) # Update the progress bar
                # with a higher percentage complete value. 
        self.trajectory = trajectory # Make trajectory a class property.
        self.state_counts = state_counts # Make state_counts a class property.
        self.energies = energies # Make energies a class property.

    def write_hbtxt(self, prefix= ""):
        THRESHOLD_TO_WRITE = 0.001
        hbtxt_file = prefix + "hb.txt"
        print("Writing data to %s" % hbtxt_file)
        with open(hbtxt_file, "w") as t:
            for i in range(len(self.residue_list)):
                for j in range(len(self.residue_list)):
                    if old_div(self.residue_hb_matrix[i, j], self.total_microstates) >= THRESHOLD_TO_WRITE:
                        t.write(
                            "%s\t%s\t%.3f\n" %
                            (self.residue_list[i],
                             self.residue_list[j],
                             old_div(self.residue_hb_matrix[i, j], self.total_microstates)))


class ConformerHbonds(object):
    """A class to represent conformer hydrogen bond matrix and associated methods.
    """
    def __init__(self, hbdat_file):
        self.hbdat_file = hbdat_file

    @function_timer
    def generate_hbmatrix_from_hbdat(self, hbdat_file=None):
        """
        Parameters
        ----------
        hbdat_file : None, optional
            Description
        """
        if hbdat_file is None:
            hbdat_file = self.hbdat_file

        with open(hbdat_file, "rb") as hbhandle:
            bytes_n_conf = hbhandle.read(4)
            n_conf = struct.unpack('i', bytes_n_conf)[0]
            #print "Total conformers: ", n_conf
            bytes_h_matrix = hbhandle.read(n_conf * n_conf)
            hb_array = struct.unpack(str(n_conf * n_conf) + "B", bytes_h_matrix)
            self.hb_matrix = np.array(hb_array).reshape(n_conf, n_conf)
            
    def generate_microstate_hb_network(self, conformer_ids=None, labels=None):
        """
        Parameters
        ----------
        conformer_ids : TYPE
            Description
        """
        if conformer_ids is None:
            G = nx.to_networkx_graph(self.hb_matrix)
            return G
        else:
            G = nx.to_networkx_graph(self.hb_matrix[conformer_ids, :][:, conformer_ids])
            if labels is None:
                nx.relabel_nodes(G, dict(list(zip(G.nodes(), conformer_ids))), copy=False)
            else:
                assert len(labels) == len(conformer_ids), "Please make sure number of nodes and labels are equal."
                nx.relabel_nodes(G, dict(list(zip(G.nodes(), labels))), copy=False)
            return G



    def check_conformer_hbonds_pairwise(self, conformer_a, conformer_b):
        """
        Parameters
        ----------
        conformer_a : TYPE
            Description
        conformer_b : TYPE
            Description
        """
        pass

    def all_hbonds_by_conformer(conformer):
        pass

        
class HydrogenBondedPath(object):
    """
    """
    def __init__(self, start_residue, end_residue):
        """
        Parameters
        ----------
        start_residue : TYPE
            Description
        end_residue : TYPE
            Description
        """
        self.start_residue = start_residue
        self.end_residue = end_residue
        self.states_analyzed = 0
        self.path_list = []
        self.path_lengths = []
        self.record_indices = []

    def search_for_path(self, hb_network, record_index, start, end):
        """
        Parameters
        ----------
        hb_network : TYPE
            Description
        start : TYPE
            Description
        end : TYPE
            Description
        """
        try:
            path = nx.shortest_path(hb_network, source=start, target=end)
            self.path_list.append(path)
            self.path_lengths.append(len(path))
            self.record_indices.append(record_index)

        except nx.NetworkXNoPath as nopath:
            self.path_list.append(None)
            self.path_lengths.append(np.nan)
            self.record_indices.append(np.nan)

        finally:
            self.states_analyzed += 1

    def filter_paths(self, conformer_data):
        """
        Parameters
        ----------
        conformer_data : TYPE
            Description
        
        """
        path_lengths_array = np.asarray(self.path_lengths)
        path_list_array = np.asarray(self.path_list)
        record_indices_array = np.asarray(self.record_indices)
        filtered_paths = []
        filtered_path_records = []

        min_length = np.nanmin(path_lengths_array)
        shortest_paths = path_list_array[np.where(path_lengths_array == min_length)]
        shortest_path_records = record_indices_array[np.where(path_lengths_array == min_length)]
        for index, path in enumerate(shortest_paths):
            non_grotthuss = [conformer_data[conf][0] for conf in path[1:-1] if not is_grotthuss(conformer_data[conf][0][0:3])]
            #print non_grotthuss
            if len(non_grotthuss) == 0:
                filtered_paths.append(list(path))
                filtered_path_records.append(int(shortest_path_records[index]))
        return filtered_paths, filtered_path_records
    

