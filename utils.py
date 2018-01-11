"""Some utility functions for microstate hydrogen bond analysis.
"""
from __future__ import print_function

# Globals
import sys
import time
from functools import wraps

GROTTHUSS = ["GLU", "ASP", "SER", "THR", "HIS", "TYR", "ASN", "GLN", "HOH", "PAA", "PDD"]
# Cofactors in the heme group

# The function below returns a boolean value based on whether a given residue is within
# the GROTTHUSS list defined above. It is very easy to see how this works. The 'residue in
# GROTTHUS' statement basically checks if the arg of the function is within the list and returns
# a True/False value. 
def is_grotthuss(residue):
    return residue in GROTTHUSS

# Functions
def function_timer(function):
    # wraps is a wrapper decorator provided by functools that basically makes the wrapped
    # function the wrapper. In the interest of clarity, the @wraps(function) statement basically
    # does this: once function_timer is used as a decorator, like @function_timer before the declaration
    # of a function, calling function() and executing print(function.__name__) will actually
    # return the name of the function, not the name of the wrapper in the function_timer decorator, which
    # is 'function_timer'.
    @wraps(function)
    # Wrapper with the name 'function_timer'
    def function_timer(*args, **kwargs):
        # Use the python 'time' module to obtain the time right before the function
        # is executed. 
        t0 = time.time()
        # Execute function and provide result to a variable that is named 'result'
        result = function(*args, **kwargs)
        # Obtain the time right after the execution of the function. This can be used 
        # to calculate the total time
        t1 = time.time()
        # Prints the name of the function executed, and the total time taken
        print("Total time running %s: %2.2f seconds" %
              (function.__name__, t1 - t0))
        return result
    # Returns an unexecuted instance of the function as an object. This is callable
    return function_timer

# Can be repeatedly called to print the progress of a routine. It prints to the terminal repeatedly
# and removes the previous bar from the terminal in the same process, and makes an aesthetic loading
# screen look.
def print_progress_bar (count, total):
    """
    Create and update progress bar during a loop.
    
    Parameters
    ----------
    iteration : int
        The number of current iteration, used to calculate current progress. 
    total : int
        Total number of iterations
    
    Notes
    -----
        Based on:
        http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console
    """
    
    # length of the bar is defined to be 60, and is essentially just the PHYSICAL 
    # length of the loading bar. 
    bar_len = 60
    # Length of the blackened (filled in) segment of the bar is mathematically calculated
    # to be count/total * 60, so that the cap is 60, the physical length of the bar. 
    filled_len = int(round(bar_len * count / float(total)))
    # The percent progress rounded to 1 decimal place
    percents = round(100.0 * count / float(total), 1)
    # \u2588 is a a filled in black bar (may also be white depending on terminal color)
    # The line below creates the bar, and used * to print a certain str multiple times. For example
    # "naman" * 3 is >> namannamannaman
    bar = "o" * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('Progress |%s| %s%3s Done\r' % (bar, percents, '%'))
    sys.stdout.flush() 
    if count == total: 
        print()


# Depending on the substate given, generates clusters of residues.
def generate_cluster_dictionary(sub_state): 
    """
    Generate a dictionary containing key clusters and corresponding residues
    for a given sub_state.  
    Parameters
    ----------
    sub_state : string
        Specifify
    """
    # Lists of residues in each cluster
    clust_dict = {  
    "BNC" : ['CUBK9500', 'HE3K9300'],
    "PLS" : ['PDDK9102', 'PAAK9301', 'PDDK9302', 'TYRA0175',
            'ASPA0407', 'HISA0411', 'ASPB0229', 'ASPA0412', 'THRA0337', 'GLUB0254',
           'SERB0218', 'PAAK9101', 'ARGA0052', 'TYRA0414', 'SERA0497', 'SERA0498'],
    "EXT" : ['HISA0093', 'SERA0156', 'THRA0187', 'SERA0168', 'THRA0100',
            'TYRB0262', 'SERA0186', 'GLUA0182', 'ASPA0188', 'ASPA0485', 'GLUA0488'],
    "GLU" : ["GLUA0286"],
    "EXT_EXP" : ['HISA0093', 'GLUA0182', 'ASPA0188', 'ASPA0485', 'GLUA0488'],
    "EXT_INT" : ['SERA0168', 'THRA0100']
    }
    # Only works if substate defined is either f1, f2 or f4. 
    if sub_state not in ["f1", "f2", "f4"]:
        sys.exit("Only f1, f2 and f4 substates are currently supported.")
    # If substate is valid, and is either f1 or f4, then simply returns the dictionary above.
    if sub_state == "f1" or sub_state == "f4":
        return clust_dict
    # However, if the substate is f3, then it will return a modified version of the dictionary
    # above. 
    else:
        # using the .keys() method for dictionaries, converts the keys of the cluster dictionary 
        # into a list. Then iterates across those keys, accessing the data in each residue 
        # cluster.
        for k in list(clust_dict.keys()):
            # Initiate updated_cluster_residue list
            updated_cluster_residues = []
            # res is a variable that holds each item in the value of each key in the dictonary.
            # Since each value is a list, then res is simply each, individual item in the list
            for res in clust_dict[k]:
                # For the given cluster, if the given reside has an "A" as its fourth character,
                # then obtain the residue number and subtract 11. Append that modified residue
                # number in front of "04d", and then append that segment after the first three 
                # characters in the original residue name. Add that value to the updated_cluster_residue 
                # list.
                if res[3] == "A":
                    resnum = int(res[4:])
                    modified_resnum = resnum - 11
                    str_resnum = "%04d" % (modified_resnum)
                    modified_res = res[0:4] + str_resnum
                    updated_cluster_residues.append(modified_res)
                else:
                    # If the 4th character is not "A", then append the original residue name
                    # to the list. 
                    updated_cluster_residues.append(res)
            # Update the new list under the original key.
            clust_dict[k] = updated_cluster_residues
        return clust_dict
    
def generate_path_pdb(conformers, source_pdb, prefix= "path"):
    """generate a pdb file from path data
    
    Parameters
    ----------
    conformers : TYPE
        Description
    source_pdb : TYPE
        Description
    prefix : str, optional
        Description
    """
    with open(source_pdb, 'r') as pdb:
        pdb_lines = pdb.readlines()
        with open(prefix + ".pdb", "w") as f:
            for conf in conformers:
                for l in pdb_lines:
                    if l.split()[3] == conf[0:3]:
                        #print l.split()[3], l.split()[4]
                        #print conf[0:3], conf[5:]    
                        if l.split()[4] == conf[5:] or l.split()[4] == conf[5:11] + "000":
                            #print l.split()[3], l.split()[4]
                            #print conf[0:3], conf[5:]
                            f.write(l)