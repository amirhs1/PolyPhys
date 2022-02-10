"""
axuiliary 
"""

import re
import math

def sort_by_int(alphanumeric):
    """
    sort_by_int split an alphanumeric into words and integers.

    Reference:
    http://code.activestate.com/recipes/135435-sort-a-string-using-numeric-order/

    Parameters:
    alphanumeric (char): an alphanumeric string.

    Return:
    a mixed list of words and integers.

    Requirement:
    re
    """

    pattern = re.compile('(\d+)')
    split = pattern.split(alphanumeric)
    mixed = [int(word) if word.isdigit() else word for word in split]
    return mixed


def file_reader(f_names,extensions = ["data",("lammpstrj","trj")], warning=False):
    """
    file_reader returns a sorted list of sorted file names based on the agg function. A trajectroy (trj and lammpstrj extention) constists of N snapshots and a topology (data extention) contains information about simulation box, particle and bond types.

    Parameters:
    f_names (A list of strings): A list of pathes to Lammps trajectroies (trj and lammpstrj extention) and topology (data extention) files.
    
    Return:
    a sorted list of tuples where each tuple has len(extensions) members.

    Requirements:
    re
    """
    f_names_sorted = [None] * len(extensions) # an empty list of len(extensions)
    # a nested list where each sublist has all the files with the same extension:
    for idx, extension in enumerate(extensions): 
        f_names_sorted[idx] = [f_name for f_name in f_names if f_name.endswith(extension)]
    
    for idx, names_same_ext in enumerate(f_names_sorted):
        f_names_sorted[idx] = sorted(names_same_ext, key = sort_by_int) # Sorting pathnames alphanumerically
    
    f_names = list(zip(*f_names_sorted))
    if warning:
        print("Total number of files is ", len(f_names))
        print("Path to the first tuple of the sorted file: ", f_names[0])
    return f_names


def round_down_first_non_zero(x):
    """
    rounds down a number to its first non-zero digit.
    
    Parameters:
    x (float or int): number which is rounded.
    
    Returns:
    a float or integer number to which to which x is rounded down to its first non-zero digit.
    """
    if x == 0:
        return x
    else:
        exponent = math.floor(math.log10(abs(x)))
        non_zero = 10 ** exponent
        return round(math.floor(x/non_zero)*non_zero, abs(exponent))
    
def round_up_nearest(dividend, diviser):
    """
    rounds up the floating-point number dividend by the diviser.
    
    Parameters:
    dividend (float): The number should be rounded up.
    diviser (float): The number used as the diviser.
    
    Return
    A floating-point number that is a ceil for dividend and is divisible by diviser. 
    """
    return math.ceil(dividend / diviser) * diviser