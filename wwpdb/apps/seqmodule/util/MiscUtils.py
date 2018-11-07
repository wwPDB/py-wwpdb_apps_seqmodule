##
# File:    MiscUtils.py
# Date:    24-Feb-2013
#
# Updates:
##
"""
Miscellaneous utilites.

"""
__docformat__ = "restructuredtext en"
__author__    = "John Westbrook"
__email__     = "jwest@rcsb.rutgers.edu"
__license__   = "Creative Commons Attribution 3.0 Unported"
__version__   = "V0.07"

import sys

def editDistance(s1, s2):
    """
    Compute the Damerau-Levenshtein distance between two given
    strings (s1 and s2)
    """    
    d = {}
    lenstr1 = len(s1)
    lenstr2 = len(s2)
    for i in xrange(-1,lenstr1+1):
        d[(i,-1)] = i+1
    for j in xrange(-1,lenstr2+1):
        d[(-1,j)] = j+1
 
    for i in xrange(0,lenstr1):
        for j in xrange(0,lenstr2):
            if s1[i] == s2[j]:
                cost = 0
            else:
                cost = 1
            d[(i,j)] = min(
                           d[(i-1,j)] + 1, # deletion.
                           d[(i,j-1)] + 1, # insertion
                           d[(i-1,j-1)] + cost, # substitution
                          )
            if i>1 and j>1 and s1[i]==s2[j-1] and s1[i-1] == s2[j]:
                d[(i,j)] = min (d[(i,j)], d[i-2,j-2] + cost) # transposition
 
    return d[lenstr1-1,lenstr2-1]

def multikeysort(items, columns):
    """
    Sort list of dictionaries on multiple keys -

    Example 
    uu = [{"a": 1, "b": 2, "c": "tiger" },
    {"a": 2, "b": 1, "c": "tiger" },
    {"a": 3, "b": 5, "c": "bear" },
    {"a": 4, "b": 4, "c": "tiger" },
    {"a": 5, "b": 1, "c": "bear" }
    ]
    result = multikeysort(uu, ['c', 'b', 'a'])
    """    
    from operator import itemgetter
    comparers = [ ((itemgetter(col[1:].strip()), -1) if col.startswith('-') else (itemgetter(col.strip()), 1)) for col in columns]
    def sign(a, b):
        if   a < b:  return -1
        elif a > b:  return 1
        else:        return 0
    def comparer(left, right):
        for fn, mult in comparers:
            result = sign(fn(left), fn(right))
            if result:
                return mult * result
        else:
            return 0
    return sorted(items, cmp=comparer)

