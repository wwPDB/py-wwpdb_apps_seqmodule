##
# File:    Autodict.py
# Date:    27-Feb-2013
#
# Updates:
##
"""
Derived dictionary classes supporting automatic initialization and 'dot' access.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"


class Autodict(dict):
    """Derived dictionary class supporting automatic initialization.

    This will support pickle serialization/deserialization.
    """

    def __getitem__(self, name):
        # if not name in self:
        #    dict.__setitem__(self, name, Autodict())
        # return dict.__getitem__(self, name)
        try:
            return dict.__getitem__(self, name)
        except KeyError:
            value = self[name] = type(self)()
            return value


class AutodictDot(dict):
    """Derived dictionary class supporting automatic initialization and 'dot' access.

    This will not support pickle'ng  jdw
    """

    def __getitem__(self, name):
        # if not name in self:
        #    dict.__setitem__(self, name, Autodict())
        # return dict.__getitem__(self, name)
        try:
            return dict.__getitem__(self, name)
        except KeyError:
            value = self[name] = type(self)()
            return value

    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            return self[name]

    def __setattr__(self, name, val):
        if name in self.__dict__:
            self.__dict__[name] = val
        else:
            self[name] = val
