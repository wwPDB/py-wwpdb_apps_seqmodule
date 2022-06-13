##
# File:  MultiProcLimit.py
# Date:  02-Dev-2018
# Updates:
##
"""
Providing a mechanism to ensure that MultiProcUtil will wait for limited time before proceeding

This software was developed as part of the World Wide Protein Data Bank
Common Deposition and Annotation System Project

Copyright (c) 2019 wwPDB

This software is provided under a Creative Commons Attribution 3.0 Unported
License described at http://creativecommons.org/licenses/by/3.0/.

"""
__docformat__ = "restructuredtext en"
__author__ = "Ezra Peisach"
__email__ = "peisach@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

from multiprocessing import Value
import time
import logging

logger = logging.getLogger()


class MultiProcLimit(object):
    """A class that will ensure that we limit requests made to service are less than rate/second.
    This is done by ensuring that each request is > 1/rate seconds from previous."""

    def __init__(self, rate):
        self.__rate = rate
        self.__delay = (1.0 / rate) + 0.0001
        self.__last = Value("d", 0.0)
        logger.debug("Instantiated with rate %r delay %r ", self.__rate, self.__delay)

    def waitnext(self):
        """Waits for the next time to fire"""

        # Acquire lock on value - will block other processes
        with self.__last.get_lock():
            current = time.time()
            wait = self.__last.value + self.__delay - current
            if wait > 0:
                # logger.debug("About to sleep %s" % wait)
                time.sleep(wait)
            self.__last.value = time.time()
