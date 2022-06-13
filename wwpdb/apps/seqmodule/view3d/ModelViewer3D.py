##
# File:  ModelViewer3D.py
# Date:  8-Feb-2010
#
# Update:
# 20-Apr-2010 jdw Port to module seqmodule.
#  2-Dec-2013 jdw Only use the cif file for Jmol here --   Astex is now obsoleted.
#  2-Dec-2013 jdw tested with the new signed applet --
##
"""
Utility methods model and map 3d visualization.

"""
__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Creative Commons Attribution 3.0 Unported"
__version__ = "V0.07"

import sys
import os.path
import os
import shutil
import traceback
from wwpdb.io.locator.PathInfo import PathInfo

from wwpdb.apps.seqmodule.io.SequenceDataStore import SequenceDataStore


class ModelViewer3D(object):
    """
    ModelViewer3D class encapsulates depiction of coordinate models and maps.

    """

    def __init__(self, reqObj=None, verbose=False, log=sys.stderr):
        self.__verbose = verbose
        self.__reqObj = reqObj
        self.__lfh = log
        #
        self.__sessionObj = None
        self.__xyzPathRel = ""
        self.__setup()

    def __setup(self):
        try:
            self.__siteId = self.__reqObj.getValue("WWPDB_SITE_ID")
            self.__sessionObj = self.__reqObj.getSessionObj()
            self.__sessionId = self.__sessionObj.getId()
            self.__sessionPath = self.__sessionObj.getPath()
            self.__identifier = str(self.__reqObj.getValue("identifier")).upper()
            self.__pI = PathInfo(siteId=self.__siteId, sessionPath=self.__sessionPath, verbose=self.__verbose, log=self.__lfh)
            pdbxFilePath = self.__pI.getModelPdbxFilePath(self.__identifier, fileSource="session")
            pdbFilePath = self.__pI.getModelPdbFilePath(self.__identifier, fileSource="session")
            simpleFileName = self.__identifier + ".cif"
            self.__xyzPathRel = os.path.join("/sessions", self.__sessionId, simpleFileName)

            if False:  # pylint: disable=using-constant-test

                if not os.access(self.__xyzPathRel, os.R_OK):
                    shutil.copyfile(pdbxFilePath, os.path.join(self.__sessionPath, simpleFileName))

                if os.access(pdbFilePath, os.R_OK):
                    simpleFileName = self.__identifier + ".pdb"
                    shutil.copyfile(pdbFilePath, os.path.join(self.__sessionPath, simpleFileName))
                    self.__xyzPathRel = os.path.join("/sessions", self.__sessionId, simpleFileName)

                sds = SequenceDataStore(reqObj=self.__reqObj, verbose=self.__verbose, log=self.__lfh)
                self.__rcsbId = str(sds.getEntryDetail("RCSB_ID")).lower()
                self.__pdbId = str(sds.getEntryDetail("PDB_ID")).lower()

                self.__xyzPathRel = ""
                if len(self.__rcsbId):
                    tPath = os.path.join(self.__sessionPath, self.__rcsbId + ".pdb")
                    self.__lfh.write("+ModelViewer3D.__setup() -- value of rcsb tPath is %s\n" % (tPath))
                    if os.access(tPath, os.R_OK):
                        self.__xyzPathRel = os.path.join("/sessions", self.__sessionId, self.__rcsbId + ".pdb")

                if len(self.__pdbId):
                    tPath = os.path.join(self.__sessionPath, self.__pdbId + ".pdb")
                    self.__lfh.write("+ModelViewer3D.__setup() -- value of pdb tPath is %s\n" % (tPath))
                    if os.access(tPath, os.R_OK):
                        self.__xyzPathRel = os.path.join("/sessions", self.__sessionId, self.__pdbId + ".pdb")
            if self.__verbose:
                self.__lfh.write("+ModelViewer3D.__setup()  sessionId %s xyz path %s\n" % (self.__sessionId, self.__xyzPathRel))
        except:  # noqa: E722 pylint: disable=bare-except
            if self.__verbose:
                self.__lfh.write("+ModelViewer3D.__setup()  sessionId %s failed\n" % (self.__sessionObj.getId()))
                traceback.print_exc(file=self.__lfh)

    def launch2(self):
        htmlL = []
        htmlL.append('<script type="text/javascript">')
        htmlL.append('jmolApplet(300, "load %s");' % self.__xyzPathRel)
        htmlL.append("</script>")
        return htmlL

    def launchJmol(self):
        setupCmds = "background black; wireframe only; wireframe 0.05; labels off; slab 100; depth 40; slab on;"

        htmlL = []
        # htmlL.append('<applet name="jmolApplet0" id="jmolApplet0" code="JmolApplet" archive="JmolApplet0.jar" codebase="/applets/jmol" mayscript="true" height="100%" width="100%">')
        htmlL.append(
            '<applet name="jmolApplet0" id="jmolApplet0" code="JmolApplet" archive="JmolAppletSigned0.jar" codebase="/applets/jmol-dev/jsmol/java" mayscript="true" height="100%" width="100%">'  # noqa: E501
        )
        htmlL.append('<param name="progressbar" value="true">')
        htmlL.append('<param name="progresscolor" value="blue">')
        htmlL.append('<param name="boxbgcolor" value="white">')
        htmlL.append('<param name="boxfgcolor" value="black">')
        htmlL.append('<param name="boxmessage" value="Downloading JmolApplet ...">')

        htmlL.append('<param name="script" value="load %s; %s">' % (self.__xyzPathRel, setupCmds))
        htmlL.append("</applet>")
        myD = {}
        myD["htmlcontent"] = str("".join(htmlL))
        return myD

    def launchAstex(self):
        astexScript = r"""
        background '0x000000';
        set symmetry off;
        molecule load mol '%s';
        display spheres off all;         display cylinders off all;         display sticks off all;
        display lines on all;
        color_by_atom;
        object display 'mol_schematic' on;
        """

        setupCmds = astexScript % self.__xyzPathRel
        htmlL = []
        htmlL.append('<applet   id="astexviewer"    name="astexviewer" width="90%" height="90%" code="MoleculeViewerApplet.class" archive="/applets/astex/OpenAstexViewer.jar">')

        htmlL.append('<param name="script" value="%s">' % setupCmds)
        htmlL.append("</applet>")
        # htmlL.extend(self.astexControls())
        myD = {}
        myD["htmlcontent"] = str("".join(htmlL))
        return myD

    def astexControls(self):
        oL = []
        oL.append('<form name="appletcontrols" onsubmit="return false">')
        oL.append('<input  name="avdisplaystyle" checked="checked"  onClick="av_wireframe(this)"  type="radio">Wireframe')
        oL.append('<input  name="avdisplaystyle"                    onClick="av_sticks(this)"     type="radio">Sticks')
        oL.append("</form>")
        return oL
