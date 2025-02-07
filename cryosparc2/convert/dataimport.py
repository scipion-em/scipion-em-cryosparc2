import os
import logging
import time

import emtable

from cryosparc2 import RELIONCOLUMNS
from cryosparc2.convert import (convertCs2Star, readSetOfParticles,
                                cryosparcToLocation)
from pwem import ALIGN_PROJ
from pwem.objects import Coordinate, SetOfCoordinates

logger = logging.getLogger(__name__)

class cryoSPARCImport:
    """ Class used to import particles from cryoSPARC projects into Scipion.
    """
    def __init__(self, protocol, csFile):
        self.protocol = protocol
        self._csFile = csFile
        self._createFilenameTemplates()

    def _createFilenameTemplates(self):
        """ Centralize how files are called. """
        myDict = {
            'output': self.protocol._getExtraPath('output.star')
        }
        self.protocol._updateFilenamesDict(myDict)

    def importParticles(self):
        """
        Import particles from a cs 'particles.cs'
        """
        try:
            csPartFile = os.path.abspath(self._csFile)
            self.outputStarFn = self.protocol._getFileName('output')
            argsList = [csPartFile, self.outputStarFn]
            csFileName = os.path.basename(self._csFile)
            csFileDir = os.path.dirname(self._csFile)
            passthroughtFileName = csFileName.split('_')[0] + '_passthrough_' + csFileName.split('_')[1]
            passthroughtFilePath = os.path.join(csFileDir, passthroughtFileName)
            if os.path.exists(passthroughtFilePath):
                argsList.append(passthroughtFilePath)
            convertCs2Star(argsList)
            # Validate the start file generated
            self._validateConvert()

            self.partSet = self.protocol._createSetOfParticles()
            self.partSet.setObjComment('Particles imported from cryosPARC .cs file:\n%s' % self._csFile)

            # Update both samplingRate and acquisition with parameters
            # selected in the protocol form
            self.protocol.setSamplingRate(self.partSet)
            self._pixelSize = self.protocol.samplingRate.get()
            self.partSet.setIsPhaseFlipped(self.protocol.haveDataBeenPhaseFlipped.get())
            self.protocol.fillAcquisition(self.partSet.getAcquisition())

            self._fillDataFromIter(self.partSet)
            self.protocol._defineOutputs(outputParticles=self.partSet)

        except Exception as e:
            raise Exception("The .cs file has not been imported: %s" % e)

    def importCoordinates(self):
        """
        Import coordinates from a cs 'particles.cs'
        """
        micSetPtr = self.protocol.inputMicrographs
        outputCoords = SetOfCoordinates.create(self.protocol.getPath())
        outputCoords.setMicrographs(micSetPtr)

        micList = {os.path.basename(mic.getFileName()): mic.clone() for mic in micSetPtr.get()}

        for fileName, _ in self.protocol.iterFiles():
            if fileName.endswith('.cs'):
                try:
                    csPartFile = os.path.abspath(fileName)
                    outputStarFn = self.protocol._getFileName('output')
                    argsList = [csPartFile, outputStarFn]
                    convertCs2Star(argsList)
                    self._fillSetOfCoordinates(outputCoords, outputStarFn, micList)

                except Exception as e:
                    logger.error("The .cs file has not been imported: %s" % fileName, exc_info=e)

        return outputCoords

    def _fillSetOfCoordinates(self, outputCoords, outputStarFn, micList):

        coord = Coordinate()
        mdFileName = '%s@%s' % ('particles', outputStarFn)
        table = emtable.Table(fileName=outputStarFn)

        for row in table.iterRows(mdFileName):
            coord.setObjId(None)
            micName = os.path.basename(row.get(RELIONCOLUMNS.rlnMicrographName.value))
            splitMicName = micName.split('_')
            if len(splitMicName) > 1:
                micName = '_'.join(splitMicName[1:])
            else:
                micName = splitMicName[-1]
            coord.setMicrograph(micList[micName])
            x = row.get(RELIONCOLUMNS.rlnCoordinateX.value)
            y = row.get(RELIONCOLUMNS.rlnCoordinateY.value)
            dim = micList[micName].getDimensions()
            flipY = dim[1] - y
            coord.setPosition(x, flipY)
            # Add it to the set
            outputCoords.append(coord)

    def _fillDataFromIter(self, imgSet):
        outImgsFn = 'particles@' + self.protocol._getFileName('output')
        readSetOfParticles(outImgsFn, imgSet,
                           postprocessImageRow=self._updateItem,
                           alignType=ALIGN_PROJ,
                           samplingRate=imgSet.getSamplingRate())

    def _updateItem(self, item, row):
        index, file = item.getLocation()
        binaryPath = self.findImagesFrom(self._csFile, file)
        item.setLocation(index, binaryPath)
        item.setSamplingRate(self._pixelSize)

    def _validateConvert(self):
        self._validateMetadata("rlnImageName", warnings=True)

    def validateParticles(self):
        errors = []
        return errors

    def _validateMetadata(self, label, warnings=True):
        # read data_particles table
        table = emtable.Table(fileName=self.outputStarFn, tableName='particles')
        row = table[0]

        if row is None:
            raise Exception("Cannot import from empty metadata: %s"
                            % self.outputStarFn)

        if not row.get(label, False):
            raise Exception("Label *%s* is missing in metadata: %s"
                            % (label, self.outputStarFn))

        index, fn = cryosparcToLocation(row.get(label))
        if fn.startswith("/"):
            raise Exception("ERROR: %s cannot be an absolute path: %s\n"
                            "Please create a symlink to an abs path instead."
                            % (label, fn))

        self._imgPath = self.findImagesFrom(self._csFile, fn)

        if warnings and self._imgPath is None:
            self.protocol.warning("WARNING: Binary data was not found from metadata: %s"
                                  % self.outputStarFn)
        return row

    def findImagesFrom(self, referenceFile, searchFile):
        separador = os.path.sep
        dirParticlesPath = os.path.dirname(os.path.abspath(referenceFile))
        dirProject = separador.join(dirParticlesPath.split(separador)[:-1])
        imageName = os.path.basename(searchFile)
        imageFolder = os.path.dirname(searchFile)
        realdataPath = dirParticlesPath
        folderContent = []

        # Case in which the binaries associated to the .cs file are in the same
        # folder of the .cs file
        filesPath = dirParticlesPath
        if os.path.exists(filesPath):
            folderContent += os.listdir(filesPath)

        # Case in which the binaries associated to the .cs file are located in
        # the cryoSPARC folders format
        filesPath2 = os.path.join(dirProject, imageFolder)
        if os.path.exists(filesPath2):
            folderContent += os.listdir(filesPath2)
            realdataPath = filesPath2

        matchFile = [f for f in folderContent if f.endswith(imageName)]
        if matchFile:
            file = os.path.join(realdataPath, matchFile[0])
            return file
        else:
            raise Exception("The images were expected to be found in one of "
                            "these locations: (%s) or (%s)" % (filesPath,
                                                               filesPath2))






