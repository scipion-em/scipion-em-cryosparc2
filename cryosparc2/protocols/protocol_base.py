# **************************************************************************
# *
# * Authors: Yunior C. Fonseca Reyna    (cfonseca@cnb.csic.es)
# *
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
import ast
import requests
from pkg_resources import parse_version

import pwem.protocols as pw
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pwem.objects import FSC

from ..constants import V3_0_0, V3_3_1
from ..convert import convertBinaryVol, writeSetOfParticles, ImageHandler
from ..utils import (getProjectPath, createEmptyProject,
                     createEmptyWorkSpace, getProjectName,
                     getCryosparcProjectsDir, createProjectDir,
                     doImportParticlesStar, doImportVolumes, killJob, clearJob,
                     get_job_streamlog, getSystemInfo, getJobStatus,
                     STOP_STATUSES, getCryosparcVersion)


class ProtCryosparcBase(pw.EMProtocol):
    """
    This class contains the common functions for all Cryosparc protocols.
    """
    _protCompatibility = []
    _className = ""
    _fscColumns = 6

    def _initializeCryosparcProject(self):
        """
        Initialize the cryoSPARC project and workspace
        """
        self._initializeUtilsVariables()
        # create empty project or load an exists one
        folderPaths = getProjectPath(self.projectPath)
        if not folderPaths:
            self.emptyProject = createEmptyProject(self.projectPath, self.projectDirName)
            self.projectName = self.emptyProject[-1].split()[-1]
        else:
            self.projectName = str(folderPaths[0])

        self.projectName = pwobj.String(self.projectName)
        self._store(self)

        # create empty workspace
        self.emptyWorkSpace = createEmptyWorkSpace(self.projectName, self.getRunName(),
                                      self.getObjComment())
        self.workSpaceName = pwobj.String(self.emptyWorkSpace[-1].split()[-1])
        self._store(self)

        self.currenJob = pwobj.String()
        self._store(self)

    def _initializeUtilsVariables(self):
        """
        Initialize all utils cryoSPARC variables
        """
        # Create a cryoSPARC project dir
        self.projectDirName = getProjectName(self.getProject().getShortName())
        self.projectPath = pw.pwutils.join(getCryosparcProjectsDir(),
                                        self.projectDirName)
        self.projectDir = createProjectDir(self.projectPath)

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self._getInputParticles()
        if imgSet is not None:
            # Create links to binary files and write the relion .star file
            writeSetOfParticles(imgSet, self._getFileName('input_particles'),
                                self._getTmpPath())
            self._importParticles()

        volume = self._getInputVolume()
        if volume is not None:
            self._importVolume()

        mask = self._getInputMask()
        if mask is not None:
            self._importMask()
        else:
            self.mask = pwobj.String()
        self._store(self)

    def _getScaledAveragesFile(self, csAveragesFile, force=False):

        # For the moment this is the best possible result, scaling from 128 to
        # 300 does not render nice results apart that the factor turns to
        # 299x299. But without this the representative subset is wrong.
        # return csAveragesFile

        scaledFile = self._getScaledAveragesFileName(csAveragesFile, force)

        if not os.path.exists(scaledFile):

            inputSize = self._getInputParticles().getDim()[0]
            csSize = ImageHandler().getDimensions(csAveragesFile)[0]

            if csSize == inputSize:
                print("No binning detected: linking averages cs file.",
                      flush=True)
                pwutils.createLink(csAveragesFile, scaledFile)
            else:
                print("Scaling CS averages file to match particle "
                      "size (%s -> %s)." % (csSize, inputSize), flush=True)
                try:
                    if force:
                        scaleFactor = inputSize/csSize
                        ImageHandler.scaleSplines(csAveragesFile, scaledFile,
                                                  scaleFactor,
                                                  finalDimension=inputSize,
                                                  forceVolume=force)
                    else:

                        ImageHandler.scale2DStack(csAveragesFile, scaledFile,
                                                  finalDimension=inputSize)
                except Exception as ex:
                    print("The CS averages could not be scaled. %s " % ex)
                    return csAveragesFile

        return scaledFile

    def _getScaledAveragesFileName(self, csAveragesFile, isVolume=False):

        extension = ".mrc" if isVolume else ".mrcs"
        return pwutils.removeExt(csAveragesFile) + "_scaled" + extension

    def setFilePattern(self, path):
        baseName = os.path.basename(path).split('.')[0]
        self.inputFileNamePattern = path.replace(baseName, '%s')

    def updateParticlePath(self, part, row):
        fn = part.getFileName()
        baseName = os.path.basename(fn).split('.')[0]
        newFileName = self.inputFileNamePattern % baseName
        part.setFileName(newFileName)

    def _getInputParticles(self):
        if self.hasAttribute('inputParticles'):
            return self.inputParticles.get()
        return None

    def _getInputVolume(self):
        if self.hasAttribute('refVolume'):
            return self.refVolume.get()
        return None

    def _getInputMask(self):
        if self.hasAttribute('refMask'):
            return self.refMask.get()
        return None

    def _initializeVolumeSuffix(self):
        """
        Create a output volume suffix depend of the CS version
        """
        cryosparcVersion = parse_version(getCryosparcVersion())
        self.outputVolumeSuffix = '.imported_volume.map'
        self.outputMaskSuffix = '.imported_mask.map'
        self.outputVolumeHalf_A = '.imported_volume.map_half_A'
        self.outputVolumeHalf_B = '.imported_volume.map_half_B'
        if cryosparcVersion >= parse_version(V3_3_1):
            self.outputVolumeSuffix = '.imported_volume_1.map'
            self.outputMaskSuffix = '.imported_mask_1.map'
            self.outputVolumeHalf_A = '.imported_volume_1.map_half_A'
            self.outputVolumeHalf_B = '.imported_volume_1.map_half_B'

    def _initializeMaskSuffix(self):
        """
        Create a output mask suffix depend of the CS version
        """
        cryosparcVersion = parse_version(getCryosparcVersion())
        self.outputMaskSuffix = '.imported_mask.map'
        if cryosparcVersion >= parse_version(V3_3_1):
            self.outputMaskSuffix = '.imported_mask_1.map'

    def _importVolume(self):
        vol = self._getInputVolume()
        self._initializeVolumeSuffix()
        vol_fn = os.path.join(os.getcwd(), convertBinaryVol(vol, self._getTmpPath()))
        importVolumeJob = doImportVolumes(self, vol_fn, vol, 'map', 'Importing volume...')
        self.volume = pwobj.String(str(importVolumeJob.get()) + self.outputVolumeSuffix)

        if vol.hasHalfMaps():
            halfMaps = vol.getHalfMaps().split(",")
            map_half_A_fn = os.path.abspath(halfMaps[0])
            importVolumeHalfAJob = doImportVolumes(self, map_half_A_fn, vol,
                                                   'map_half_A', 'Importing half volume A...')
            self.importVolumeHalfA = pwobj.String(str(importVolumeHalfAJob.get()) + self.outputVolumeHalf_A)

            map_half_B_fn = os.path.abspath(halfMaps[1])
            importVolumeHalfBJob = doImportVolumes(self, map_half_B_fn, vol,
                                                   'map_half_B', 'Importing half volume B...')
            self.importVolumeHalfB = pwobj.String(str(importVolumeHalfBJob.get()) + self.outputVolumeHalf_B)

        self.currenJob.set(importVolumeJob.get())

    def _importMask(self):
        self._initializeMaskSuffix()
        maskFn = os.path.join(os.getcwd(), convertBinaryVol(self._getInputMask(),
                                                            self._getTmpPath()))

        importMaskJob = doImportVolumes(self, maskFn, self._getInputMask(),
                                        'mask', 'Importing mask... ')
        self.currenJob.set(importMaskJob.get())
        self.mask = pwobj.String(str(importMaskJob.get()) + self.outputMaskSuffix)

    def _importParticles(self):

        # import_particles_star
        importedParticlesJob = doImportParticlesStar(self)
        self.currenJob = pwobj.String(str(importedParticlesJob.get()))
        self.particles = pwobj.String(str(importedParticlesJob.get()) +
                                      '.imported_particles')

    def setAborted(self):
        """ Set the status to aborted and updated the endTime. """
        pw.EMProtocol.setAborted(self)
        if hasattr(self, 'projectName') and hasattr(self, 'currenJob') and self.currenJob.get() is not None:
            status = getJobStatus(self.projectName.get(), self.currenJob.get())
            if status not in STOP_STATUSES:
                killJob(str(self.projectName.get()), str(self.currenJob.get()))
                clearJob(str(self.projectName.get()), str(self.currenJob.get()))

    def createFSC(self, idd, imgSet, vol):
        # Need to get the cryosparc master address
        system_info = getSystemInfo()
        status_errors = system_info[0]
        if not status_errors:
            system_info = eval(system_info[1])
            master_hostname = system_info.get('master_hostname')
            port_webapp = system_info.get('port_webapp')

            url = "http://%s:%s/file/%s" % (master_hostname, port_webapp, idd)
            fscRequest = requests.get(url, allow_redirects=True)
            fscFile = "fsc.txt"
            fscFilePath = os.path.join(self._getExtraPath(), fscFile)

            # Convert into scipion fsc format
            open(fscFilePath, 'wb').write(fscRequest.content)
            f = open(fscFilePath, 'r')
            lines = f.readlines()
            wv = []
            corr = []
            factor = self._getInputParticles().getDim()[0] * imgSet.getSamplingRate()
            for x in lines[1:-1]:
                wv.append(str(float(x.split('\t')[0])/factor))
                corr.append(x.split('\t')[self._fscColumns])
            f.close()
            fsc = FSC(objLabel=self.getRunName())
            fsc.setData(wv, corr)
            self._defineOutputs(outputFSC=fsc)
            self._defineSourceRelation(vol, fsc)

    def findLastIteration(self, jobName):
        get_job_streamlog(self.projectName.get(),
                          jobName,
                          self._getFileName('stream_log'))

        # Get the metadata information from stream.log
        with open(self._getFileName('stream_log')) as f:
            data = f.readlines()

        x = ast.literal_eval(data[0])

        # Find the ID of last iteration and the map resolution
        for y in x:
            if 'text' in y:
                z = str(y['text'])

                if z.startswith('FSC Iteration') or z.startswith('FSC iIteration'):
                    idd = y['imgfiles'][2]['fileid']
                    itera = z.split(',')[0][-3:]
                    self._store(self)
                elif 'Estimated Bfactor' in z:
                    estBFactor = str(y['text']).split(':')[1].replace('\n',
                                                                      '')
                    self.estBFactor = pwobj.String(estBFactor)
                    self._store(self)
        return idd, itera

    def _createModelFile(self):
        pass
