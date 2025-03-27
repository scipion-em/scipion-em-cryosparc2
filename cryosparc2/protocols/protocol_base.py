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
import itertools
import os
import ast
import time

import requests
import logging
logger = logging.getLogger(__name__)

from pkg_resources import parse_version

import pwem.protocols as pw
import pyworkflow.object as pwobj
import pyworkflow.utils as pwutils
from pwem.objects import FSC

from ..constants import V3_3_1, excludedFSCValues, fscValues, V4_0_0, V4_1_0, RELIONCOLUMNS
from ..convert import convertBinaryVol, writeSetOfParticles, ImageHandler
from ..utils import (getProjectPath, createEmptyProject,
                     createEmptyWorkSpace, getProjectName,
                     getCryosparcProjectsDir, createProjectContainerDir,
                     doImportParticlesStar, doImportVolumes, killJob, clearJob,
                     get_job_streamlog, getSystemInfo, getJobStatus,
                     STOP_STATUSES, getCryosparcVersion, getProjectInformation,
                     getCryosparcProjectId, _getLicenceFromFile, doImportMicrographs, getCryosparcProjectsList,
                     getCryosparcWorkSpaces)


class ProtCryosparcBase(pw.EMProtocol):
    """
    This class contains the common functions for all Cryosparc protocols.
    """
    _protCompatibility = []
    _className = ""
    _fscColumns = 6
    _logLastLine = 0

    def _initializeCryosparcProject(self):
        """
        Initialize the cryoSPARC project and workspace
        """
        self._initializeUtilsVariables()
        projectsList = getCryosparcProjectsList()
        matchProjects = [project for project in projectsList if project.get('title') == self.projectDirName]
        folderPaths = getProjectPath(self.projectContainerDir)
        # create an empty project or load an exists one
        if not matchProjects or not folderPaths:
            # create an empty project
            self.emptyProject = createEmptyProject(self.projectPath, self.projectDirName)
            self.projectName = pwobj.String(self.emptyProject[-1].split()[-1])
            self.projectDir = pwobj.String(getProjectInformation(self.projectName,
                                           info='project_dir'))
            # create an empty workspace
            self.emptyWorkSpace = createEmptyWorkSpace(self.projectName, self.getRunName(),
                                                       self.getObjComment())
            self.workSpaceName = pwobj.String(self.emptyWorkSpace[-1].split()[-1])
            self._store(self)
        else:
            self.projectDir = pwobj.String(matchProjects[-1]['project_dir'])
            cryosparcVersion = getCryosparcVersion()
            if parse_version(cryosparcVersion) < parse_version(V4_0_0):
                self.projectName = pwobj.String(matchProjects[-1]['title'])
            else:
                self.projectName = pwobj.String(matchProjects[-1]['uid'])

            workspacesList = getCryosparcWorkSpaces(str(self.projectName))
            self.workSpaceName = pwobj.String(workspacesList[-1]['uid'])

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
        self.projectContainerDir = createProjectContainerDir(self.projectPath)[1]

    def convertInputStep(self):
        """ Create the input file in STAR format as expected by Relion.
        If the input particles comes from Relion, just link the file.
        """
        imgSet = self._getInputParticles()
        if imgSet is not None:
            # Create links to binary files and write the relion .star file
            writeSetOfParticles(imgSet, self._getFileName('input_particles'),
                                self._getPath())
            self._importParticles()

        volume = self._getInputVolume()
        if volume is not None:
            self._importVolume()

        mask = self._getInputMask()
        if mask is not None:
            self._importMask()
        else:
            self.mask = pwobj.String()

        focusMask = self._getInputFocusMask()
        if focusMask is not None:
            self._importFocusMask()
        else:
            self.focusMask = pwobj.String()

        micrographs = self._getInputMicrographs()
        if micrographs is not None:
            self._importMicrographs()

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
                self.info("No binning detected: linking averages cs file.")
                pwutils.createLink(csAveragesFile, scaledFile)
            else:
                self.info("Scaling CS averages file to match particle "
                      "size (%s -> %s)." % (csSize, inputSize))
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
                    self._log.error("The CS averages could not be scaled. %s ", exc_info=ex)
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

    def _getInputParticlesPointer(self):
        if self.hasAttribute('inputParticles'):
            return self.inputParticles
        return None

    def _getInputVolume(self):
        if self.hasAttribute('refVolume'):
            return self.refVolume.get()
        return None

    def _getInputMask(self):
        if self.hasAttribute('refMask'):
            return self.refMask.get()
        return None

    def _getInputFocusMask(self):
        if self.hasAttribute('refFocusMask'):
            return self.refFocusMask.get()
        return None

    def _getInputMicrographs(self):
        if self.hasAttribute('inputMicrographs'):
            return self.inputMicrographs.get()
        return None

    def _initializeVolumeSuffix(self):
        """
        Create an output volume suffix depend on the CS version
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

    def _initializeMaskSuffix(self, sufix='.imported_mask_1.map'):
        """
        Create a output mask suffix depend of the CS version
        """
        cryosparcVersion = parse_version(getCryosparcVersion())
        self.outputMaskSuffix = '.imported_mask.map'
        if cryosparcVersion >= parse_version(V3_3_1):
            self.outputMaskSuffix = sufix

    def _importVolume(self):
        vol = self._getInputVolume()
        self._initializeVolumeSuffix()
        vol_fn = os.path.join(os.getcwd(), convertBinaryVol(vol, self._getTmpPath()))
        importVolumeJob = doImportVolumes(self, vol_fn, vol, 'map', 'Importing volume...')
        self.volume = pwobj.String(str(importVolumeJob.get()) + self.outputVolumeSuffix)

        if vol.hasHalfMaps():
            halfMaps = vol.getHalfMaps().split(",")
            map_half_A_fn = os.path.abspath(halfMaps[0].split(':mrc')[0])
            importVolumeHalfAJob = doImportVolumes(self, map_half_A_fn, vol,
                                                   'map_half_A', 'Importing half volume A...')
            self.importVolumeHalfA = pwobj.String(str(importVolumeHalfAJob.get()) + self.outputVolumeHalf_A)

            map_half_B_fn = os.path.abspath(halfMaps[1].split(':mrc')[0])
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

    def _importFocusMask(self):
        self._initializeMaskSuffix()
        maskFn = os.path.join(os.getcwd(), convertBinaryVol(self._getInputFocusMask(),
                                                            self._getTmpPath()))

        importFocusMaskJob = doImportVolumes(self, maskFn, self._getInputFocusMask(),
                                             'mask', 'Importing focus mask... ')
        self.currenJob.set(importFocusMaskJob.get())
        self.focusMask = pwobj.String(str(importFocusMaskJob.get()) + self.outputMaskSuffix)

    def _importParticles(self):
        # import_particles_star
        importedParticlesJob = doImportParticlesStar(self)
        self.currenJob = pwobj.String(str(importedParticlesJob.get()))
        self.particles = pwobj.String(str(importedParticlesJob.get()) +
                                      '.imported_particles')

    def _importMicrographs(self):
        importedMicrographsJob = doImportMicrographs(self)
        self.currenJob = pwobj.String(str(importedMicrographsJob.get()))
        self.micrographs = pwobj.String(str(importedMicrographsJob.get()) +
                                      '.imported_micrographs')

    def setAborted(self):
        """ Set the status to aborted and updated the endTime. """
        pw.EMProtocol.setAborted(self)
        if hasattr(self, 'projectName') and hasattr(self, 'currenJob') and self.currenJob.get() is not None:
            job = str(self.currenJob.get())
            project = str(self.projectName.get())
            status = getJobStatus(project, job)
            if status not in STOP_STATUSES:
                try:
                    killJob(project, job)
                    clearJob(project, job)
                except Exception as e:
                    logger.error("Can't kill job %s from project %s" % (job, project), exc_info=e)

    def createFSC(self, idd, imgSet, vol):
        # Need to get the cryosparc master address
        system_info = getSystemInfo()
        status_errors = system_info[0]

        if not status_errors:
            cryosparcVersion = getCryosparcVersion()
            system_info = eval(system_info[1])
            master_hostname = system_info.get('master_hostname')
            if parse_version(cryosparcVersion) < parse_version(V4_1_0):
                port_webapp = system_info.get('port_webapp')
                url = "http://%s:%s/file/%s" % (master_hostname, port_webapp, idd)
                fscRequest = requests.get(url, allow_redirects=True)
            else:
                port_webapp = system_info.get('port_command_vis')
                url = "http://%s:%s/get_job_file" % (master_hostname, port_webapp)
                jsonParam = {'fileid': idd}
                licence_id = _getLicenceFromFile()
                headers = {'License-ID': licence_id}
                fscRequest = requests.post(url, json=jsonParam, headers=headers,
                                          allow_redirects=True)
            fscFile = "fsc.txt"
            fscFilePath = os.path.join(self._getExtraPath(), fscFile)
            factor = self._getInputParticles().getDim()[0] * imgSet.getSamplingRate()

            # Convert into scipion fsc format
            open(fscFilePath, 'wb').write(fscRequest.content)
            fscSet = self.getSetOfFCSsFromFile(fscFilePath, factor)
            self._defineOutputs(outputFSC=fscSet)
            self._defineSourceRelation(vol, fscSet)

    def getSetOfFCSsFromFile(self, file, factor):
        f = open(file, 'r')
        lines = f.readlines()
        fscSet = self._createSetOfFSCs()
        columns = lines[0].strip().split('\t')[1:]
        col = 1
        fsc_t = None
        fsc_nt = None

        for column in columns:
            if column == 'fsc_tightmask':
                fsc_t = \
                self.getFSCFromRawData(lines, column, col, factor).getData()[1]
            if column == 'fsc_noisesub_true':
                fsc_nt = \
                self.getFSCFromRawData(lines, column, col, factor).getData()[1]
            if column not in excludedFSCValues:
                fsc = self.getFSCFromRawData(lines, column, col, factor)
                fscSet.append(fsc)
            col += 1
        f.close()
        corr = []

        if fsc_t is not None and fsc_nt is not None:  # Phase Randomized Masket Map can be calculated
            for i in range(len(fsc_t)):
                fsc_nt_value = 0.99
                if fsc_nt[i] != 1:
                     fsc_nt_value = fsc_nt[i]
                corr.append((fsc_t[i] - fsc_nt[i]) / (1.0 - fsc_nt_value))
            fsc = FSC(objLabel=fscValues['fsc_prmm'])
            fsc_wv = fscSet.getFirstItem().getData()[0]
            fsc.setData(fsc_wv, corr)
            fscSet.append(fsc)
        fscSet.write()
        return fscSet

    def getFSCFromRawData(self, lines, label, col, factor):
        wv = []
        corr = []
        for x in lines[1:]:
            wv_value = float(x.strip().split('\t')[0])
            coor_value = x.strip().split('\t')[col]
            wv.append(str(wv_value / factor))
            corr.append(coor_value)
        fsc = FSC(objLabel=fscValues[label])
        fsc.setData(wv, corr)
        return fsc

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
                    for imgfile in y['imgfiles']:
                        if imgfile['filetype'] == 'txt':
                            idd = imgfile['fileid']
                            break
                    itera = z.split(',')[0][-3:]
                    self._store(self)
                elif 'Using Filter Radius' in z:
                    nomRes = str(y['text']).split('(')[1].split(')')[
                        0].replace(
                        'A', 'Å')
                    self.mapResolution = pwobj.String(nomRes)
                    self._store(self)
                elif 'Estimated Bfactor' in z:
                    estBFactor = str(y['text']).split(':')[1].replace('\n',
                                                                      '')
                    self.estBFactor = pwobj.String(estBFactor)
                    self._store(self)
        return idd, itera

    def _createModelFile(self):
        pass

    def getLogLine(self):
        return self._logLastLine

    def setLogLine(self, lastLine: int):
        self._logLastLine = lastLine
