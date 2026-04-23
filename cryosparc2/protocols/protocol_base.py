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
import logging
import re
from pathlib import Path
import tempfile
logger = logging.getLogger(__name__)

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
                     getCryosparcWorkSpaces, parse_version)


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

    def _getJobDirectory(self, jobUid=None):
        """
        Return the absolute path to the CryoSPARC job directory on disk.
        """
        if jobUid is None:
            if hasattr(self, 'currenJob') and self.currenJob.get() is not None:
                jobUid = str(self.currenJob.get())
            else:
                return None

        if hasattr(self, 'projectDir') and self.projectDir.get():
            return os.path.join(self.projectDir.get(), str(jobUid))

        if hasattr(self, 'projectName') and self.projectName.get():
            try:
                projectDir = getProjectInformation(self.projectName.get(), info='project_dir')
                return os.path.join(projectDir, str(jobUid))
            except Exception:
                return None

        return None

    def _looksLikeFscTxt(self, filePath):
        """
        Heuristic check for CryoSPARC FSC raw-data txt files.
        """
        try:
            with open(filePath, "r", encoding="utf-8", errors="replace") as fh:
                firstLine = fh.readline().strip().lower()
        except Exception:
            return False

        if not firstLine:
            return False

        columns = [c.strip() for c in firstLine.split("\t")]

        if "wavenumber" not in columns:
            return False

        return any(
            c.startswith("fsc_") or
            "mask" in c or
            "resolution" in c or
            "noisesub" in c
            for c in columns
        )

    def _findFscTextFilesInJobDir(self, jobUid=None):
        """
        Return candidate FSC txt files found in the job directory, ordered by relevance.
        """
        jobDir = self._getJobDirectory(jobUid)
        if not jobDir or not os.path.isdir(jobDir):
            return []

        candidates = []

        for root, _, files in os.walk(jobDir):
            for name in files:
                lowerName = name.lower()
                if not lowerName.endswith(".txt"):
                    continue

                fullPath = os.path.join(root, name)

                # First filter by name hints
                nameScore = 0
                if "fsc" in lowerName:
                    nameScore += 100
                if "validation" in lowerName:
                    nameScore += 20
                if "plot" in lowerName:
                    nameScore += 10

                # Then verify by content
                if self._looksLikeFscTxt(fullPath):
                    try:
                        mtime = os.path.getmtime(fullPath)
                    except Exception:
                        mtime = 0
                    candidates.append((nameScore, mtime, fullPath))

        # Prefer higher score, then more recent file
        candidates.sort(key=lambda item: (item[0], item[1]), reverse=True)
        return [path for _, _, path in candidates]

    def _downloadLegacyFscFileById(self, fileId):
        """
        Legacy fallback: download the FSC text payload using fileid.
        Returns a local path or None.
        """
        if not fileId:
            return None

        system_info = getSystemInfo()
        status_errors = system_info[0]

        if status_errors:
            return None

        system_info = eval(system_info[1])
        cryosparcVersion = getCryosparcVersion()

        try:
            if parse_version(cryosparcVersion) < parse_version(V4_1_0):
                master_hostname = system_info.get('master_hostname')
                port_webapp = system_info.get('port_webapp')
                url = "http://%s:%s/file/%s" % (master_hostname, port_webapp, fileId)
                response = requests.get(url, allow_redirects=True)
            else:
                master_hostname = system_info.get('master_hostname')
                port_webapp = system_info.get('port_command_vis')
                url = "http://%s:%s/get_job_file" % (master_hostname, port_webapp)
                jsonParam = {'fileid': fileId}
                licence_id = _getLicenceFromFile()
                headers = {'License-ID': licence_id}
                response = requests.post(url, json=jsonParam, headers=headers, allow_redirects=True)

            response.raise_for_status()

            tmpDir = tempfile.mkdtemp(prefix="scipion_fsc_")
            fscFilePath = os.path.join(tmpDir, "fsc.txt")
            with open(fscFilePath, "wb") as fh:
                fh.write(response.content)

            if self._looksLikeFscTxt(fscFilePath):
                return fscFilePath

            return None
        except Exception:
            return None

    def _resolveFscTextFile(self, fileIdHint=None, jobUid=None):
        """
        Resolve the FSC text file with this order:
          1. filesystem lookup inside the job directory
          2. legacy HTTP download via fileid
        """
        # Prefer on-disk job outputs
        localCandidates = self._findFscTextFilesInJobDir(jobUid=jobUid)
        if localCandidates:
            return localCandidates[0]

        # Legacy fallback
        return self._downloadLegacyFscFileById(fileIdHint)

    def createFSC(self, idd, imgSet, vol):
        """
        Build a Scipion FSC set from a CryoSPARC FSC txt file.

        New strategy:
          1. Resolve the FSC txt file directly from the job directory on disk.
          2. If not found, fallback to legacy download through fileid.
        """
        jobUid = None
        if hasattr(self, 'currenJob') and self.currenJob.get() is not None:
            jobUid = str(self.currenJob.get())

        fscTxtPath = self._resolveFscTextFile(fileIdHint=idd, jobUid=jobUid)

        if not fscTxtPath or not os.path.exists(fscTxtPath):
            raise Exception(
                "Could not resolve the FSC txt file for job %s in project %s."
                % (jobUid, self.projectName.get() if hasattr(self, 'projectName') else "unknown")
            )

        factor = self._getInputParticles().getDim()[0] * imgSet.getSamplingRate()
        fscSet = self.getSetOfFCSsFromFile(fscTxtPath, factor)

        self._defineOutputs(outputFSC=fscSet)
        self._defineSourceRelation(vol, fscSet)

    def _splitFscLine(self, line):
        line = line.strip()
        if not line:
            return []

        if "\t" in line:
            return [token.strip() for token in line.split("\t") if token.strip()]

        return [token.strip() for token in re.split(r"\s+", line) if token.strip()]

    def _safeFloat(self, value):
        if value is None:
            return None

        text = str(value).strip()
        if text == "" or text.lower() in {"none", "nan", "null"}:
            return None

        try:
            return float(text)
        except Exception:
            return None

    def _canonicalizeFscColumnName(self, columnName):
        raw = str(columnName).strip()
        key = raw.lower()
        key = key.replace("%", "pct")
        key = re.sub(r"[^a-z0-9]+", "_", key)
        key = re.sub(r"_+", "_", key).strip("_")
        return key

    def _resolveFscColumn(self, columnName):
        """
        Map legacy and v5 CryoSPARC FSC column names to a stable internal key and
        a display label for Scipion.

        Returns: (internalKey, displayLabel) or (None, None) to skip a column.
        """
        key = self._canonicalizeFscColumnName(columnName)

        # Legacy
        if key in {"fsc_nomask", "nomask", "no_mask"}:
            return "nomask", "No mask"

        if key in {"fsc_sphericalmask", "sphericalmask", "spherical"}:
            return "spherical", "Spherical"

        if key in {"fsc_loosemask", "loosemask", "loose_mask"}:
            return "loose", "Loose"

        if key in {"fsc_tightmask", "tightmask", "tight_mask"}:
            return "tight", "Tight"

        if key in {"fsc_noisesub_raw", "noisesub_raw", "noise_sub_raw"}:
            return None, None

        if key in {"fsc_noisesub_true", "noisesub_true", "noise_sub_true"}:
            return "legacy_noise_sub_true", None

        if key in {"fsc_noisesub", "noisesub", "noise_sub", "corrected"}:
            return "corrected", "Corrected"

        # v5+: resolution mask
        if (
                "resolution" in key and
                "mask" in key and
                "auto" in key and
                "tight" in key and
                ("correct" in key or "noise" in key or "sub" in key)
        ):
            return "resolution_auto_tightened_corrected", "Auto-tightened corrected"

        if (
                "input" in key and
                "mask" in key and
                ("correct" in key or "noise" in key or "sub" in key)
        ):
            return "input_mask_corrected", "Input mask corrected"

        if (
                "resolution" in key and
                "mask" in key and
                "auto" in key and
                "tight" in key
        ):
            return "resolution_auto_tightened", "Auto-tightened resolution mask"

        if "input" in key and "mask" in key:
            return "input_mask", "Input mask"

        if "resolution" in key and "mask" in key:
            return "resolution_mask", "Resolution mask"

        # Generic corrected fallback for modern columns with unclear exact naming
        if "correct" in key or ("noise" in key and "sub" in key):
            return "corrected", "Corrected"

        # Unknown columns are skipped silently
        return None, None

    def _makeUniqueFscLabel(self, label, usedLabels, rawColumnName):
        if label not in usedLabels:
            usedLabels.add(label)
            return label

        newLabel = "%s (%s)" % (label, rawColumnName)
        usedLabels.add(newLabel)
        return newLabel

    def getSetOfFCSsFromFile(self, file, factor):
        with open(file, "r", encoding="utf-8", errors="replace") as f:
            lines = [line.rstrip("\n") for line in f if line.strip()]

        if not lines:
            raise Exception("Empty FSC file: %s" % file)

        header = self._splitFscLine(lines[0])
        if len(header) < 2:
            raise Exception("Invalid FSC header in file: %s" % file)

        columns = header[1:]
        fscSet = self._createSetOfFSCs()

        usedLabels = set()
        legacyTight = None
        legacyNoiseSubTrue = None
        explicitCorrectedPresent = False

        for colIndex, rawColumnName in enumerate(columns, start=1):
            internalKey, displayLabel = self._resolveFscColumn(rawColumnName)

            if internalKey is None:
                continue

            fsc = self.getFSCFromRawData(
                lines=lines,
                col=colIndex,
                factor=factor,
                objLabel=displayLabel or rawColumnName
            )

            # Skip empty/broken curves
            xData, yData = fsc.getData()
            if not xData or not yData:
                continue

            if internalKey == "legacy_noise_sub_true":
                legacyNoiseSubTrue = yData
                continue

            if internalKey == "tight":
                legacyTight = yData

            if internalKey in {
                "corrected",
                "resolution_auto_tightened_corrected",
                "input_mask_corrected",
            }:
                explicitCorrectedPresent = True

            finalLabel = self._makeUniqueFscLabel(
                displayLabel or rawColumnName,
                usedLabels,
                rawColumnName
            )
            fsc.setObjLabel(finalLabel)
            fscSet.append(fsc)

        # Legacy fallback: reconstruct "Corrected" if old files provide
        # tight + noisesub_true but no explicit corrected curve.
        if (
                not explicitCorrectedPresent and
                legacyTight is not None and
                legacyNoiseSubTrue is not None and
                len(legacyTight) == len(legacyNoiseSubTrue)
        ):
            corrected = []

            for tightVal, noiseVal in zip(legacyTight, legacyNoiseSubTrue):
                if noiseVal is None or tightVal is None:
                    corrected.append(None)
                    continue

                denom = 1.0 - noiseVal

                # Keep behavior numerically safe near 1.0
                if abs(denom) < 1e-6:
                    denom = 1e-6

                corrected.append((tightVal - noiseVal) / denom)

            firstItem = fscSet.getFirstItem()
            if firstItem is not None:
                xAxis = list(firstItem.getData()[0])
                yAxis = [v for v in corrected]

                correctedFsc = FSC(objLabel=self._makeUniqueFscLabel(
                    "Corrected",
                    usedLabels,
                    "legacy_corrected"
                ))
                correctedFsc.setData(xAxis, yAxis)
                fscSet.append(correctedFsc)

        fscSet.write()
        return fscSet

    def getFSCFromRawData(self, lines, col, factor, objLabel="FSC"):
        xAxis = []
        yAxis = []

        for line in lines[1:]:
            parts = self._splitFscLine(line)
            if len(parts) <= col:
                continue

            waveNumber = self._safeFloat(parts[0])
            corrValue = self._safeFloat(parts[col])

            if waveNumber is None or corrValue is None:
                continue

            xAxis.append(waveNumber / factor)
            yAxis.append(corrValue)

        fsc = FSC(objLabel=objLabel)
        fsc.setData(xAxis, yAxis)
        return fsc

    def findLastIteration(self, jobName):
        get_job_streamlog(
            self.projectName.get(),
            jobName,
            self._getFileName('stream_log')
        )

        with open(self._getFileName('stream_log'), encoding="utf-8", errors="replace") as f:
            rawText = f.read()

        parsedEvents = None

        # First, try the old structured representation.
        try:
            parsed = ast.literal_eval(rawText.strip())
            if isinstance(parsed, dict):
                parsedEvents = [parsed]
            elif isinstance(parsed, list):
                parsedEvents = parsed
        except Exception:
            parsedEvents = None

        idd = None
        itera = None

        if parsedEvents is not None:
            for event in parsedEvents:
                if not isinstance(event, dict):
                    continue

                text = str(event.get("text", ""))

                if text.startswith("FSC Iteration") or text.startswith("FSC iIteration"):
                    match = re.search(r"FSC i?Iteration\s+(\d+)", text)
                    if match:
                        itera = match.group(1)

                    for imgfile in event.get("imgfiles", []) or []:
                        if isinstance(imgfile, dict) and imgfile.get("filetype") == "txt":
                            idd = imgfile.get("fileid")
                            break

                elif "Using Filter Radius" in text:
                    try:
                        nomRes = text.split("(")[1].split(")")[0].replace("A", "Å")
                        self.mapResolution = pwobj.String(nomRes)
                        self._store(self)
                    except Exception:
                        pass

                elif "Estimated Bfactor" in text:
                    try:
                        estBFactor = text.split(":", 1)[1].strip()
                        self.estBFactor = pwobj.String(estBFactor)
                        self._store(self)
                    except Exception:
                        pass

            return idd, itera

        # Fallback for plain-text v5 event output.
        for line in rawText.splitlines():
            if "FSC Iteration" in line or "FSC iIteration" in line:
                match = re.search(r"FSC i?Iteration\s+(\d+)", line)
                if match:
                    itera = match.group(1)

            if "Using Filter Radius" in line:
                match = re.search(r"\(([^)]+)\)", line)
                if match:
                    self.mapResolution = pwobj.String(match.group(1).replace("A", "Å"))
                    self._store(self)

            if "Estimated Bfactor" in line:
                match = re.search(r"Estimated Bfactor\s*:?\s*(.+)$", line)
                if match:
                    self.estBFactor = pwobj.String(match.group(1).strip())
                    self._store(self)

            # Best-effort file id extraction if the raw output still contains it.
            if idd is None:
                match = re.search(r"fileid['\"]?\s*[:=]\s*['\"]([^'\"]+)['\"]", line)
                if match:
                    idd = match.group(1)

        return idd, itera

    def _createModelFile(self):
        pass

    def getLogLine(self):
        return self._logLastLine

    def setLogLine(self, lastLine: int):
        self._logLastLine = lastLine
