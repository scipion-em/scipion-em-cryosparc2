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
import os.path
import zipfile

from pyworkflow import BETA
from pyworkflow.protocol.params import (PointerParam, FloatParam,
                                        LEVEL_ADVANCED, EnumParam, IntParam,
                                        Positive, BooleanParam)
from pwem.protocols import ProtRefine3D

from . import ProtCryosparcBase
from ..convert import *
from ..utils import *
from ..constants import *


class ProtCryoSparc3DVariabilityAnalisys(ProtCryosparcBase, ProtRefine3D):
    """
    Protocol to create various versions of a 3D variability result that can be
    used for display
    """
    _label = '3D variability Display'
    _devStatus = BETA

    def _defineFileNames(self):
        """ Centralize how files are called. """
        myDict = {
            'out_particles': self._getExtraPath('output_particle.star'),
            'stream_log': self._getPath() + '/stream.log',
            'out_class': self._getExtraPath() + '/output_class.star'
        }
        self._updateFilenamesDict(myDict)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('input3DVariablityAnalisysProt', PointerParam,
                      label="3D variability analysis protocol",
                      pointerClass='ProtCryoSparc3DVariability',
                      help='Particle stacks to use.')
        form.addParallelSection(threads=1, mpi=1)

        # --------------[3D Variability Display]---------------------------
        form.addSection(label='3D Variability Output')

        form.addParam('var_output_mode', EnumParam,
                      choices=['cluster', 'simple', 'intermediates'],
                      default=0,
                      label="Output mode",
                      help='simple mode: output a simple linear "movie" of '
                           'volumes along each dimension. Number of frames is '
                           'the next param. '
                           'cluster mode: fit clusters to reaction coordinates '
                           'and output volumes and particles from each cluster. '
                           'Number of clusters is the next param. '
                           'intermediates mode: - reconstruct multiple '
                           'intermediate volumes along each variability '
                           'dimension, useful for better visualizing '
                           'non-linear changes. Number of intermediate frames '
                           'is the next param.')
        form.addParam('var_num_frames', IntParam, default=20,
                      validators=[Positive],
                      label="Number of frames/clusters",
                      help='Number of clusters for cluster mode, or number of '
                           'frames for simple and intermediates mode.')

        form.addParam('var_range_percentile', FloatParam, default=3,
                      validators=[Positive],
                      label="Min/Max Range Percentile (%)",
                      help='Percentile for computing the min and max for '
                           'each component.')

        form.addParam('var_N', IntParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Downsample to box size",
                      help='Downsample the output volumes to this size.')

        form.addParam('var_M', IntParam, default=None,
                      expertLevel=LEVEL_ADVANCED,
                      allowsNull=True,
                      label="Crop to size (after downsample)",
                      help='Crop the output volumes to this size after '
                           'downsampling.')

        # form.addParam('var_num_particles', IntParam, default=None,
        #               validators=[Positive],
        #               label="Only use this many particles",
        #               help='Only use this many particles for reconstructions.')

        form.addParam('var_filter_res', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Filter resolution (A)",
                      help='Resolution at which results are filtered')

        form.addParam('var_filter_order', FloatParam, default=1.5,
                      validators=[Positive],
                      label="Filter order",
                      help='Order of filter')

        form.addParam('var_highpass_res', FloatParam, default=None,
                      allowsNull=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Highpass resolution (A)",
                      help='Resolution at which results are filtered')

        form.addParam('var_highpass_order', FloatParam, default=8,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="Highpass order",
                      help='Order of filter')

        form.addParam('var_vol_flip', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Flip hand of outputs",
                      help='Flip the hand out output volumes (including series). '
                           'Note that this does not flip output '
                           'alignments - only the output volumes.')

        form.addParam('var_skip_reconstruction', BooleanParam, default=False,
                      expertLevel=LEVEL_ADVANCED,
                      label="Skip reconstruction",
                      help='Skip reconstructions in cluster and '
                           'intermediate mode - helpful to quickly inspect '
                           'whether a clustering or division of particles will'
                           ' be helpful.')

        form.addParam('var_cluster_3D_plots', BooleanParam, default=True,
                      expertLevel=LEVEL_ADVANCED,
                      label="Cluster mode: 3D plots",
                      help='Make 3D instead of 2D plots to show clusters.')

        form.addParam('var_intermediate_width', IntParam, default=2,
                      validators=[Positive],
                      expertLevel=LEVEL_ADVANCED,
                      label="Intermediates: window (frames)",
                      help='Intermediates: window (frames)')
        form.addSection(label="Compute settings")
        addComputeSectionParams(form, allowMultipleGPUs=False)

    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._defineFileNames()
        self._defineParamsName()
        self._initializeCryosparcProject()
        self._insertFunctionStep(self.processStep)
        self._insertFunctionStep(self.generateStartFiles)
        self._insertFunctionStep(self.createOutputStep)

    # --------------------------- STEPS functions ------------------------------
    def processStep(self):
        print(pwutils.yellowStr("3D Variability started..."), flush=True)
        self.doRun3DVariabilityDisplay()

    def generateStartFiles(self):
        print(pwutils.yellowStr("Copying files from CS to Scipion folder..."),
              flush=True)
        csOutputFolder = os.path.join(self.projectDir.get(),
                                      self.run3DVariabilityDisplay.get())
        import time
        time.sleep(10)
        if self.var_output_mode.get() == 0:  # Cluster mode

            # Copy the CS output to extra folder
            outputFolder = os.path.join(self._getExtraPath(), 'outputs')
            pwutils.cleanPath(outputFolder)
            copyFiles(csOutputFolder, outputFolder)

            # Creating the folder where .star files will be create
            starFilesPath = os.path.join(outputFolder, 'clustersFiles')
            os.makedirs(starFilesPath)

            # Generating the .star file por all clusters
            allClusterSeries = '%s%s_series_all_clusters.cs' % (getOutputPreffix(self.projectName.get()),
                                                             self.run3DVariabilityDisplay.get())
            allClusterSeriesStar = '%s%s_series_all_clusters.star' % (getOutputPreffix(self.projectName.get()),
                                                             self.run3DVariabilityDisplay.get())

            filePath = os.path.join(outputFolder, allClusterSeries)
            outputStarFn = os.path.join(starFilesPath, allClusterSeriesStar)
            argsList = [filePath, outputStarFn]
            convertCs2Star(argsList)

            # Creating the .star cluster per cluster
            numOfClusters = self.var_num_frames.get()
            print(pwutils.yellowStr("Generating .star files for all clusters"), flush=True)
            for i in range(numOfClusters):
                print(pwutils.yellowStr("Processing cluster %d ..." % i), flush=True)

                clusterParticlesPattern = '%s%s_cluster_%03d_particles.cs' % (getOutputPreffix(self.projectName.get()),
                                                                              self.run3DVariabilityDisplay.get(), i)

                passthroughParticlesPattern = '%s%s_passthrough_particles_cluster_%d.cs' % (getOutputPreffix(self.projectName.get()),
                                                                                            self.run3DVariabilityDisplay.get(), i)
                patterns = [clusterParticlesPattern, passthroughParticlesPattern]
                outputs = ['output_cluster_particle%03d.star' % i,
                           'output_passthrough_particle%03d.star' % i]
                # Generating .star files
                for j in range(len(patterns)):
                    filePath = os.path.join(outputFolder, patterns[j])
                    outputStarFn = os.path.join(starFilesPath, outputs[j])
                    argsList = [filePath, outputStarFn]
                    convertCs2Star(argsList)
        else:
            # Copy the CS output to extra folder
            outputFolder = self._getExtraPath()
            pwutils.cleanPath(outputFolder)
            copyFiles(csOutputFolder, outputFolder)
            # Creating the .star cluster particles
            numOfComponets = int(self.input3DVariablityAnalisysProt.get().var_K)

            for i in range(numOfComponets):
                print(pwutils.yellowStr("Extracting component %d ..." % i),
                      flush=True)
                # Creating the folder where clusters will be unzip
                componetFilesPath = self._getExtraPath('component%03d' % i)
                os.makedirs(componetFilesPath)
                componentPattern = "%s%s_component_%03d.zip" % (
                    getOutputPreffix(self.projectName.get()),
                    self.run3DVariabilityDisplay.get(),
                    i)
                cmd = "unzip %s -d %s" % (
                self._getExtraPath(componentPattern),
                componetFilesPath)
                os.system(cmd)

    def createOutputStep(self):
        self._initializeUtilsVariables()
        print(pwutils.yellowStr("Creating the output..."), flush=True)

        pass

    # --------------------------- UTILS functions ------------------------------

    def _defineParamsName(self):
        """ Define a list with all protocol parameters names"""
        self._paramsName = ['var_output_mode',
                            'var_num_frames',
                            'var_range_percentile',
                            'var_N',
                            'var_M',
                            'var_filter_res',
                            'var_filter_order',
                            'var_highpass_res',
                            'var_highpass_order',
                            'var_vol_flip',
                            'var_skip_reconstruction',
                            'var_cluster_3D_plots',
                            'var_intermediate_width']
        self.lane = str(self.getAttributeValue('compute_lane'))

    def _validate(self):
        validateMsgs = cryosparcValidate()
        if not validateMsgs:
            validateMsgs = gpusValidate(self.getGpuList(), checkSingleGPU=True)
        return validateMsgs

    def doRun3DVariabilityDisplay(self):
        """
        :return:
        """
        className = "var_3D_disp"
        varAnalysisJob = str(self.input3DVariablityAnalisysProt.get().run3DVariability)
        input_group_conect = {"particles": str('%s.particles' %varAnalysisJob) ,
                              "volume": str('%s.volume'%varAnalysisJob)}
        params = {}

        for paramName in self._paramsName:
            if (paramName != 'var_output_mode' and paramName != 'var_N' and
                    paramName != 'var_M' and paramName != 'var_filter_res' and
                    paramName != 'var_highpass_res'):
                params[str(paramName)] = str(self.getAttributeValue(paramName))
            elif paramName == 'var_output_mode':
                params[str(paramName)] = str(VAR_OUTPUT_MODE[self.var_output_mode.get()])
            elif paramName == 'var_N' and self.var_N.get() is not None:
                params[str(paramName)] = str(self.var_N.get())
            elif paramName == 'var_M' and self.var_M.get() is not None:
                params[str(paramName)] = str(self.var_M.get())
            elif paramName == 'var_filter_res' and self.var_filter_res.get() is not None:
                params[str(paramName)] = str(self.var_filter_res.get())
            elif paramName == 'var_highpass_res' and self.var_highpass_res.get() is not None:
                params[str(paramName)] = str(self.var_highpass_res.get())

        # Determinate the GPUs to use (in dependence of the cryosparc version)
        try:
            gpusToUse = self.getGpuList()
        except Exception:
            gpusToUse = False

        self.run3DVariabilityDisplay = enqueueJob(className,
                                                  self.projectName.get(),
                                                  self.workSpaceName.get(),
                                                  str(params).replace('\'',
                                                                      '"'),
                                                  str(input_group_conect).replace(
                                                      '\'',
                                                      '"'),
                                                  self.lane, gpusToUse)

        self.currenJob.set(self.run3DVariabilityDisplay.get())
        self._store(self)

        waitForCryosparc(self.projectName.get(),
                         self.run3DVariabilityDisplay.get(),
                         "An error occurred in the 3D Variability process. "
                         "Please, go to cryosPARC software for more "
                         "details.")
