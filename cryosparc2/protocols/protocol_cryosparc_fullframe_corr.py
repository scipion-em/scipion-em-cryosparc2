# **************************************************************************
# *
# *  Authors:     Szu-Chi Chung (phonchi@stat.sinica.edu.tw)
# *
# * SABID Laboratory, Institute of Statistical Science, Academia Sinica
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

import sys
import os
import commands
import numpy as np
from math import ceil
from itertools import izip

import pyworkflow.em as em
import pyworkflow.em.metadata as md
import pyworkflow as pw

import pyworkflow.protocol.constants as cons
from pyworkflow.protocol.params import (BooleanParam, IntParam, Positive, StringParam, EnumParam)
from pyworkflow.em.data import String
from pyworkflow.em.protocol import ProtAlignMovies
from pyworkflow.gui.plotter import Plotter
from pyworkflow.protocol import STEPS_SERIAL

from cryosparc2.convert import *
from cryosparc2.utils import *
from cryosparc2.constants import *

relionPlugin = pwutils.importFromPlugin("relion.convert", doRaise=True)


class ProtCryoFullFrame(ProtAlignMovies):
    """ Wrapper to CryoSparc Full-frame motion correction.
        Estimate and correct for full-frame motion (e.g., stage drift) from movie data.
    """
    _label = 'perform Full-frame motion correction'
    
    def __init__(self, **kwargs):
        ProtAlignMovies.__init__(self, **kwargs)
        self.stepsExecutionMode = STEPS_SERIAL
        if self.numberOfMpi.get() < 2:
            self.numberOfMpi.set(2)

    def _getConvertExtension(self, filename):
        """ Check whether it is needed to convert to .mrc or not """
        ext = pwutils.getExt(filename).lower()
        return None if ext in ['.mrc', '.mrcs', '.tiff', '.tif'] else 'mrc'

    # --------------------------- DEFINE param functions -----------------------
    def _defineAlignmentParams(self, form):
        
        form.addParam('smooth_lambda_cal', IntParam, default=2,
                      validators=[Positive],
                      label='Calibrated smoothing',
                      help='Calibrated smoothing constant applied to trajectories.')
        form.addParam('res_max_align', IntParam, default=5,
                      validators=[Positive],
                      label='Maximum alignment resolution (A)',
                      help='Maximum resolution (in A) to consider when aligning frames.'
                           'Generally, betwen 5A and 3A is best.')
        form.addParam('bfactor', IntParam, default=500,
                      validators=[Positive],
                      label='B-factor during alignment',
                      help='B-factor that blurs frames before aligning.'
                           'Generally 500 to 100 is best.')
        form.addParam('integer_shifts', BooleanParam, default=False,
                      label='Force integer shifts',
                      help='Force shifts to be an integer number of pixels.'
                           'This ensures that during frame averaging, no sub-pixel shifts are used.'
                           'Sub-pixel shifts can sometimes cause aliasing artefacts.')
        line = form.addLine('Frames for corrected SUM',
                             help='First and last frames to use in corrected '
                                  'average (starts counting at 1 and 0 as last '
                                  'means util the last frame in the movie). ')
        line.addParam('frame_start', IntParam, default=1,
                      label='from')
        line.addParam('frame_end', IntParam, default=0,
                      label='to')
        form.addParam('zero_shift_frame', IntParam, default=-1, expertLevel=cons.LEVEL_ADVANCED,
                       label='Zero-shift frame',
                       help='Default to central frame. This frame number is relative to the start of the movie,'
                            'not the frame_start parameter above.')
        form.addParam('do_plots', BooleanParam, default=True, expertLevel=cons.LEVEL_ADVANCED,
                       label='Make motion diagnostic plots',
                       help="Whether or not to make plots of motion trajectories.")
        form.addParam('num_plots', IntParam, default=10, expertLevel=cons.LEVEL_ADVANCED,
                       label='Number of movies to plot',
                       help=" Only make plots for the first this many movies.")

        form.addParam('random_num', IntParam, default=-1, expertLevel=cons.LEVEL_ADVANCED,
                       label='Only process this many movies',
                       help="Randomly select this many movies to process."
                            " -1 means use all movies. "
                            "Helpful for tweaking params.")

        form.addParam('gainref_rotate_num', EnumParam, default=0, expertLevel=cons.LEVEL_ADVANCED,
                       choices=['No rotation (0)',
                                ' 90 degrees (1)',
                                '180 degrees (2)',
                                '270 degrees (3)'],
                       label='Gain rotation',
                       help="Rotate gain ref counter-clockwise by 90 degrees this many times.")

        form.addParam('gainref_flip_x', BooleanParam, default=False, expertLevel=cons.LEVEL_ADVANCED,
                       label='Gain flip',
                       help="Flip gain ref left-to-right (in X axis).")
        form.addParam('gainref_flip_y', BooleanParam, default=False, expertLevel=cons.LEVEL_ADVANCED,
                       label='Gain flip',
                       help="Flip gain ref top-to-bottom (in Y axis).")
        form.addParallelSection(threads=1, mpi=1)
        # ----------- [Compute settings] --------------------------------
        form.addSection(label="Compute settings")
        form.addParam('numberGPU', IntParam, default=1,
                      validators=[Positive],
                      label='Number of GPUs to parallelize',
                      help='Number of GPUs to use during motion correction.')



    # --------------------------- INSERT steps functions -----------------------
    def _insertAllSteps(self):
        self._initializeCryosparcProject()
        self._insertFunctionStep("convertInputAndProcessStep")
        self._insertFunctionStep('createOutputStep')

    # --------------------------- STEPS functions ------------------------------
    def convertInputAndProcessStep(self):
        print(pwutils.greenStr("Importing Movies..."))
        inputMovies = self.inputMovies.get()
        pw.utils.makePath(self._getExtraPath("process"))
        relionPlugin.writeSetOfMovies(
            inputMovies, self._getPath("movies.star"),
            postprocessImageRow=self._postprocessImageRow)
        self._importMovies()
        while getJobStatus(self.projectName.get(), self.importedMovies.get()) not in STOP_STATUSES:
            waitJob(self.projectName.get(), self.importedMovies.get())

        print(pwutils.greenStr("Motion Correction Started..."))
        print(str(self.mov))
        self.runFullFrame = String(self.doRunFullFrame()[-1].split()[-1])
        self.currenJob.set(self.runFullFrame.get())
        #self._store(self)
        while getJobStatus(self.projectName.get(), self.runFullFrame.get()) not in STOP_STATUSES:
            waitJob(self.projectName.get(), self.runFullFrame.get())

        self._initializeUtilsVariables()
        print (pwutils.greenStr("Creating the output..."))

        # Link the folder on SSD to scipion directory
        os.system("ln -s " + self.projectPath + "/" + self.projectName.get() + '/' +
                  self.runFullFrame.get() + " " + self._getExtraPath())
        self.cryolink = self.runFullFrame.get() + "/motioncorrected/" 
        self._loadInputList()
        for idx, movie in enumerate(self.listOfMovies):
            # Mark this movie as finished
            if self.frame_end.get() == 0:
                end = self.inputMovies.get().getDim()[2]
            else:
                end = self.frame_end.get()
            if idx < end and idx >= self.frame_start.get() -1:
                self._moveFiles(movie)
                self._computeExtra(movie)
                movieDoneFn = self._getMovieDone(movie)
                open(movieDoneFn, 'w').close() # Trigger the createoutput!

    def createOutputStep(self):
        # validate that we have some output movies
        output = self.outputMicrographs

        if output.getSize() == 0 and len(self.listOfMovies) != 0:
            raise Exception(redStr("All movies failed, didn't create outputMicrographs."
                                   "Please review movie processing steps above."))
        elif output.getSize() < len(self.listOfMovies):
            self.warning(pwutils.yellowStr("WARNING - Failed to align %d movies."
                                   % (len(self.listOfMovies) - self.outputMicrographs.getSize())))
   
    def setAborted(self):
        """ Set the status to aborted and updated the endTime. """
        ProtCryoFullFrame.setAborted(self)
        killJob(str(self.projectName.get()), str(self.currenJob.get()))

    # --------------------------- INFO functions -------------------------------
    def _summary(self):
        summary = []
        return summary

    def _validate(self):
        # Check base validation before the specific ones for Motioncor
        errors = ProtAlignMovies._validate(self)
        return errors
    
    # --------------------------- UTILS functions ------------------------------
    def _getInputParticles(self):
        return self.inputParticles.get()

    def _initializeUtilsVariables(self):
        """
        Initialize all utils cryoSPARC variables
        """
        self._program = getCryosparcProgram()
        self._user = getCryosparcUser()
        self._ssd = getCryosparcSSD()

        # Create a cryoSPARC project dir
        self.projectDirName = suffix + self.getProject().getShortName()
        self.projectPath = pwutils.join(self._ssd, self.projectDirName)
        self.projectDir = createProjectDir(self.projectPath)

    def _initializeCryosparcProject(self):
        """
        Initialize the cryoSPARC project and workspace
        """
        self._initializeUtilsVariables()
        # create empty project or load an exists one
        folderPaths = getProjectPath(self.projectPath)
        if not folderPaths:
            self.a = createEmptyProject(self.projectPath, self.projectDirName)
            self.projectName = self.a[-1].split()[-1]
        else:
            self.projectName = str(folderPaths[0])

        self.projectName = String(self.projectName)
        self._store(self)

        # create empty workspace
        self.b = createEmptyWorkSpace(self.projectName, self.getRunName(),
                                      self.getObjComment())
        self.workSpaceName = String(self.b[-1].split()[-1])
        self._store(self)

    def _importMovies(self):

        print("Importing Movies")

        # import_particles_star
        self.c = self.doImportMovies()

        self.importedMovies = String(self.c[-1].split()[-1])
        self._store(self)

        self.currenJob = String(self.importedMovies.get())
        self._store(self)

        self.mov = String(self.importedMovies.get() + '.imported_movies')

    def doImportMovies(self):
        """
        do_import_particles_star(puid, wuid, uuid, abs_star_path,
                                 abs_blob_path=None, psize_A=None)
        returns the new uid of the job that was created
        """
        className = "import_movies"
        # {'particles' : 'JXX.imported_particles' }
        preExp, dose = self._getCorrectedDose(self.inputMovies.get())
        params = {"blob_paths": str(os.getcwd()+"/"+self._getExtraPath("process")+"/*"),
                  "gainref_path": str(self.inputMovies.get().getGain()),
                  "psize_A": str(self.inputMovies.get().getSamplingRate()),
                  "accel_kv": str(self.inputMovies.get().getAcquisition().getVoltage()),
                  "cs_mm": str(self.inputMovies.get().getAcquisition().getSphericalAberration()),
                  "total_dose_e_per_A2": str(dose*self.inputMovies.get().getDim()[2]),
                  "gainref_flip_x": str(self.gainref_flip_x.get()),
                  "gainref_flip_y": str(self.gainref_flip_y.get()),
                  "gainref_rotate_num": str(self.gainref_rotate_num.get())}

        return doJob(className, self.projectName.get(), self.workSpaceName.get(),
                     str(params).replace('\'', '"'), {})

    def doRunFullFrame(self):
        """
        doRunFullFrame:  do_job(job_type, puid='P1', wuid='W1',
                                 uuid='devuser', params={},
                                 input_group_connects={})
        returns: the new uid of the job that was created
        """
        if self.numberGPU.get() > 1:
            className = "rigid_motion_correction_multi"
        else:
            className = "rigid_motion_correction"
        
        input_group_conect = {"movies": str(self.mov)}
        # {'movies' : 'JXX.imported_movies' }

        params = {"smooth_lambda_cal": str(self.smooth_lambda_cal.get()),
                  "res_max_align": str(self.res_max_align.get()),
                  "bfactor": str(self.bfactor.get()),
                  "integer_shifts": str(self.integer_shifts.get()),
                  "frame_start": str(self.frame_start.get()-1),
                  "num_plots": str(self.num_plots.get()),
                  "do_plots": str(self.do_plots.get()),
                  "compute_num_gpus": str(self.numberGPU.get())}

        if self.frame_end.get() != 0:
            params['frame_end'] = self.frame_end.get()
        if self.zero_shift_frame.get() != -1:
            params['zero'] = self.zero_shift_frame.get()
        if self.random_num.get() != -1:
            params['random_num'] = self.random_num.get()

        return doJob(className, self.projectName.get(), self.workSpaceName.get(),
                     str(params).replace('\'', '"'),
                     str(input_group_conect).replace('\'', '"'))

    def _postprocessImageRow(self, img, row):
        """ Simple way to link movie to extra. """
        # link movie to extra!!!!!!!
        rlnMicrographMovieName = row.getValue('rlnMicrographMovieName')
        linkName = self._getExtraPath("process/"+rlnMicrographMovieName.split('/')[-1])
        os.system("ln -s "+os.getcwd()+'/'+rlnMicrographMovieName+' '+linkName)
        row.setValue('rlnMicrographMovieName', linkName)

    def _getMovieOutFn(self, movie, suffix):
        movieBase = pwutils.removeBaseExt(movie.getFileName())#.replace('.', '_')
        return os.path.join(self._getExtraPath(), self.cryolink ,
                            '%s%s' % (movieBase, suffix))

    def _getMovieExtraFn(self, movie, suffix):
        """ Return filenames in the extra directory with the prefix of this movie.
        Used to keep files associated with each micrograph.
        """
        movieBase = pwutils.removeBaseExt(movie.getFileName())
        return self._getExtraPath('%s%s' % (movieBase, suffix))

    def _moveFiles(self, movie):
        # It really annoying that Relion default names changes if you use DW or not
        # if use DW, the default name are DW and the others noDW
        pwutils.createLink(self._getMovieOutFn(movie, '_rigid_aligned.mrc'),
                             self._getExtraPath(self._getOutputMicName(movie)))
        pwutils.createLink(self._getMovieOutFn(movie, '_traj.npy'), self._getMovieExtraFn(movie, '_traj.npy'))
    def _stepsCheck(self):
        #self._checkNewInput()
        # we will handle check input here for batch processing, the streaming functionality may be more adequate using cryosparc live
        self._checkNewOutput()

    def _getMovieDone(self, movie):
        return self._getExtraPath('DONE_movie_%06d.TXT' % movie.getObjId())

    def _computeExtra(self, movie):
        """ Compute thumbnail, PSD and plots. """
        #inputMovies = self.inputMovies.get()
        #movieFolder = self._getExtraPath(self.cryolink)
        #outMicFn = self._getMovieOutFn(movie, '_rigid_aligned.mrc')

        #if self.doComputeMicThumbnail:
        #    self.computeThumbnail(outMicFn,
        #                          outputFn=self._getOutputMicThumbnail(movie))
        #
        #if self.doComputePSD:
        #    #fakeShiftsFn = self.writeZeroShifts(movie)
        #    movieFn = movie.getFileName()
        #    aveMicFn = os.path.join(movieFolder,
        #                            pwutils.removeBaseExt(movieFn) + "_tmp.mrc")
        #    self.averageMovie(movie, movieFn, aveMicFn,
        #                      binFactor=self.binFactor.get(),
        #                      dark=inputMovies.getDark(),
        #                      gain=inputMovies.getGain())
        #
        #    self.computePSDs(movie, aveMicFn, outMicFn,
        #                     outputFnCorrected=self._getPsdJpeg(movie))

        self._saveAlignmentPlots(movie)

    def _saveAlignmentPlots(self, movie):
        # Create plots and save as an image
        shiftsX, shiftsY = self._getMovieShifts(movie)
        first, _ = self._getFrameRange(movie.getNumberOfFrames(), 'align')
        plotter = createGlobalAlignmentPlot(shiftsX, shiftsY, first)
        plotter.savefig(self._getPlotGlobal(movie))
        plotter.close()

    def _getPlotGlobal(self, movie):
        return self._getNameExt(movie, '_global_shifts', 'png', extra=True)

    def _getNameExt(self, movie, postFix, ext, extra=False):
        fn = self._getMovieRoot(movie) + postFix + '.' + ext
        return self._getExtraPath(fn) if extra else fn

    def _createOutputMicrographs(self):
        # To create the unweighted average micrographs
        # we only consider the 'doSaveUnweightedMic' flag if the
        # weighted ones should be created.
        return True

    def _createOutputWeightedMicrographs(self):
        return False

    def _createOutputMovies(self):
        return False

    def _getRange(self, movie):
        n = self._getNumberOfFrames(movie)
        iniFrame, _, indxFrame = movie.getFramesRange()
        first, last = self._getFrameRange(n, 'sum')

        if iniFrame != indxFrame:
            first -= iniFrame
            last -= iniFrame

        return first, last

    def _getNumberOfFrames(self, movie):
        _, lstFrame, _ = movie.getFramesRange()

        if movie.hasAlignment():
            _, lastFrmAligned = movie.getAlignment().getRange()
            if lastFrmAligned != lstFrame:
                return lastFrmAligned
        return movie.getNumberOfFrames()

    def _getMovieShifts(self, movie):
        outnpy = self._getMovieExtraFn(movie, '_traj.npy')
        first, last = self._getRange(movie)
        n = last - first + 1
        shift_array = np.load(outnpy)
        shift_array = shift_array[0]
        xShifts, yShifts = [], []

        for j in xrange(shift_array.shape[0]):
            xShifts.append(float(shift_array[j][0]))
            yShifts.append(float(shift_array[j][1]))
            if len(xShifts) == n:
                break

        return xShifts, yShifts

    def _preprocessOutputMicrograph(self, mic, movie):
        self._setPlotInfo(movie, mic)

    def _setPlotInfo(self, movie, mic):
        mic.plotGlobal = em.Image(location=self._getPlotGlobal(movie))
        #if self.doComputePSD:
        #    mic.psdCorr = em.Image(location=self._getPsdCorr(movie))
        #    mic.psdJpeg = em.Image(location=self._getPsdJpeg(movie))
        #if self.doComputeMicThumbnail:
        #    mic.thumbnail = em.Image(
        #        location=self._getOutputMicThumbnail(movie))

def createGlobalAlignmentPlot(meanX, meanY, first):
    """ Create a plotter with the shift per frame. """
    figureSize = (6, 4)
    plotter = Plotter(*figureSize)
    figure = plotter.getFigure()
    ax = figure.add_subplot(111)
    ax.grid()
    ax.set_title('Global shift')
    ax.set_xlabel('Shift x (pixels)')
    ax.set_ylabel('Shift y (pixels)')

    i = first
    skipLabels = ceil(len(meanX)/10.0)
    labelTick = 1

    for x, y in izip(meanX, meanY):
        if labelTick == 1:
            ax.text(x - 0.02, y + 0.02, str(i))
            labelTick = skipLabels
        else:
            labelTick -= 1
        i += 1

    ax.plot(meanX, meanY, color='b')
    ax.plot(meanX, meanY, 'yo')

    plotter.tightLayout()

    return plotter

