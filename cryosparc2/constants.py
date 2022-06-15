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
import enum
from pwem.constants import (
    SYM_CYCLIC, SYM_TETRAHEDRAL, SYM_OCTAHEDRAL, SYM_I222,
    SYM_I222r)


from pwem.constants import SYM_DIHEDRAL_Y
# Root folder where cryosparc is installed, we will look here for the client
CRYOSPARC_HOME = 'CRYOSPARC_HOME'
CRYOSPARC_DIR = 'CRYOSPARC_DIR'  # Legacy, replaced by CRYOSPARC_HOME
# Optional: Folder where cryosparc projects will be created
CRYO_PROJECTS_DIR = 'CRYO_PROJECTS_DIR'
CRYOSPARC_USER = 'CRYOSPARC_USER'
CRYOSPARC_USE_SSD = 'CRYOSPARC_USE_SSD'
CRYOSPARC_MASTER = 'cryosparc_master'
CRYOSPARC_STANDALONE_INSTALLATION = 'CRYOSPARC_STANDALONE_INSTALLATION'
CRYOSPARC_DEFAULT_LANE = 'CRYOSPARC_DEFAULT_LANE'

# Supported versions:
V_UNKNOWN ='v0.0.0'
V2_5_0 = 'v2.5.0'
V2_8_0 = 'v2.8.0'
V2_9_0 = 'v2.9.0'
V2_11_0 = 'v2.11.0'
V2_12_0 = 'v2.12.0'
V2_12_2 = 'v2.12.2'
V2_12_4 = 'v2.12.4'
V2_13_0 = 'v2.13.0'
V2_13_2 = 'v2.13.2'
V2_14_0 = 'v2.14.0'
V2_14_2 = 'v2.14.2'
V2_15_0 = 'v2.15.0'
V3_0_0 = 'v3.0.0'
V3_0_1 = 'v3.0.1'
V3_1_0 = 'v3.1.0'
V3_2_0 = 'v3.2.0'
V3_3_0 = 'v3.3.0'
V3_3_1 = 'v3.3.1'

# Symmetry dict
CS_SYM_NAME = dict()
CS_SYM_NAME[SYM_CYCLIC] = 'Cn'
CS_SYM_NAME[SYM_DIHEDRAL_Y] = 'Dn'
CS_SYM_NAME[SYM_TETRAHEDRAL] = 'T'
CS_SYM_NAME[SYM_OCTAHEDRAL] = 'O'
CS_SYM_NAME[SYM_I222] = 'I1'
CS_SYM_NAME[SYM_I222r] = 'I2'

NOISE_MODEL_CHOICES = ['symmetric',
                       'white',
                       'coloured']
REFINE_MASK_CHOICES = ['dynamic',
                       'static',
                       'null']

REFINE_FILTER_TYPE = ['butterworth',
                     'rect',
                     'gaussian']

REFINE_FULCRUM_LOCATION = ['mask_center',
                           'box_center']

COMPUTE_FACILITY_CHOICES = ['GPU',
                            'CPU']

CLASS_3D_INIT_MODE = ['simple', 'PCA', 'input']

EWS_CURVATURE_SIGN = ['positive', 'negative']

EWS_CORRECTION_METHOD = ['simple', 'iterative']

# Viewer constants
LAST_ITER = 0
ALL_ITERS = 1
SELECTED_ITERS = 2

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1


# VOLUME_SLICES = 0
VOLUME_CHIMERA = 1
VOLUME_CRYOSPARC = 0
DATA_VIEWER = 0

fscValues = dict()
fscValues['fsc_nomask'] = 'No mask'
fscValues['fsc_loosemask'] = 'Loose'
fscValues['fsc_tightmask'] = 'Tight'
fscValues['fsc_noisesub_raw'] = 'Noise sub_raw'
fscValues['fsc_noisesub_true'] = 'Noise sub_true'
fscValues['fsc_noisesub'] = 'Corrected'
fscValues['fsc_sphericalmask'] = 'Spherical'
fscValues['fsc_prmm'] = 'Phase randomized'

excludedFSCValues = ['fsc_noisesub_raw', 'fsc_noisesub_true']

HALF_EVEN = 0
HALF_ODD = 1
FULL_MAP = 2
ALL_MAPS = 3

OBJCMD_CLASSAVG_PROJS = 'Show class-averages/projections'
OBJCMD_PROJS = 'Show only projections'
OBJCMD_INITVOL = 'Show initial random volume'


# METADATA

class RELIONCOLUMNS(enum.Enum):
    rlnOriginX = 'rlnOriginX'                        # RLN_ORIENT_ORIGIN_X
    rlnOriginY = 'rlnOriginY'                        # RLN_ORIENT_ORIGIN_Y
    rlnAngleRot = 'rlnAngleRot'                      # RLN_ORIENT_ROT
    rlnAnglePsi = 'rlnAnglePsi'                      # RLN_ORIENT_PSI
    rlnOriginZ = 'rlnOriginZ'                        # RLN_ORIENT_ORIGIN_Z
    rlnClassNumber = 'rlnClassNumber'                # RLN_PARTICLE_CLASS
    rlnImageName = 'rlnImageName'                    # RLN_IMAGE_NAME
    rlnRandomSubset = 'rlnRandomSubset'              # RLN_PARTICLE_RANDOM_SUBSET
    rlnAngleTilt = 'rlnAngleTilt'                    # RLN_ORIENT_TILT
    rlnDefocusU = 'rlnDefocusU'                      # RLN_CTF_DEFOCUSU
    rlnDefocusV = 'rlnDefocusV'                      # RLN_CTF_DEFOCUSV
    rlnDefocusAngle = 'rlnDefocusAngle'              # RLN_CTF_DEFOCUS_ANGLE
    rlnPhaseShift = 'rlnPhaseShift'                  # RLN_CTF_PHASESHIFT
    rlnBfactor = 'rlnBfactor'                        # RLN_CTF_BFACTOR
    rlnOpticsGroup = 'rlnOpticsGroup'
    rlnVoltage = 'rlnVoltage'                        # RLN_CTF_VOLTAGE
    rlnSphericalAberration = 'rlnSphericalAberration'  # RLN_CTF_CS
    rlnAmplitudeContrast = 'rlnAmplitudeContrast'    # RLN_CTF_Q0
    rlnImageSize = 'rlnImageSize'                    # RLN_IMAGE_SIZE
    rlnEnabled = 'rlnEnabled'                        # RLN_IMAGE_ENABLED
    rlnCtfFigureOfMerit = 'rlnCtfFigureOfMerit'     # RLN_CTF_FOM
    rlnMagnification = 'rlnMagnification'           # RLN_CTF_MAGNIFICATION
    rlnDetectorPixelSize = 'rlnDetectorPixelSize'    # RLN_CTF_DETECTOR_PIXEL_SIZE
    rlnCtfImage = 'rlnCtfImage'
    rlnAreaId = 'rlnAreaId'                        # RLN_AREA_ID
    rlnAreaName = 'rlnAreaName'                    # RLN_AREA_NAME
    rlnCtfScalefactor = 'rlnCtfScalefactor'         # RLN_CTF_SCALEFACTOR\
    rlnChromaticAberration = 'rlnChromaticAberration'  # RLN_CTF_CA
    rlnEnergyLoss = 'rlnEnergyLoss'                   # RLN_CTF_ENERGY_LOSS
    rlnLensStability = 'rlnLensStability'             # RLN_CTF_LENS_STABILITY
    rlnCtfMaxResolution = 'rlnCtfMaxResolution'       # RLN_CTF_MAXRES
    rlnConvergenceCone = 'rlnConvergenceCone'         # RLN_CTF_CONVERGENCE_CONE
    rlnLongitudinalDisplacement = 'rlnLongitudinalDisplacement'  # RLN_CTF_LONGITUDINAL_DISPLACEMENT
    rlnTransversalDisplacement = 'rlnTransversalDisplacement'  # RLN_CTF_TRANSVERSAL_DISPLACEMENT
    rlnCtfValidationScore = 'rlnCtfValidationScore'       # RLN_CTF_VALIDATIONSCORE
    rlnCtfValue = 'rlnCtfValue'                      # RLN_CTF_VALUE
    rlnReconstructImageName = 'rlnReconstructImageName'  # RLN_IMAGE_RECONSTRUCT_NAME
    rlnImageId = 'rlnImageId'           # RLN_IMAGE_ID
    rlnDataType = 'rlnDataType'         # RLN_IMAGE_DATATYPE
    rlnAutopickFigureOfMerit = 'rlnAutopickFigureOfMerit'   # RLN_PARTICLE_AUTOPICK_FOM
    rlnCoordinateX = 'rlnCoordinateX'     # RLN_IMAGE_COORD_X
    rlnCoordinateY = 'rlnCoordinateY'     # RLN_IMAGE_COORD_Y
    rlnCoordinateZ = 'rlnCoordinateZ'     # RLN_IMAGE_COORD_Z
    rlnMicrographName = 'rlnMicrographName' # RLN_MICROGRAPH_NAME
    rlnParticleSelectZScore = 'rlnParticleSelectZScore'  # RLN_SELECT_PARTICLES_ZSCORE
    rlnMovieFrameNumber = 'rlnMovieFrameNumber'  # RLN_IMAGE_FRAME_NR
    rlnReferenceImage = 'rlnReferenceImage'     # RLN_MLMODEL_REF_IMAGE
    rlnMicrographId = 'rlnMicrographId'          # RLN_MICROGRAPH_ID
    rlnParticleId = 'rlnParticleId'      # RLN_PARTICLE_ID



"""
        // Label rlnImageDimensionality is originally rlnDataDimensionality, which is
        // duplicated for other label. A relion bug???
        MDL::addLabel(RLN_IMAGE_DIMENSIONALITY, LABEL_INT, "rlnImageDimensionality");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_X, LABEL_DOUBLE, "rlnBeamTiltX");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_Y, LABEL_DOUBLE, "rlnBeamTiltY");
        MDL::addLabel(RLN_IMAGE_BEAMTILT_GROUP, LABEL_STRING, "rlnBeamTiltGroupName");
        MDL::addLabel(RLN_IMAGE_MAGNIFICATION_CORRECTION, LABEL_DOUBLE, "rlnMagnificationCorrection");
        MDL::addLabel(RLN_IMAGE_NORM_CORRECTION, LABEL_DOUBLE, "rlnNormCorrection");
        MDL::addLabel(RLN_IMAGE_ORI_NAME, LABEL_STRING, "rlnImageOriginalName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE, LABEL_DOUBLE, "rlnSamplingRate");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_X, LABEL_DOUBLE, "rlnSamplingRateX");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_Y, LABEL_DOUBLE, "rlnSamplingRateY");
        MDL::addLabel(RLN_IMAGE_SAMPLINGRATE_Z, LABEL_DOUBLE, "rlnSamplingRateZ");
        MDL::addLabel(RLN_IMAGE_SIZEX, LABEL_INT, "rlnImageSizeX");
        MDL::addLabel(RLN_IMAGE_SIZEY, LABEL_INT, "rlnImageSizeY");
        MDL::addLabel(RLN_IMAGE_SIZEZ, LABEL_INT, "rlnImageSizeZ");
        MDL::addLabel(RLN_IMAGE_STATS_MIN, LABEL_DOUBLE, "rlnMinimumValue");
        MDL::addLabel(RLN_IMAGE_STATS_MAX, LABEL_DOUBLE, "rlnMaximumValue");
        MDL::addLabel(RLN_IMAGE_STATS_AVG, LABEL_DOUBLE, "rlnAverageValue");
        MDL::addLabel(RLN_IMAGE_STATS_STDDEV, LABEL_DOUBLE, "rlnStandardDeviationValue");
        MDL::addLabel(RLN_IMAGE_STATS_SKEW, LABEL_DOUBLE, "rlnSkewnessValue");
        MDL::addLabel(RLN_IMAGE_STATS_KURT, LABEL_DOUBLE, "rlnKurtosisExcessValue");
        MDL::addLabel(RLN_IMAGE_WEIGHT, LABEL_DOUBLE, "rlnImageWeight");

        MDL::addLabel(RLN_MASK_NAME, LABEL_STRING, "rlnMaskName");

        MDL::addLabel(RLN_MATRIX_1_1, LABEL_DOUBLE, "rlnMatrix_1_1");
        MDL::addLabel(RLN_MATRIX_1_2, LABEL_DOUBLE, "rlnMatrix_1_2");
        MDL::addLabel(RLN_MATRIX_1_3, LABEL_DOUBLE, "rlnMatrix_1_3");
        MDL::addLabel(RLN_MATRIX_2_1, LABEL_DOUBLE, "rlnMatrix_2_1");
        MDL::addLabel(RLN_MATRIX_2_2, LABEL_DOUBLE, "rlnMatrix_2_2");
        MDL::addLabel(RLN_MATRIX_2_3, LABEL_DOUBLE, "rlnMatrix_2_3");
        MDL::addLabel(RLN_MATRIX_3_1, LABEL_DOUBLE, "rlnMatrix_3_1");
        MDL::addLabel(RLN_MATRIX_3_2, LABEL_DOUBLE, "rlnMatrix_3_2");
        MDL::addLabel(RLN_MATRIX_3_3, LABEL_DOUBLE, "rlnMatrix_3_3");

        MDL::addLabel(RLN_MICROGRAPH_MOVIE_NAME, LABEL_STRING, "rlnMicrographMovieName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MICROGRAPH_NAME_WODOSE, LABEL_STRING, "rlnMicrographNameNoDW", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_MICROGRAPH_TILT_ANGLE, LABEL_DOUBLE, "rlnMicrographTiltAngle");
        MDL::addLabel(RLN_MICROGRAPH_TILT_AXIS_DIRECTION, LABEL_DOUBLE, "rlnMicrographTiltAxisDirection");
        MDL::addLabel(RLN_MICROGRAPH_TILT_AXIS_OUTOFPLANE, LABEL_DOUBLE, "rlnMicrographTiltAxisOutOfPlane");

        MDL::addLabel(RLN_MLMODEL_ACCURACY_ROT, LABEL_DOUBLE, "rlnAccuracyRotations");
        MDL::addLabel(RLN_MLMODEL_ACCURACY_TRANS, LABEL_DOUBLE, "rlnAccuracyTranslations");
        MDL::addLabel(RLN_MLMODEL_AVE_PMAX, LABEL_DOUBLE, "rlnAveragePmax");
        MDL::addLabel(RLN_MLMODEL_CURRENT_RESOLUTION, LABEL_DOUBLE, "rlnCurrentResolution");
        MDL::addLabel(RLN_MLMODEL_CURRENT_SIZE, LABEL_INT, "rlnCurrentImageSize");
        MDL::addLabel(RLN_MLMODEL_DATA_VS_PRIOR_REF, LABEL_DOUBLE, "rlnSsnrMap");
        MDL::addLabel(RLN_MLMODEL_DIMENSIONALITY, LABEL_INT, "rlnReferenceDimensionality");
        MDL::addLabel(RLN_MLMODEL_DIMENSIONALITY_DATA, LABEL_INT, "rlnDataDimensionality");
        MDL::addLabel(RLN_MLMODEL_DIFF2_HALVES_REF, LABEL_DOUBLE, "rlnDiff2RandomHalves");
        MDL::addLabel(RLN_MLMODEL_ESTIM_RESOL_REF, LABEL_DOUBLE, "rlnEstimatedResolution");
        MDL::addLabel(RLN_MLMODEL_FOURIER_COVERAGE_REF, LABEL_DOUBLE, "rlnFourierCompleteness");
        MDL::addLabel(RLN_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, LABEL_DOUBLE, "rlnOverallFourierCompleteness");
        MDL::addLabel(RLN_MLMODEL_FSC_HALVES_REF, LABEL_DOUBLE, "rlnGoldStandardFsc");
        MDL::addLabel(RLN_MLMODEL_GROUP_NAME, LABEL_STRING, "rlnGroupName");
        MDL::addLabel(RLN_MLMODEL_GROUP_NO, LABEL_SIZET, "rlnGroupNumber");
        MDL::addLabel(RLN_MLMODEL_GROUP_NR_PARTICLES, LABEL_SIZET, "rlnGroupNrParticles");
        MDL::addLabel(RLN_MLMODEL_GROUP_SCALE_CORRECTION, LABEL_DOUBLE, "rlnGroupScaleCorrection");
        MDL::addLabel(RLN_MLMODEL_HELICAL_NR_ASU, LABEL_INT, "rlnNrHelicalAsymUnits");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST, LABEL_DOUBLE, "rlnHelicalTwist");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST_MIN, LABEL_DOUBLE, "rlnHelicalTwistMin");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST_MAX, LABEL_DOUBLE, "rlnHelicalTwistMax");
        MDL::addLabel(RLN_MLMODEL_HELICAL_TWIST_INITIAL_STEP, LABEL_DOUBLE, "rlnHelicalTwistInitialStep");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE, LABEL_DOUBLE, "rlnHelicalRise");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE_MIN, LABEL_DOUBLE, "rlnHelicalRiseMin");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE_MAX, LABEL_DOUBLE, "rlnHelicalRiseMax");
        MDL::addLabel(RLN_MLMODEL_HELICAL_RISE_INITIAL_STEP, LABEL_DOUBLE, "rlnHelicalRiseInitialStep");
        MDL::addLabel(RLN_MLMODEL_INTERPOLATOR, LABEL_INT, "rlnFourierSpaceInterpolator");
        MDL::addLabel(RLN_MLMODEL_IS_HELIX, LABEL_BOOL, "rlnIsHelix");
        MDL::addLabel(RLN_MLMODEL_LL, LABEL_DOUBLE, "rlnLogLikelihood");
        MDL::addLabel(RLN_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,  LABEL_INT, "rlnMinRadiusNnInterpolation");
        MDL::addLabel(RLN_MLMODEL_NORM_CORRECTION_AVG, LABEL_DOUBLE, "rlnNormCorrectionAverage");
        MDL::addLabel(RLN_MLMODEL_NR_BODIES, LABEL_INT, "rlnNrBodies");
        MDL::addLabel(RLN_MLMODEL_NR_CLASSES, LABEL_INT, "rlnNrClasses");
        MDL::addLabel(RLN_MLMODEL_NR_GROUPS, LABEL_INT, "rlnNrGroups");
        MDL::addLabel(RLN_MLMODEL_ORIENTABILITY_CONTRIBUTION, LABEL_DOUBLE, "rlnSpectralOrientabilityContribution");
        MDL::addLabel(RLN_MLMODEL_ORIGINAL_SIZE, LABEL_INT, "rlnOriginalImageSize");
        MDL::addLabel(RLN_MLMODEL_PADDING_FACTOR, LABEL_DOUBLE, "rlnPaddingFactor");
        MDL::addLabel(RLN_MLMODEL_PDF_CLASS, LABEL_DOUBLE, "rlnClassDistribution");
        MDL::addLabel(RLN_MLMODEL_PRIOR_OFFX_CLASS, LABEL_DOUBLE, "rlnClassPriorOffsetX");
        MDL::addLabel(RLN_MLMODEL_PRIOR_OFFY_CLASS, LABEL_DOUBLE, "rlnClassPriorOffsetY");
        MDL::addLabel(RLN_MLMODEL_PDF_ORIENT, LABEL_DOUBLE, "rlnOrientationDistribution");
        MDL::addLabel(RLN_MLMODEL_PIXEL_SIZE, LABEL_DOUBLE, "rlnPixelSize");
        MDL::addLabel(RLN_MLMODEL_POWER_REF, LABEL_DOUBLE, "rlnReferenceSpectralPower");
        MDL::addLabel(RLN_MLMODEL_PRIOR_MODE, LABEL_INT, "rlnOrientationalPriorMode");
        MDL::addLabel(RLN_MLMODEL_SGD_GRADIENT_IMAGE, LABEL_STRING, "rlnSGDGradientImage");
        MDL::addLabel(RLN_MLMODEL_SIGMA_OFFSET, LABEL_DOUBLE, "rlnSigmaOffsets");
        MDL::addLabel(RLN_MLMODEL_SIGMA2_NOISE, LABEL_DOUBLE, "rlnSigma2Noise");
        MDL::addLabel(RLN_MLMODEL_SIGMA2_REF, LABEL_DOUBLE, "rlnReferenceSigma2");
        MDL::addLabel(RLN_MLMODEL_SIGMA_ROT, LABEL_DOUBLE, "rlnSigmaPriorRotAngle");
        MDL::addLabel(RLN_MLMODEL_SIGMA_TILT, LABEL_DOUBLE, "rlnSigmaPriorTiltAngle");
        MDL::addLabel(RLN_MLMODEL_SIGMA_PSI, LABEL_DOUBLE, "rlnSigmaPriorPsiAngle");
        MDL::addLabel(RLN_MLMODEL_SSNR_REF, LABEL_DOUBLE, "rlnSignalToNoiseRatio");
        MDL::addLabel(RLN_MLMODEL_TAU2_FUDGE_FACTOR, LABEL_DOUBLE, "rlnTau2FudgeFactor");
        MDL::addLabel(RLN_MLMODEL_TAU2_REF, LABEL_DOUBLE, "rlnReferenceTau2");

        MDL::addLabel(RLN_OPTIMISER_ACCURACY_ROT, LABEL_DOUBLE, "rlnOverallAccuracyRotations");
        MDL::addLabel(RLN_OPTIMISER_ACCURACY_TRANS, LABEL_DOUBLE, "rlnOverallAccuracyTranslations");
        MDL::addLabel(RLN_OPTIMISER_ADAPTIVE_FRACTION, LABEL_DOUBLE, "rlnAdaptiveOversampleFraction");
        MDL::addLabel(RLN_OPTIMISER_ADAPTIVE_OVERSAMPLING, LABEL_INT, "rlnAdaptiveOversampleOrder");
        MDL::addLabel(RLN_OPTIMISER_AUTO_LOCAL_HP_ORDER, LABEL_INT, "rlnAutoLocalSearchesHealpixOrder");
        MDL::addLabel(RLN_OPTIMISER_AVAILABLE_MEMORY, LABEL_DOUBLE, "rlnAvailableMemory");
        MDL::addLabel(RLN_OPTIMISER_BEST_RESOL_THUS_FAR, LABEL_DOUBLE, "rlnBestResolutionThusFar");
        MDL::addLabel(RLN_OPTIMISER_COARSE_SIZE, LABEL_INT, "rlnCoarseImageSize");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, LABEL_DOUBLE, "rlnChangesOptimalOffsets");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, LABEL_DOUBLE, "rlnChangesOptimalOrientations");
        MDL::addLabel(RLN_OPTIMISER_CHANGES_OPTIMAL_CLASSES, LABEL_DOUBLE, "rlnChangesOptimalClasses");
        MDL::addLabel(RLN_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, LABEL_BOOL, "rlnCtfDataArePhaseFlipped");
        MDL::addLabel(RLN_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, LABEL_BOOL, "rlnCtfDataAreCtfPremultiplied");
        MDL::addLabel(RLN_OPTIMISER_DATA_STARFILE, LABEL_STRING, "rlnExperimentalDataStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_DO_AUTO_REFINE, LABEL_BOOL, "rlnDoAutoRefine");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_CTF, LABEL_BOOL, "rlnDoCorrectCtf");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_MAGNIFICATION, LABEL_BOOL, "rlnDoCorrectMagnification");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_NORM, LABEL_BOOL, "rlnDoCorrectNorm");
        MDL::addLabel(RLN_OPTIMISER_DO_CORRECT_SCALE, LABEL_BOOL, "rlnDoCorrectScale");
        MDL::addLabel(RLN_OPTIMISER_DO_HELICAL_REFINE, LABEL_BOOL, "rlnDoHelicalRefine");
        MDL::addLabel(RLN_OPTIMISER_DO_MAP, LABEL_BOOL, "rlnDoMapEstimation");
        MDL::addLabel(RLN_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, LABEL_BOOL, "rlnDoOnlyFlipCtfPhases");
        MDL::addLabel(RLN_OPTIMISER_DO_REALIGN_MOVIES, LABEL_BOOL, "rlnDoRealignMovies");
        MDL::addLabel(RLN_OPTIMISER_DO_SGD, LABEL_BOOL, "rlnDoStochasticGradientDescent");
        MDL::addLabel(RLN_OPTIMISER_DO_SOLVENT_FLATTEN, LABEL_BOOL, "rlnDoSolventFlattening");
        MDL::addLabel(RLN_OPTIMISER_DO_SKIP_ALIGN, LABEL_BOOL, "rlnDoSkipAlign");
        MDL::addLabel(RLN_OPTIMISER_DO_SKIP_ROTATE, LABEL_BOOL, "rlnDoSkipRotate");
        MDL::addLabel(RLN_OPTIMISER_DO_SPLIT_RANDOM_HALVES, LABEL_BOOL, "rlnDoSplitRandomHalves");
        MDL::addLabel(RLN_OPTIMISER_DO_ZERO_MASK, LABEL_BOOL, "rlnDoZeroMask");
        MDL::addLabel(RLN_OPTIMISER_FIX_SIGMA_NOISE, LABEL_BOOL, "rlnFixSigmaNoiseEstimates");
        MDL::addLabel(RLN_OPTIMISER_FIX_SIGMA_OFFSET ,LABEL_BOOL, "rlnFixSigmaOffsetEstimates");
        MDL::addLabel(RLN_OPTIMISER_FIX_TAU, LABEL_BOOL, "rlnFixTauEstimates");
        MDL::addLabel(RLN_OPTIMISER_HAS_CONVERGED, LABEL_BOOL, "rlnHasConverged");
        MDL::addLabel(RLN_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, LABEL_BOOL, "rlnHasHighFscAtResolLimit");
        MDL::addLabel(RLN_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, LABEL_INT, "rlnHasLargeSizeIncreaseIterationsAgo");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_TWIST_INITIAL, LABEL_DOUBLE, "rlnHelicalTwistInitial");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_RISE_INITIAL, LABEL_DOUBLE, "rlnHelicalRiseInitial");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_Z_PERCENTAGE, LABEL_DOUBLE, "rlnHelicalCentralProportion");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, LABEL_DOUBLE, "rlnHelicalMaskTubeInnerDiameter");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, LABEL_DOUBLE, "rlnHelicalMaskTubeOuterDiameter");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, LABEL_BOOL, "rlnHelicalSymmetryLocalRefinement");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_SIGMA_DISTANCE, LABEL_DOUBLE, "rlnHelicalSigmaDistance");
        MDL::addLabel(RLN_OPTIMISER_IGNORE_HELICAL_SYMMETRY, LABEL_BOOL, "rlnIgnoreHelicalSymmetry");
        MDL::addLabel(RLN_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, LABEL_BOOL, "rlnHelicalKeepTiltPriorFixed");
        MDL::addLabel(RLN_OPTIMISER_HIGHRES_LIMIT_EXP, LABEL_DOUBLE, "rlnHighresLimitExpectation");
        MDL::addLabel(RLN_OPTIMISER_HIGHRES_LIMIT_SGD, LABEL_DOUBLE, "rlnHighresLimitSGD");
        MDL::addLabel(RLN_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, LABEL_BOOL, "rlnDoIgnoreCtfUntilFirstPeak");
        MDL::addLabel(RLN_OPTIMISER_INCR_SIZE, LABEL_INT, "rlnIncrementImageSize");
        MDL::addLabel(RLN_OPTIMISER_ITERATION_NO, LABEL_INT, "rlnCurrentIteration");
        MDL::addLabel(RLN_OPTIMISER_LOCAL_SYMMETRY_FILENAME, LABEL_STRING, "rlnLocalSymmetryFile");
        MDL::addLabel(RLN_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, LABEL_DOUBLE, "rlnJoinHalvesUntilThisResolution");
        MDL::addLabel(RLN_OPTIMISER_MAGNIFICATION_RANGE, LABEL_DOUBLE, "rlnMagnificationSearchRange");
        MDL::addLabel(RLN_OPTIMISER_MAGNIFICATION_STEP, LABEL_DOUBLE, "rlnMagnificationSearchStep");
        MDL::addLabel(RLN_OPTIMISER_MAX_COARSE_SIZE, LABEL_INT, "rlnMaximumCoarseImageSize");
        MDL::addLabel(RLN_OPTIMISER_MAX_NR_POOL, LABEL_INT, "rlnMaxNumberOfPooledParticles");
        MDL::addLabel(RLN_OPTIMISER_MODEL_STARFILE, LABEL_STRING, "rlnModelStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_MODEL_STARFILE2, LABEL_STRING, "rlnModelStarFile2", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_NR_ITERATIONS, LABEL_INT, "rlnNumberOfIterations");
        MDL::addLabel(RLN_OPTIMISER_NR_ITER_WO_RESOL_GAIN, LABEL_INT, "rlnNumberOfIterWithoutResolutionGain");
        MDL::addLabel(RLN_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, LABEL_INT, "rlnNumberOfIterWithoutChangingAssignments");
        MDL::addLabel(RLN_OPTIMISER_OUTPUT_ROOTNAME, LABEL_STRING, "rlnOutputRootName");
        MDL::addLabel(RLN_OPTIMISER_PARTICLE_DIAMETER, LABEL_DOUBLE, "rlnParticleDiameter");
        MDL::addLabel(RLN_OPTIMISER_RADIUS_MASK_3D_MAP, LABEL_INT, "rlnRadiusMaskMap");
        MDL::addLabel(RLN_OPTIMISER_RADIUS_MASK_EXP_PARTICLES, LABEL_INT, "rlnRadiusMaskExpImages");
        MDL::addLabel(RLN_OPTIMISER_RANDOM_SEED, LABEL_INT, "rlnRandomSeed");
        MDL::addLabel(RLN_OPTIMISER_REFS_ARE_CTF_CORRECTED, LABEL_BOOL, "rlnRefsAreCtfCorrected");
        MDL::addLabel(RLN_OPTIMISER_SGD_MU, LABEL_DOUBLE, "rlnSgdMuFactor");
        MDL::addLabel(RLN_OPTIMISER_SGD_SIGMA2FUDGE_INI, LABEL_DOUBLE, "rlnSgdSigma2FudgeInitial");
        MDL::addLabel(RLN_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, LABEL_SIZET, "rlnSgdSigma2FudgeHalflife");
        MDL::addLabel(RLN_OPTIMISER_SGD_SUBSET_START, LABEL_INT, "rlnSgdNextSubset");
        MDL::addLabel(RLN_OPTIMISER_SGD_SUBSET_SIZE, LABEL_SIZET, "rlnSgdSubsetSize");
        MDL::addLabel(RLN_OPTIMISER_SGD_WRITE_EVERY_SUBSET, LABEL_INT, "rlnSgdWriteEverySubset");
        MDL::addLabel(RLN_OPTIMISER_SGD_MAX_SUBSETS, LABEL_SIZET, "rlnSgdMaxSubsets");
        MDL::addLabel(RLN_OPTIMISER_SGD_STEPSIZE, LABEL_DOUBLE, "rlnSgdStepsize");


        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, LABEL_INT, "rlnSmallestChangesClasses");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, LABEL_DOUBLE, "rlnSmallestChangesOffsets");
        MDL::addLabel(RLN_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, LABEL_DOUBLE, "rlnSmallestChangesOrientations");
        MDL::addLabel(RLN_OPTIMISER_SAMPLING_STARFILE, LABEL_STRING, "rlnOrientSamplingStarFile", TAGLABEL_METADATA);
        MDL::addLabel(RLN_OPTIMISER_SOLVENT_MASK_NAME, LABEL_STRING, "rlnSolventMaskName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_OPTIMISER_SOLVENT_MASK2_NAME, LABEL_STRING, "rlnSolventMask2Name", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_OPTIMISER_TAU_SPECTRUM_NAME, LABEL_STRING, "rlnTauSpectrumName");
        MDL::addLabel(RLN_OPTIMISER_USE_TOO_COARSE_SAMPLING, LABEL_BOOL, "rlnUseTooCoarseSampling");
        MDL::addLabel(RLN_OPTIMISER_WIDTH_MASK_EDGE, LABEL_INT, "rlnWidthMaskEdge");

        MDL::addLabel(RLN_ORIENT_FLIP, LABEL_BOOL, "rlnIsFlip");
        MDL::addLabel(RLN_ORIENT_ID, LABEL_SIZET, "rlnOrientationsID");
        MDL::addLabel(RLN_ORIENT_ORIGIN_X_PRIOR, LABEL_DOUBLE, "rlnOriginXPrior");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Y_PRIOR, LABEL_DOUBLE, "rlnOriginYPrior");
        MDL::addLabel(RLN_ORIENT_ORIGIN_Z_PRIOR, LABEL_DOUBLE, "rlnOriginZPrior");
        MDL::addLabel(RLN_ORIENT_ROT_PRIOR, LABEL_DOUBLE, "rlnAngleRotPrior");
        MDL::addLabel(RLN_ORIENT_TILT_PRIOR, LABEL_DOUBLE, "rlnAngleTiltPrior");
        MDL::addLabel(RLN_ORIENT_PSI_PRIOR, LABEL_DOUBLE, "rlnAnglePsiPrior");
        MDL::addLabel(RLN_ORIENT_PSI_PRIOR_FLIP_RATIO, LABEL_DOUBLE, "rlnAnglePsiFlipRatio");

        MDL::addLabel(RLN_PARTICLE_AUTOPICK_FOM, LABEL_DOUBLE, "rlnAutopickFigureOfMerit");
        MDL::addLabel(RLN_PARTICLE_DLL, LABEL_DOUBLE, "rlnLogLikeliContribution");
        MDL::addLabel(RLN_PARTICLE_FOM, LABEL_DOUBLE, "rlnParticleFigureOfMerit");
        MDL::addLabel(RLN_PARTICLE_HELICAL_TUBE_ID, LABEL_INT, "rlnHelicalTubeID");
        MDL::addLabel(RLN_PARTICLE_HELICAL_TUBE_PITCH, LABEL_DOUBLE, "rlnHelicalTubePitch");
        MDL::addLabel(RLN_PARTICLE_HELICAL_TRACK_LENGTH, LABEL_DOUBLE, "rlnHelicalTrackLength");
        MDL::addLabel(RLN_PARTICLE_KL_DIVERGENCE, LABEL_DOUBLE, "rlnKullbackLeibnerDivergence");
        MDL::addLabel(RLN_PARTICLE_MOVIE_RUNNING_AVG, LABEL_INT, "rlnMovieFramesRunningAverage");
        MDL::addLabel(RLN_PARTICLE_NAME, LABEL_STRING, "rlnParticleName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_PARTICLE_NR_SIGNIFICANT_SAMPLES, LABEL_INT, "rlnNrOfSignificantSamples"); /**< particle, Number of orientations contributing to weights*/
        MDL::addLabel(RLN_PARTICLE_NR_FRAMES, LABEL_INT, "rlnNrOfFrames");
        MDL::addLabel(RLN_PARTICLE_NR_FRAMES_AVG, LABEL_INT, "rlnAverageNrOfFrames");
        MDL::addLabel(RLN_PARTICLE_ORI_NAME, LABEL_STRING, "rlnOriginalParticleName", TAGLABEL_IMAGE);
        MDL::addLabel(RLN_PARTICLE_PMAX, LABEL_DOUBLE, "rlnMaxValueProbDistribution"); /**< particle, Maximum value of probability distribution */

        MDL::addLabel(RLN_PIPELINE_EDGE_FROM, LABEL_STRING , "rlnPipeLineEdgeFromNode");
        MDL::addLabel(RLN_PIPELINE_EDGE_TO, LABEL_STRING ,"rlnPipeLineEdgeToNode");
        MDL::addLabel(RLN_PIPELINE_EDGE_PROCESS, LABEL_STRING ,"rlnPipeLineEdgeProcess");
        MDL::addLabel(RLN_PIPELINE_JOB_COUNTER, LABEL_INT, "rlnPipeLineJobCounter");
        MDL::addLabel(RLN_PIPELINE_NODE_NAME, LABEL_STRING , "rlnPipeLineNodeName");
        MDL::addLabel(RLN_PIPELINE_NODE_TYPE, LABEL_INT, "rlnPipeLineNodeType");
        MDL::addLabel(RLN_PIPELINE_PROCESS_ALIAS, LABEL_STRING , "rlnPipeLineProcessAlias");
        MDL::addLabel(RLN_PIPELINE_PROCESS_NAME, LABEL_STRING , "rlnPipeLineProcessName");
        MDL::addLabel(RLN_PIPELINE_PROCESS_TYPE, LABEL_INT, "rlnPipeLineProcessType");
        MDL::addLabel(RLN_PIPELINE_PROCESS_STATUS, LABEL_INT, "rlnPipeLineProcessStatus");

        MDL::addLabel(RLN_POSTPROCESS_AMPLCORR_MASKED, LABEL_DOUBLE, "rlnAmplitudeCorrelationMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_AMPLCORR_UNMASKED,  LABEL_DOUBLE, "rlnAmplitudeCorrelationUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_BFACTOR, LABEL_DOUBLE, "rlnBfactorUsedForSharpening");
        MDL::addLabel(RLN_POSTPROCESS_DPR_MASKED, LABEL_DOUBLE, "rlnDifferentialPhaseResidualMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_DPR_UNMASKED,  LABEL_DOUBLE, "rlnDifferentialPhaseResidualUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FINAL_RESOLUTION, LABEL_DOUBLE, "rlnFinalResolution");
        MDL::addLabel(RLN_POSTPROCESS_FSC_GENERAL, LABEL_DOUBLE, "rlnFourierShellCorrelation");
        MDL::addLabel(RLN_POSTPROCESS_FSC_TRUE, LABEL_DOUBLE, "rlnFourierShellCorrelationCorrected");
        MDL::addLabel(RLN_POSTPROCESS_FSC_MASKED, LABEL_DOUBLE, "rlnFourierShellCorrelationMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FSC_UNMASKED, LABEL_DOUBLE, "rlnFourierShellCorrelationUnmaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_FSC_RANDOM_MASKED, LABEL_DOUBLE, "rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_INTERCEPT, LABEL_DOUBLE, "rlnFittedInterceptGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_SLOPE, LABEL_DOUBLE, "rlnFittedSlopeGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_FIT_CORRELATION, LABEL_DOUBLE, "rlnCorrelationFitGuinierPlot");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_IN, LABEL_DOUBLE, "rlnLogAmplitudesOriginal");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_INVMTF, LABEL_DOUBLE, "rlnLogAmplitudesMTFCorrected");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_WEIGHTED, LABEL_DOUBLE, "rlnLogAmplitudesWeighted");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_SHARPENED, LABEL_DOUBLE, "rlnLogAmplitudesSharpened");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_VALUE_INTERCEPT, LABEL_DOUBLE, "rlnLogAmplitudesIntercept");
        MDL::addLabel(RLN_POSTPROCESS_GUINIER_RESOL_SQUARED, LABEL_DOUBLE, "rlnResolutionSquared");
        MDL::addLabel(RLN_POSTPROCESS_MTF_VALUE, LABEL_DOUBLE, "rlnMtfValue");

        MDL::addLabel(RLN_SAMPLING_IS_3D, LABEL_BOOL, "rlnIs3DSampling");
        MDL::addLabel(RLN_SAMPLING_IS_3D_TRANS, LABEL_BOOL, "rlnIs3DTranslationalSampling");
        MDL::addLabel(RLN_SAMPLING_HEALPIX_ORDER, LABEL_INT, "rlnHealpixOrder");
        MDL::addLabel(RLN_SAMPLING_HELICAL_OFFSET_STEP, LABEL_DOUBLE, "rlnHelicalOffsetStep");
        MDL::addLabel(RLN_SAMPLING_LIMIT_TILT, LABEL_DOUBLE, "rlnTiltAngleLimit");
        MDL::addLabel(RLN_SAMPLING_OFFSET_RANGE, LABEL_DOUBLE, "rlnOffsetRange");
        MDL::addLabel(RLN_SAMPLING_OFFSET_STEP, LABEL_DOUBLE, "rlnOffsetStep");
        MDL::addLabel(RLN_SAMPLING_PERTURB, LABEL_DOUBLE, "rlnSamplingPerturbInstance");
        MDL::addLabel(RLN_SAMPLING_PERTURBATION_FACTOR, LABEL_DOUBLE, "rlnSamplingPerturbFactor");
        MDL::addLabel(RLN_SAMPLING_PSI_STEP, LABEL_DOUBLE, "rlnPsiStep");
        MDL::addLabel(RLN_SAMPLING_SYMMETRY, LABEL_STRING, "rlnSymmetryGroup");

        MDL::addLabel(RLN_SELECTED, LABEL_BOOL, "rlnSelected");
        MDL::addLabel(RLN_SORTED_IDX, LABEL_SIZET, "rlnSortedIndex");
        MDL::addLabel(RLN_STARFILE_MOVIE_PARTICLES, LABEL_STRING, "rlnStarFileMovieParticles");
        MDL::addLabel(RLN_PERFRAME_CUMULATIVE_WEIGHT, LABEL_DOUBLE, "rlnPerFrameCumulativeWeight");
        MDL::addLabel(RLN_PERFRAME_RELATIVE_WEIGHT, LABEL_DOUBLE, "rlnPerFrameRelativeWeight");
        MDL::addLabel(RLN_RESOLUTION, LABEL_DOUBLE, "rlnResolution");
        MDL::addLabel(RLN_RESOLUTION_ANGSTROM, LABEL_DOUBLE, "rlnAngstromResolution");
        MDL::addLabel(RLN_RESOLUTION_INVPIXEL, LABEL_DOUBLE, "rlnResolutionInversePixel");
        MDL::addLabel(RLN_SPECTRAL_IDX, LABEL_INT, "rlnSpectralIndex");

        MDL::addLabelAlias(RLN_CTF_BFACTOR, "rlnCtfBfactor"); //Relion-2.0

"""

COOR_EXTRA_LABELS = [RELIONCOLUMNS.rlnAutopickFigureOfMerit.value,
                     RELIONCOLUMNS.rlnClassNumber.value,
                     RELIONCOLUMNS.rlnAnglePsi.value
]

# COOR_EXTRA_LABELS = [ # Additional autopicking-related metadata
#     md.RLN_PARTICLE_AUTOPICK_FOM,
#     md.RLN_PARTICLE_CLASS,
#     md.RLN_ORIENT_PSI
# ]

COOR_DICT = {"_x": RELIONCOLUMNS.rlnCoordinateX.value,
             "_y": RELIONCOLUMNS.rlnCoordinateY.value}

# Some extra labels
IMAGE_EXTRA_LABELS = [RELIONCOLUMNS.rlnParticleSelectZScore.value,
                      RELIONCOLUMNS.rlnMovieFrameNumber.value]

CTF_DICT = {"_defocusU": RELIONCOLUMNS.rlnDefocusU.value,
            "_defocusV": RELIONCOLUMNS.rlnDefocusV.value,
            "_defocusAngle": RELIONCOLUMNS.rlnDefocusAngle.value}

CTF_PSD_DICT = {"_psdFile": RELIONCOLUMNS.rlnCtfImage.value}

CTF_EXTRA_LABELS = [
    RELIONCOLUMNS.rlnCtfFigureOfMerit.value,
    RELIONCOLUMNS.rlnPhaseShift.value,
    RELIONCOLUMNS.rlnAmplitudeContrast.value,
    RELIONCOLUMNS.rlnVoltage.value,
    RELIONCOLUMNS.rlnMagnification.value,
    RELIONCOLUMNS.rlnDetectorPixelSize.value
]

# This dictionary will be used to map
# between CTFModel properties and Xmipp labels

ACQUISITION_DICT = {
    "_amplitudeContrast": RELIONCOLUMNS.rlnAmplitudeContrast.value,
    "_sphericalAberration": RELIONCOLUMNS.rlnSphericalAberration.value,
    "_voltage": RELIONCOLUMNS.rlnVoltage.value,
    "_magnification": RELIONCOLUMNS.rlnMagnification.value}

ALIGNMENT_DICT = {
    "_rlnOriginX": RELIONCOLUMNS.rlnOriginX.value,
    "_rlnOriginY": RELIONCOLUMNS.rlnOriginY.value,
    "_rlnOriginZ": RELIONCOLUMNS.rlnOriginZ.value,
    "_rlnAngleRot": RELIONCOLUMNS.rlnAngleRot.value,
    "_rlnAngleTilt": RELIONCOLUMNS.rlnAngleTilt.value,
    "_rlnAnglePsi": RELIONCOLUMNS.rlnAnglePsi.value}


