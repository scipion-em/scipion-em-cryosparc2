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
from collections import OrderedDict
import pwem.emlib.metadata as md
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

# Supported versions:
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

COMPUTE_FACILITY_CHOICES = ['GPU',
                            'CPU']

# Viewer constants
LAST_ITER = 0
ALL_ITERS = 1
SELECTED_ITERS = 2

ANGDIST_2DPLOT = 0
ANGDIST_CHIMERA = 1


# VOLUME_SLICES = 0
VOLUME_CHIMERA = 0
VOLUME_CRYOSPARC = 1
DATA_VIEWER = 0

FSC_UNMASK = 0
FSC_SPHERICALMASK = 1
FSC_LOOSEMASK = 2
FSC_TIGHTMASK = 3
FSC_CORRECTEDMASK = 4
FSC_ALL = 5

HALF_EVEN = 0
HALF_ODD = 1
FULL_MAP = 2
ALL_MAPS = 3

OBJCMD_CLASSAVG_PROJS = 'Show class-averages/projections'
OBJCMD_PROJS = 'Show only projections'
OBJCMD_INITVOL = 'Show initial random volume'

COOR_EXTRA_LABELS = [
    # Additional autopicking-related metadata
    md.RLN_PARTICLE_AUTOPICK_FOM,
    md.RLN_PARTICLE_CLASS,
    md.RLN_ORIENT_PSI
]

COOR_DICT = OrderedDict([
    ("_x", md.RLN_IMAGE_COORD_X),
    ("_y", md.RLN_IMAGE_COORD_Y)
])

# Some extra labels
IMAGE_EXTRA_LABELS = [
    md.RLN_SELECT_PARTICLES_ZSCORE,
    md.RLN_IMAGE_FRAME_NR
]


CTF_DICT = OrderedDict([
    ("_defocusU", md.RLN_CTF_DEFOCUSU),
    ("_defocusV", md.RLN_CTF_DEFOCUSV),
    ("_defocusAngle", md.RLN_CTF_DEFOCUS_ANGLE)
])

CTF_PSD_DICT = OrderedDict([
    ("_psdFile", md.RLN_CTF_IMAGE)
])

CTF_EXTRA_LABELS = [
    md.RLN_CTF_FOM,
    md.RLN_CTF_PHASESHIFT,
    # In Relion the ctf also contains acquisition information
    md.RLN_CTF_Q0,
    md.RLN_CTF_CS,
    md.RLN_CTF_VOLTAGE,
    md.RLN_CTF_MAGNIFICATION,
    md.RLN_CTF_DETECTOR_PIXEL_SIZE
]

# This dictionary will be used to map
# between CTFModel properties and Xmipp labels

ACQUISITION_DICT = OrderedDict([
    ("_amplitudeContrast", md.RLN_CTF_Q0),
    ("_sphericalAberration", md.RLN_CTF_CS),
    ("_voltage", md.RLN_CTF_VOLTAGE),
    ("_magnification", md.RLN_CTF_MAGNIFICATION)
])

ALIGNMENT_DICT = OrderedDict([
    ("_rlnOriginX", md.RLN_ORIENT_ORIGIN_X),
    ("_rlnOriginY", md.RLN_ORIENT_ORIGIN_Y),
    ("_rlnOriginZ", md.RLN_ORIENT_ORIGIN_Z),
    ("_rlnAngleRot", md.RLN_ORIENT_ROT),
    ("_rlnAngleTilt", md.RLN_ORIENT_TILT),
    ("_rlnAnglePsi", md.RLN_ORIENT_PSI),
])
