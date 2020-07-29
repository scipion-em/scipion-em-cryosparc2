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
import pwem as em
from pyworkflow.utils import Environ

from .constants import *

__version__ = '3.1.0'
_references = ['Punjani2017', 'Brubaker2017', 'daniel_asarnow_2019_3576630']
_logo = 'cryosparc2_logo.png'


class Plugin(em.Plugin):
    _homeVar = CRYOSPARC_HOME
    _pathVars = [CRYOSPARC_HOME]
    _supportedVersions = [V2_5_0, V2_8_0, V2_9_0, V2_11_0, V2_12_0, V2_12_2,
                          V2_12_4, V2_13_0, V2_13_2, V2_14_0, V2_14_2, V2_15_0]

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(CRYOSPARC_HOME, os.environ.get(CRYOSPARC_DIR) or "")
        cls._defineVar(CRYO_PROJECTS_DIR, "scipion_projects")

    @classmethod
    def getEnviron(cls):
        """ Setup the environment variables needed to launch cryoSparc. """
        environ = Environ(os.environ)

        environ.update({
            'PATH': Plugin.getHome(),
            'LD_LIBRARY_PATH': str.join(cls.getHome(), 'cryosparclib')
                               + ":" + cls.getHome(),
        }, position=Environ.BEGIN)

        return environ

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(V2_15_0)

    @classmethod
    def defineBinaries(cls, env):
        pyemLibcmd = 'pip install git+https://github.com/asarnow/pyem.git@d46691bcacae63043346e98cec9ff7b621ca1427'
        env.addPipModule('pyem', version='0.4', pipCmd=pyemLibcmd)
