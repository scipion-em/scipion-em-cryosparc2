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
import pyworkflow.em
from pyworkflow.utils import Environ
from constants import *


_references = ['Punjani2017', 'Brubaker2017']
_logo = 'cryosparc2_logo.png'


class Plugin(pyworkflow.em.Plugin):
    _homeVar = CRYOSPARC_HOME

    @classmethod
    def _defineVariables(cls):
        pass

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
        return cls.getActiveVersion().startswith(V2_5_0)

    @classmethod
    def defineBinaries(cls, env):

        pyemLibcmd = "git+https://github.com/asarnow/pyem"

        env.addPipModule('pandas', version='0.24.2')
        env.addPipModule('pyfftw', version='0.11.1')


pyworkflow.em.Domain.registerPlugin(__name__)
