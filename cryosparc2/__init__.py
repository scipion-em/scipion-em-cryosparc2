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
import pyworkflow.utils as pwutils
from pyworkflow import VarTypes

from .constants import *

__version__ = '4.2.1'
_references = ['Punjani2017', 'Brubaker2017', 'daniel_asarnow_2019_3576630']
_logo = 'cryosparc2_logo.png'


class Plugin(em.Plugin):
    _url = "https://github.com/scipion-em/scipion-em-cryosparc2"
    _homeVar = CRYOSPARC_HOME
    _pathVars = [CRYOSPARC_HOME]
    _supportedVersions = SUPORTED_VERSIONS

    @classmethod
    def _defineVariables(cls):
        cls._defineVar(CRYOSPARC_HOME, os.environ.get(CRYOSPARC_DIR, ""),
                       description="The root directory where cryoSPARC is installed. Inside it you may have cryosparc_master folder.",
                       var_type=VarTypes.FOLDER)

        cls._defineVar(CRYO_PROJECTS_DIR, "scipion_projects",
                       description="Folder (available to all workers) where scipion will create cryosparc projects. "
                                   "(default to <CRYOSPARC_HOME>/scipion_projects)",
                       var_type=VarTypes.FOLDER)

        cls._defineVar(CRYOSPARC_PASSWORD, None,
                       description='The password with which cryoSPARC was installed. This is only required for the use '
                                   'of the Flexutils plugin and its connection to the 3D flex training protocol.')
        cls._defineVar(CRYOSPARC_USER, None, description='This is the email with which cryoSPARC was installed.')

    @classmethod
    def getPyemEnvActivation(cls):
        return PYEM_ACTIVATION_CMD

    @classmethod
    def getEnviron(cls):
        """ Set up the environment variables needed to launch cryoSparc. """
        environ = pwutils.Environ(os.environ)
        environ.update({'PATH': cls.getHome()},
                       position=pwutils.Environ.BEGIN)

        return environ

    @classmethod
    def getDependencies(cls):
        """ Return a list of dependencies. Include conda if
            activation command was not found. """
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = ['wget']
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def getUserPassword(cls):
        """Get the user password taking into account the environment variable"""
        return cls.getVar(CRYOSPARC_PASSWORD)

    @classmethod
    def getUser(cls):
        """Get the user email taking into account the environment variable"""
        return cls.getVar(CRYOSPARC_USER)

    @classmethod
    def addPyemPackage(cls, env):
        PYEM_INSTALLED = f"pyem_{PYEM_VERSION}_installed"
        ENV_NAME = getPyemEnvName(PYEM_VERSION)

        installCmd = ["pip uninstall -y pyem && ",
                      cls.getCondaActivationCmd(),
                      f'conda create -y -n {ENV_NAME} python=3.8 -c conda-forge -c anaconda && ',
                      f'conda activate {ENV_NAME} && pip install git+https://github.com/asarnow/pyem.git@0394d5bf6096377ca7cc7b6dd74484f1f40f37a8 && pip install numpy==1.23.5']

        # install pyem
        #installCmd.append('pip install git+https://github.com/asarnow/pyem.git@47cf8f70488500be5988b4db1b6ef7002916e0e0')

        # Flag installation finished
        installCmd.append(f'&& touch {PYEM_INSTALLED}')

        pyem_commands = [(" ".join(installCmd), PYEM_INSTALLED)]

        envPath = os.environ.get('PATH', "")
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage('pyem', version=PYEM_VERSION,
                       tar='void.tgz',
                       commands=pyem_commands,
                       neededProgs=cls.getDependencies(),
                       default=True,
                       createBuildDir=True,
                       buildDir='pyem-%s' % PYEM_VERSION,
                       vars=installEnvVars)

    @classmethod
    def defineBinaries(cls, env):
        cls.addPyemPackage(env)
